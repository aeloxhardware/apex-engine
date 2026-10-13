// Voir apex/bird_lsq.hpp. Transcription de scipy (trf_no_bounds, lsmr, common.py).
#include "apex/bird_lsq.hpp"

#include "apex/internal/parallel_pool.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace apex {

// Compteurs du dernier appel (phase 0 du plan Bird C++/Vulkan) : temps en
// produits matrice-vecteur, residus, jacobienne, iterations LSMR. Mesure seule,
// aucun effet sur le calcul.
static LsqProfile g_last_profile;

LsqProfile bird_lsq_last_profile() { return g_last_profile; }

// Parallelisme (2026-10-13) : produits J x et J^T u, residus et jacobienne repartis
// sur un pool persistant. Resultat IDENTIQUE AU BIT a 1 fil : chaque ligne de J x
// garde sa somme dans l'ordre des colonnes, chaque composante de J^T u est sommee
// dans l'ordre CROISSANT des lignes (transposee construite une fois), residus et
// jacobienne sont independants par tetraedre / paire. Les reductions (produits
// scalaires, normes) restent sequentielles.
static std::atomic<int> g_lsq_threads{0};   // 0 = tous les coeurs, 1 = sequentiel
static std::atomic<std::int64_t> g_lsq_min_parallel{4096};   // en dessous : un seul bloc

void bird_lsq_set_threads(int threads, std::int64_t min_parallel) {
    g_lsq_threads = threads < 0 ? 0 : threads;
    g_lsq_min_parallel = min_parallel < 1 ? 1 : min_parallel;
}

namespace {

internal::ParallelPool& lsq_pool() {
    static internal::ParallelPool* pool =
        new internal::ParallelPool(std::max(1u, std::thread::hardware_concurrency()));
    return *pool;
}

// fn(begin, end) sur [0, count) : pool si autorise, sinon en un bloc.
void parallel_for(std::int64_t count, const std::function<void(std::int64_t, std::int64_t)>& fn) {
    if (g_lsq_threads.load() == 1) {
        fn(0, count);
        return;
    }
    lsq_pool().run(count, fn, g_lsq_min_parallel.load());
}

using Vec = std::vector<double>;
using Clock = std::chrono::steady_clock;

thread_local LsqProfile g_prof;

double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

double dot(const Vec& a, const Vec& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double norm2(const Vec& a) { return std::sqrt(dot(a, a)); }

double sign(double v) { return v > 0 ? 1.0 : (v < 0 ? -1.0 : 0.0); }

// Jacobienne CSR a structure fixe, valeurs en float32 (comme `.astype(np.float32)`).
struct Csr {
    std::int64_t m = 0, n = 0;
    std::vector<std::int64_t> indptr, indices;
    std::vector<float> values;
    // Transposee (structure seulement) : pour chaque colonne j, positions k dans
    // `values` des coefficients de cette colonne, par ligne CROISSANTE ; t_rows = leurs lignes.
    std::vector<std::int64_t> t_indptr, t_pos, t_rows;

    void build_transpose() {
        t_indptr.assign(n + 1, 0);
        for (std::int64_t k = 0; k < static_cast<std::int64_t>(indices.size()); ++k) ++t_indptr[indices[k] + 1];
        for (std::int64_t j = 0; j < n; ++j) t_indptr[j + 1] += t_indptr[j];
        t_pos.assign(indices.size(), 0);
        t_rows.assign(indices.size(), 0);
        std::vector<std::int64_t> fill(t_indptr.begin(), t_indptr.end() - 1);
        for (std::int64_t r = 0; r < m; ++r)           // lignes croissantes : ordre conserve par colonne
            for (std::int64_t k = indptr[r]; k < indptr[r + 1]; ++k) {
                const std::int64_t q = fill[indices[k]]++;
                t_pos[q] = k;
                t_rows[q] = r;
            }
    }

    // y = J x
    void matvec(const double* x, double* y) const {
        const auto t0 = Clock::now();
        parallel_for(m, [&](std::int64_t begin, std::int64_t end) {
            for (std::int64_t r = begin; r < end; ++r) {
                double s = 0.0;
                for (std::int64_t k = indptr[r]; k < indptr[r + 1]; ++k) s += double(values[k]) * x[indices[k]];
                y[r] = s;
            }
        });
        g_prof.matvec_s += seconds_since(t0);
        ++g_prof.matvec_calls;
    }
    // y = J^T u : y[j] = somme sur les lignes CROISSANTES de J[r, j] u[r], comme
    // l'accumulation sequentielle d'origine (0 + a1 + a2 + ...).
    void rmatvec(const double* u, double* y) const {
        const auto t0 = Clock::now();
        parallel_for(n, [&](std::int64_t begin, std::int64_t end) {
            for (std::int64_t j = begin; j < end; ++j) {
                double s = 0.0;
                for (std::int64_t q = t_indptr[j]; q < t_indptr[j + 1]; ++q)
                    s += double(values[t_pos[q]]) * u[t_rows[q]];
                y[j] = s;
            }
        });
        g_prof.rmatvec_s += seconds_since(t0);
        ++g_prof.rmatvec_calls;
    }
};

class Problem {
public:
    virtual ~Problem() = default;
    virtual void fun(const Vec& x, Vec& f) const = 0;
    virtual void jac(const Vec& x, Csr& J) const = 0;     // remplit J.values (structure deja posee)
};

// compute_jac_scale : normes de colonnes, accumulees en float32 comme scipy.
void jac_scale(const Csr& J, std::vector<float>& scale_inv, bool has_old) {
    std::vector<float> s(J.n, 0.0f);
    for (std::int64_t r = 0; r < J.m; ++r)
        for (std::int64_t k = J.indptr[r]; k < J.indptr[r + 1]; ++k) {
            const float v = J.values[k];
            s[J.indices[k]] += v * v;
        }
    for (std::int64_t j = 0; j < J.n; ++j) {
        float si = std::sqrt(s[j]);
        if (!has_old) {
            if (si == 0.0f) si = 1.0f;
            scale_inv[j] = si;
        } else {
            scale_inv[j] = std::max(si, scale_inv[j]);
        }
    }
}

void sym_ortho(double a, double b, double& c, double& s, double& r) {
    if (b == 0) { c = sign(a); s = 0; r = std::fabs(a); }
    else if (a == 0) { c = 0; s = sign(b); r = std::fabs(b); }
    else if (std::fabs(b) > std::fabs(a)) {
        const double tau = a / b;
        s = sign(b) / std::sqrt(1 + tau * tau);
        c = s * tau;
        r = b / s;
    } else {
        const double tau = b / a;
        c = sign(a) / std::sqrt(1 + tau * tau);
        s = c * tau;
        r = a / c;
    }
}

// Operateur J_h = J diag(d).
struct ScaledOp {
    const Csr& J;
    const Vec& d;
    mutable Vec tmp;
    ScaledOp(const Csr& J_, const Vec& d_) : J(J_), d(d_), tmp(J_.n) {}
    void matvec(const Vec& x, Vec& y) const {
        for (std::int64_t j = 0; j < J.n; ++j) tmp[j] = x[j] * d[j];
        y.assign(J.m, 0.0);
        J.matvec(tmp.data(), y.data());
    }
    void rmatvec(const Vec& u, Vec& y) const {
        y.assign(J.n, 0.0);
        J.rmatvec(u.data(), y.data());
        for (std::int64_t j = 0; j < J.n; ++j) y[j] *= d[j];
    }
};

// scipy.sparse.linalg.lsmr (x0 = None), solution seulement.
Vec lsmr(const ScaledOp& A, const Vec& b, double damp) {
    const double atol = 1e-6, btol = 1e-6, conlim = 1e8;
    const std::int64_t m = A.J.m, n = A.J.n;
    const std::int64_t maxiter = std::min(m, n);
    Vec u = b;
    const double normb = norm2(b);
    Vec x(n, 0.0);
    double beta = normb;
    Vec v(n, 0.0);
    double alpha = 0.0;
    if (beta > 0) {
        for (auto& e : u) e *= 1.0 / beta;
        A.rmatvec(u, v);
        alpha = norm2(v);
    }
    if (alpha > 0) for (auto& e : v) e *= 1.0 / alpha;
    double zetabar = alpha * beta, alphabar = alpha, rho = 1, rhobar = 1, cbar = 1, sbar = 0;
    Vec h = v, hbar(n, 0.0);
    double betadd = beta, betad = 0, rhodold = 1, tautildeold = 0, thetatilde = 0, zeta = 0, d = 0;
    double normA2 = alpha * alpha, maxrbar = 0, minrbar = 1e100;
    double normA = std::sqrt(normA2), condA = 1, normx = 0;
    const double ctol = conlim > 0 ? 1.0 / conlim : 0.0;
    double normr = beta, normar = alpha * beta;
    if (normar == 0 || normb == 0) return x;
    Vec Av(m), Atu(n);
    std::int64_t itn = 0;
    ++g_prof.lsmr_calls;
    while (itn < maxiter) {
        ++itn;
        ++g_prof.lsmr_iterations;
        A.matvec(v, Av);
        for (std::int64_t i = 0; i < m; ++i) u[i] = u[i] * -alpha + Av[i];
        beta = norm2(u);
        if (beta > 0) {
            for (auto& e : u) e *= 1.0 / beta;
            A.rmatvec(u, Atu);
            for (std::int64_t i = 0; i < n; ++i) v[i] = v[i] * -beta + Atu[i];
            alpha = norm2(v);
            if (alpha > 0) for (auto& e : v) e *= 1.0 / alpha;
        }
        double chat, shat, alphahat;
        sym_ortho(alphabar, damp, chat, shat, alphahat);
        const double rhoold = rho;
        double c, s;
        sym_ortho(alphahat, beta, c, s, rho);
        const double thetanew = s * alpha;
        alphabar = c * alpha;
        const double rhobarold = rhobar, zetaold = zeta;
        const double thetabar = sbar * rho, rhotemp = cbar * rho;
        sym_ortho(cbar * rho, thetanew, cbar, sbar, rhobar);
        zeta = cbar * zetabar;
        zetabar = -sbar * zetabar;
        const double fh = -(thetabar * rho / (rhoold * rhobarold));
        const double fx = zeta / (rho * rhobar);
        for (std::int64_t i = 0; i < n; ++i) {
            hbar[i] = hbar[i] * fh + h[i];
            x[i] += fx * hbar[i];
        }
        const double fhh = -(thetanew / rho);
        for (std::int64_t i = 0; i < n; ++i) h[i] = h[i] * fhh + v[i];
        const double betaacute = chat * betadd, betacheck = -shat * betadd;
        const double betahat = c * betaacute;
        betadd = -s * betaacute;
        const double thetatildeold = thetatilde;
        double ctildeold, stildeold, rhotildeold;
        sym_ortho(rhodold, thetabar, ctildeold, stildeold, rhotildeold);
        thetatilde = stildeold * rhobar;
        rhodold = ctildeold * rhobar;
        betad = -stildeold * betad + ctildeold * betahat;
        tautildeold = (zetaold - thetatildeold * tautildeold) / rhotildeold;
        const double taud = (zeta - thetatilde * tautildeold) / rhodold;
        d = d + betacheck * betacheck;
        normr = std::sqrt(d + (betad - taud) * (betad - taud) + betadd * betadd);
        normA2 = normA2 + beta * beta;
        normA = std::sqrt(normA2);
        normA2 = normA2 + alpha * alpha;
        maxrbar = std::max(maxrbar, rhobarold);
        if (itn > 1) minrbar = std::min(minrbar, rhobarold);
        condA = std::max(maxrbar, rhotemp) / std::min(minrbar, rhotemp);
        normar = std::fabs(zetabar);
        normx = norm2(x);
        const double test1 = normr / normb;
        const double test2 = (normA * normr) != 0 ? normar / (normA * normr) : std::numeric_limits<double>::infinity();
        const double test3 = 1 / condA;
        const double t1 = test1 / (1 + normA * normx / normb);
        const double rtol = btol + atol * normA * normx / normb;
        int istop = 0;
        if (itn >= maxiter) istop = 7;
        if (1 + test3 <= 1) istop = 6;
        if (1 + test2 <= 1) istop = 5;
        if (1 + t1 <= 1) istop = 4;
        if (test3 <= ctol) istop = 3;
        if (test2 <= atol) istop = 2;
        if (test1 <= rtol) istop = 1;
        if (istop > 0) break;
    }
    return x;
}

// QR economique de [a b] (Householder), colonnes orthonormees de Q.
void qr2(const Vec& a, const Vec& b, Vec& q1, Vec& q2) {
    const std::size_t n = a.size();
    Eigen::MatrixXd A(n, 2);
    for (std::size_t i = 0; i < n; ++i) { A(i, 0) = a[i]; A(i, 1) = b[i]; }
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(A);
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(n, 2);
    q1.resize(n); q2.resize(n);
    for (std::size_t i = 0; i < n; ++i) { q1[i] = Q(i, 0); q2[i] = Q(i, 1); }
}

// solve_trust_region_2d (Cholesky, sinon quartique par valeurs propres du compagnon).
void trust_region_2d(const double B[2][2], const double g[2], double Delta, double p[2]) {
    if (B[0][0] > 0) {
        const double l11 = std::sqrt(B[0][0]);
        const double l21 = B[0][1] / l11;
        const double d = B[1][1] - l21 * l21;
        if (d > 0) {
            const double l22 = std::sqrt(d);
            const double y0 = -g[0] / l11;
            const double y1 = (-g[1] - l21 * y0) / l22;
            p[1] = y1 / l22;
            p[0] = (y0 - l21 * p[1]) / l11;
            if (p[0] * p[0] + p[1] * p[1] <= Delta * Delta) return;
        }
    }
    const double a = B[0][0] * Delta * Delta, b = B[0][1] * Delta * Delta, c = B[1][1] * Delta * Delta;
    const double d = g[0] * Delta, f = g[1] * Delta;
    std::vector<double> coeffs = {-b + d, 2 * (a - c + f), 6 * b, 2 * (-a + c + f), -b - d};
    // np.roots : zeros de tete retires, zeros de queue -> racines nulles.
    std::size_t lead = 0;
    while (lead < coeffs.size() && coeffs[lead] == 0.0) ++lead;
    coeffs.erase(coeffs.begin(), coeffs.begin() + static_cast<std::ptrdiff_t>(lead));
    std::size_t trailing = 0;
    while (!coeffs.empty() && coeffs.back() == 0.0) { coeffs.pop_back(); ++trailing; }
    std::vector<double> roots(trailing, 0.0);
    const int deg = static_cast<int>(coeffs.size()) - 1;
    if (deg >= 1) {
        Eigen::MatrixXd C = Eigen::MatrixXd::Zero(deg, deg);
        for (int j = 0; j < deg; ++j) C(0, j) = -coeffs[j + 1] / coeffs[0];
        for (int i = 1; i < deg; ++i) C(i, i - 1) = 1.0;
        Eigen::EigenSolver<Eigen::MatrixXd> es(C, false);
        for (int i = 0; i < deg; ++i)
            if (es.eigenvalues()[i].imag() == 0.0) roots.push_back(es.eigenvalues()[i].real());
    }
    double best = std::numeric_limits<double>::infinity();
    p[0] = p[1] = 0.0;
    for (double t : roots) {
        const double q0 = Delta * 2 * t / (1 + t * t), q1 = Delta * (1 - t * t) / (1 + t * t);
        const double value = 0.5 * (q0 * (B[0][0] * q0 + B[0][1] * q1) + q1 * (B[1][0] * q0 + B[1][1] * q1))
                             + g[0] * q0 + g[1] * q1;
        if (value < best) { best = value; p[0] = q0; p[1] = q1; }
    }
}

// trf_no_bounds, tr_solver="lsmr", regularize=True, damp=0, perte lineaire.
LsqResult trf_impl(const Problem& prob, Csr J, Vec x, double ftol, double xtol, double gtol, int max_nfev,
                   bool jac_scale_on);

// Enveloppe de mesure : compteurs remis a zero, temps total, copie dans g_last_profile.
LsqResult trf(const Problem& prob, Csr J, Vec x, double ftol, double xtol, double gtol, int max_nfev, bool jac_scale_on) {
    g_prof = LsqProfile{};
    g_prof.rows = J.m;
    g_prof.cols = J.n;
    g_prof.nnz = static_cast<std::int64_t>(J.indices.size());
    const auto t0 = Clock::now();
    LsqResult out = trf_impl(prob, std::move(J), std::move(x), ftol, xtol, gtol, max_nfev, jac_scale_on);
    g_prof.total_s = seconds_since(t0);
    g_prof.nfev = out.nfev;
    g_last_profile = g_prof;
    return out;
}

void timed_fun(const Problem& prob, const Vec& x, Vec& f) {
    const auto t0 = Clock::now();
    prob.fun(x, f);
    g_prof.fun_s += seconds_since(t0);
    ++g_prof.fun_calls;
}

void timed_jac(const Problem& prob, const Vec& x, Csr& J) {
    const auto t0 = Clock::now();
    prob.jac(x, J);
    g_prof.jac_s += seconds_since(t0);
    ++g_prof.jac_calls;
}

LsqResult trf_impl(const Problem& prob, Csr J, Vec x, double ftol, double xtol, double gtol, int max_nfev,
                   bool jac_scale_on) {
    const std::int64_t m = J.m, n = J.n;
    Vec f(m), f_new(m), g(n), x_new(n), step(n), step_h(n), g_h(n), Jv(m);
    timed_fun(prob, x, f);
    for (double v : f) if (!std::isfinite(v)) throw std::runtime_error("residus non finis au point initial");
    timed_jac(prob, x, J);
    int nfev = 1;
    double cost = 0.5 * dot(f, f);
    J.rmatvec(f.data(), g.data());
    std::vector<float> scale_inv(n, 1.0f);
    Vec scale(n, 1.0);
    if (jac_scale_on) {
        jac_scale(J, scale_inv, false);
        for (std::int64_t j = 0; j < n; ++j) scale[j] = double(1.0f / scale_inv[j]);
    }
    double Delta = 0.0;
    for (std::int64_t j = 0; j < n; ++j) { const double v = x[j] * double(scale_inv[j]); Delta += v * v; }
    Delta = std::sqrt(Delta);
    if (Delta == 0) Delta = 1.0;
    int status = -1;
    while (true) {
        double g_norm = 0.0;
        for (double v : g) g_norm = std::max(g_norm, std::fabs(v));
        if (g_norm < gtol) status = 1;
        if (status != -1 || nfev == max_nfev) break;

        const Vec& d = scale;
        for (std::int64_t j = 0; j < n; ++j) g_h[j] = d[j] * g[j];
        ScaledOp Jh(J, d);
        // regularisation (build_quadratic_1d / minimize_quadratic_1d)
        Vec minus_g(n);
        for (std::int64_t j = 0; j < n; ++j) minus_g[j] = -g_h[j];
        Jh.matvec(minus_g, Jv);
        const double qa = 0.5 * dot(Jv, Jv), qb = dot(g_h, minus_g);
        const double to_tr = Delta / norm2(g_h);
        double ts[3] = {0.0, to_tr, 0.0};
        int nt = 2;
        if (qa != 0) {
            const double ext = -0.5 * qb / qa;
            if (0 < ext && ext < to_tr) ts[nt++] = ext;
        }
        double ag_value = std::numeric_limits<double>::infinity();
        for (int i = 0; i < nt; ++i) {
            const double y = ts[i] * (qa * ts[i] + qb);
            if (y < ag_value) ag_value = y;
        }
        const double reg_term = -ag_value / (Delta * Delta);
        const double damp_full = std::sqrt(reg_term);
        const Vec gn_h = lsmr(Jh, f, damp_full);
        Vec s1, s2;
        qr2(g_h, gn_h, s1, s2);
        Vec Js1(m), Js2(m);
        Jh.matvec(s1, Js1);
        Jh.matvec(s2, Js2);
        const double B[2][2] = {{dot(Js1, Js1), dot(Js1, Js2)}, {dot(Js2, Js1), dot(Js2, Js2)}};
        const double gS[2] = {dot(s1, g_h), dot(s2, g_h)};

        double actual_reduction = -1, cost_new = cost;
        while (actual_reduction <= 0 && nfev < max_nfev) {
            double pS[2];
            trust_region_2d(B, gS, Delta, pS);
            for (std::int64_t j = 0; j < n; ++j) step_h[j] = s1[j] * pS[0] + s2[j] * pS[1];
            Jh.matvec(step_h, Jv);
            const double predicted_reduction = -(0.5 * dot(Jv, Jv) + dot(step_h, g_h));
            for (std::int64_t j = 0; j < n; ++j) { step[j] = d[j] * step_h[j]; x_new[j] = x[j] + step[j]; }
            timed_fun(prob, x_new, f_new);
            ++nfev;
            const double step_h_norm = norm2(step_h);
            bool finite = true;
            for (double v : f_new) if (!std::isfinite(v)) { finite = false; break; }
            if (!finite) { Delta = 0.25 * step_h_norm; continue; }
            cost_new = 0.5 * dot(f_new, f_new);
            actual_reduction = cost - cost_new;
            // update_tr_radius
            double ratio;
            if (predicted_reduction > 0) ratio = actual_reduction / predicted_reduction;
            else if (predicted_reduction == actual_reduction && actual_reduction == 0) ratio = 1;
            else ratio = 0;
            double Delta_new = Delta;
            if (ratio < 0.25) Delta_new = 0.25 * step_h_norm;
            else if (ratio > 0.75 && step_h_norm > 0.95 * Delta) Delta_new = Delta * 2.0;
            const double step_norm = norm2(step);
            // check_termination
            const bool ftol_ok = actual_reduction < ftol * cost && ratio > 0.25;
            const bool xtol_ok = step_norm < xtol * (xtol + norm2(x));
            if (ftol_ok && xtol_ok) status = 4;
            else if (ftol_ok) status = 2;
            else if (xtol_ok) status = 3;
            if (status != -1) break;
            Delta = Delta_new;
        }
        if (actual_reduction > 0) {
            x = x_new;
            f = f_new;
            cost = cost_new;
            timed_jac(prob, x, J);
            J.rmatvec(f.data(), g.data());
            if (jac_scale_on) {
                jac_scale(J, scale_inv, true);
                for (std::int64_t j = 0; j < n; ++j) scale[j] = double(1.0f / scale_inv[j]);
            }
        }
        if (status != -1) break;
    }
    LsqResult out;
    out.x = std::move(x);
    out.cost = cost;
    out.nfev = nfev;
    out.status = status == -1 ? 0 : status;
    return out;
}

constexpr double kN[4][4] = {{0.75, -0.25, -0.25, -0.25}, {-0.25, 0.75, -0.25, -0.25},
                             {-0.25, -0.25, 0.75, -0.25}, {-0.25, -0.25, -0.25, 0.75}};

class DeformProblem : public Problem {
public:
    DeformProblem(const double* base, std::int64_t n_points, const std::int64_t* cells, std::int64_t n_cells,
                  const std::int64_t* column, const double* old)
        : base_(base, base + n_points * 3), cells_(cells), n_cells_(n_cells), column_(column, column + n_points),
          old_(old), n_points_(n_points) {}

    // Structure : ligne = cellule ; colonnes (column[v]*3 + dim) triees.
    Csr pattern(std::int64_t n_free) {
        Csr J;
        J.m = n_cells_;
        J.n = n_free * 3;
        J.indptr.assign(n_cells_ + 1, 0);
        slot_.clear();
        for (std::int64_t c = 0; c < n_cells_; ++c) {
            std::vector<std::pair<std::int64_t, int>> entries;   // (colonne, i*3+dim)
            for (int i = 0; i < 4; ++i) {
                const std::int64_t col = column_[cells_[c * 4 + i]];
                if (col < 0) continue;
                for (int dim = 0; dim < 3; ++dim) entries.emplace_back(col * 3 + dim, i * 3 + dim);
            }
            std::sort(entries.begin(), entries.end());
            for (const auto& e : entries) { J.indices.push_back(e.first); slot_.push_back(e.second); }
            J.indptr[c + 1] = static_cast<std::int64_t>(J.indices.size());
        }
        J.values.assign(J.indices.size(), 0.0f);
        return J;
    }

    void diff(const Vec& x, std::int64_t c, double out[3][4]) const {
        double V[4][3];
        for (int i = 0; i < 4; ++i) {
            const std::int64_t v = cells_[c * 4 + i];
            const std::int64_t col = column_[v];
            for (int k = 0; k < 3; ++k) V[i][k] = col >= 0 ? x[col * 3 + k] : base_[v * 3 + k];
        }
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 3; ++k) {
                double s = 0.0;
                for (int i = 0; i < 4; ++i) s += kN[j][i] * V[i][k];
                out[k][j] = s - old_[(c * 3 + k) * 4 + j];
            }
    }

    void fun(const Vec& x, Vec& f) const override {
        parallel_for(n_cells_, [&](std::int64_t begin, std::int64_t end) {
            for (std::int64_t c = begin; c < end; ++c) {
                double dd[3][4];
                diff(x, c, dd);
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    for (int j = 0; j < 4; ++j) s += dd[k][j] * dd[k][j];
                const double nrm = std::sqrt(s);
                f[c] = nrm * nrm;                       // np.linalg.norm(...) ** 2
            }
        });
    }

    void jac(const Vec& x, Csr& J) const override {
        parallel_for(n_cells_, [&](std::int64_t begin, std::int64_t end) {
            for (std::int64_t c = begin; c < end; ++c) {
                double dd[3][4];
                diff(x, c, dd);
                for (std::int64_t k = J.indptr[c]; k < J.indptr[c + 1]; ++k) {
                    const int slot = slot_[k];
                    J.values[k] = static_cast<float>(2.0 * dd[slot % 3][slot / 3]);
                }
            }
        });
    }

private:
    Vec base_;
    const std::int64_t* cells_;
    std::int64_t n_cells_;
    std::vector<std::int64_t> column_;
    const double* old_;
    std::int64_t n_points_;
    std::vector<int> slot_;
};

class RotationProblem : public Problem {
public:
    RotationProblem(const std::int64_t* pairs, std::int64_t n_pairs, const std::int64_t* targets,
                    std::int64_t n_targets, const double* initial, std::int64_t n_cells, const std::int64_t* column,
                    double weight)
        : pairs_(pairs), n_pairs_(n_pairs), targets_(targets), n_targets_(n_targets), initial_(initial),
          n_cells_(n_cells), column_(column), w_(weight) {}

    Csr pattern(std::int64_t n_free) {
        Csr J;
        J.m = n_pairs_ + n_targets_;
        J.n = n_free;
        J.indptr.assign(J.m + 1, 0);
        sign_.clear();
        for (std::int64_t r = 0; r < n_pairs_; ++r) {
            std::vector<std::pair<std::int64_t, double>> e;
            const std::int64_t a = column_[pairs_[r * 2]], b = column_[pairs_[r * 2 + 1]];
            if (a >= 0) e.emplace_back(a, 1.0);
            if (b >= 0) e.emplace_back(b, -1.0);
            std::sort(e.begin(), e.end());
            for (const auto& q : e) { J.indices.push_back(q.first); sign_.push_back(q.second); }
            J.indptr[r + 1] = static_cast<std::int64_t>(J.indices.size());
        }
        for (std::int64_t t = 0; t < n_targets_; ++t) {
            J.indices.push_back(column_[targets_[t]]);
            sign_.push_back(0.0);
            J.indptr[n_pairs_ + t + 1] = static_cast<std::int64_t>(J.indices.size());
        }
        J.values.assign(J.indices.size(), 0.0f);
        return J;
    }

    double value(const Vec& x, std::int64_t c) const { return column_[c] >= 0 ? x[column_[c]] : 0.0; }

    void fun(const Vec& x, Vec& f) const override {
        for (std::int64_t r = 0; r < n_pairs_; ++r) {
            const double d = value(x, pairs_[r * 2]) - value(x, pairs_[r * 2 + 1]);
            f[r] = w_ * (d * d);
        }
        for (std::int64_t t = 0; t < n_targets_; ++t) {
            const std::int64_t c = targets_[t];
            const double d = value(x, c) - initial_[c];
            f[n_pairs_ + t] = d * d;
        }
    }

    void jac(const Vec& x, Csr& J) const override {
        for (std::int64_t r = 0; r < n_pairs_; ++r) {
            const double d = value(x, pairs_[r * 2]) - value(x, pairs_[r * 2 + 1]);
            for (std::int64_t k = J.indptr[r]; k < J.indptr[r + 1]; ++k)
                J.values[k] = static_cast<float>(sign_[k] * 2 * w_ * d);
        }
        for (std::int64_t t = 0; t < n_targets_; ++t) {
            const std::int64_t c = targets_[t];
            J.values[J.indptr[n_pairs_ + t]] = static_cast<float>(2 * (value(x, c) - initial_[c]));
        }
    }

private:
    const std::int64_t* pairs_;
    std::int64_t n_pairs_;
    const std::int64_t* targets_;
    std::int64_t n_targets_;
    const double* initial_;
    std::int64_t n_cells_;
    const std::int64_t* column_;
    double w_;
    std::vector<double> sign_;
};

}  // namespace

LsqResult bird_lsq_deform(const double* base_vertices, std::int64_t n_points, const std::int64_t* cells,
                          std::int64_t n_cells, const std::int64_t* column, std::int64_t n_free,
                          const double* old_transformed, int max_nfev, double ftol, double xtol, double gtol) {
    DeformProblem prob(base_vertices, n_points, cells, n_cells, column, old_transformed);
    Csr J = prob.pattern(n_free);
    J.build_transpose();
    Vec x0(n_free * 3);
    for (std::int64_t v = 0; v < n_points; ++v)
        if (column[v] >= 0)
            for (int k = 0; k < 3; ++k) x0[column[v] * 3 + k] = base_vertices[v * 3 + k];
    return trf(prob, std::move(J), std::move(x0), ftol, xtol, gtol, max_nfev, true);
}

LsqResult bird_lsq_rotation(const std::int64_t* pairs, std::int64_t n_pairs, const std::int64_t* targets,
                            std::int64_t n_targets, const double* initial, std::int64_t n_cells,
                            const std::int64_t* column, std::int64_t n_free, double weight, int max_nfev,
                            double ftol, double xtol, double gtol) {
    RotationProblem prob(pairs, n_pairs, targets, n_targets, initial, n_cells, column, weight);
    Csr J = prob.pattern(n_free);
    J.build_transpose();
    return trf(prob, std::move(J), Vec(n_free, 0.0), ftol, xtol, gtol, max_nfev, false);
}

}  // namespace apex

namespace apex {

void viterbi_choice(const double* mapped, const std::uint8_t* valid, std::int64_t n, std::int64_t k,
                    std::int64_t* choice) {
    if (n <= 0) return;
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> cost(k), next(k);
    std::vector<std::int64_t> back(static_cast<std::size_t>(n * k), 0);
    for (std::int64_t j = 0; j < k; ++j) cost[j] = valid[j] ? 0.0 : inf;
    for (std::int64_t i = 1; i < n; ++i) {
        const double* prev = mapped + (i - 1) * k * 3;
        const double* cur = mapped + i * k * 3;
        for (std::int64_t b = 0; b < k; ++b) {
            // total[a, b] = cost[a] + |prev[a] - cur[b]| ; argmin sur a (premier minimum)
            double best = 0.0;
            std::int64_t arg = 0;
            for (std::int64_t a = 0; a < k; ++a) {
                const double dx = prev[a * 3] - cur[b * 3];
                const double dy = prev[a * 3 + 1] - cur[b * 3 + 1];
                const double dz = prev[a * 3 + 2] - cur[b * 3 + 2];
                const double total = cost[a] + std::sqrt(dx * dx + dy * dy + dz * dz);
                if (a == 0 || total < best || (std::isnan(total) && !std::isnan(best))) {
                    best = total;
                    arg = a;
                }
            }
            back[i * k + b] = arg;
            next[b] = valid[i * k + b] ? best : inf;
        }
        cost.swap(next);
    }
    // cost.argmin() : premier minimum, et un NaN l'emporte (regle de numpy).
    std::int64_t last = 0;
    for (std::int64_t j = 1; j < k; ++j) {
        if (std::isnan(cost[last])) break;
        if (std::isnan(cost[j]) || cost[j] < cost[last]) last = j;
    }
    choice[n - 1] = last;
    for (std::int64_t i = n - 1; i > 0; --i) choice[i - 1] = back[i * k + choice[i]];
}

}  // namespace apex
