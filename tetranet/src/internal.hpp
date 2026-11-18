// TetraNet — types et outils internes partagés entre les fichiers de src/ (pas d'API publique ici).
#pragma once
#define _USE_MATH_DEFINES  // M_PI sous MSVC
#include "tetranet/tetranet.hpp"
#include "tetranet/features.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

namespace tetranet::detail {

using Clock = std::chrono::steady_clock;
inline double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

struct V3 {
    double x = 0, y = 0, z = 0;
    V3() = default;
    V3(double a, double b, double c) : x(a), y(b), z(c) {}
    V3 operator+(const V3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V3 operator-(const V3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V3 operator*(double s) const { return {x * s, y * s, z * s}; }
};
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(const V3& a, const V3& b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }
inline V3 unit(const V3& a) { double n = norm(a); return n > 0 ? a * (1.0 / n) : V3(); }
inline double tet_volume(const V3& a, const V3& b, const V3& c, const V3& d) { return dot(cross(b - a, c - a), d - a) / 6.0; }

template <class F>
void parallel_for(int64_t n, unsigned threads, F fn) {  // fn(begin, end, thread)
    if (threads <= 1 || n < 2048) { fn(0, n, 0u); return; }
    std::vector<std::thread> pool;
    int64_t chunk = (n + threads - 1) / threads;
    for (unsigned t = 0; t < threads; ++t) {
        int64_t b = t * chunk, e = std::min<int64_t>(n, b + chunk);
        if (b >= e) break;
        pool.emplace_back([=, &fn] { fn(b, e, t); });
    }
    for (auto& th : pool) th.join();
}

// ---------------------------------------------------------------- progression réelle
// Sérialise les rappels, garde la fraction croissante dans une phase et limite le débit (0,5 %).
class Progress {
public:
    explicit Progress(const std::function<void(int, double)>& fn) : fn_(fn) {}
    void report(int phase, double f) {
        if (!fn_) return;
        std::lock_guard<std::mutex> g(m_);
        f = std::min(1.0, std::max(0.0, f));
        if (phase == phase_) { if (f <= last_ || (f < 1.0 && f - last_ < 0.005)) return; }
        else phase_ = phase;
        last_ = f;
        fn_(phase, f);
    }
private:
    const std::function<void(int, double)>& fn_;
    std::mutex m_;
    int phase_ = -1;
    double last_ = -1.0;
};

// Tranche [lo, hi] d'une phase pour une boucle de `total` tours : les fils appellent `add` par paquets.
class Span {
public:
    Span(Progress& p, int phase, double lo, double hi, int64_t total) : p_(p), phase_(phase), lo_(lo), hi_(hi), total_(std::max<int64_t>(1, total)) {}
    void add(int64_t n) {
        const int64_t d = done_.fetch_add(n) + n;
        p_.report(phase_, lo_ + (hi_ - lo_) * std::min(1.0, double(d) / double(total_)));
    }
    // À appeler en tête de chaque tour : un paquet de BATCH tours terminés est signalé d'un coup.
    struct Tick {
        Span& s; int64_t n = 0;
        void operator()() { if (++n == BATCH) { s.add(n); n = 0; } }
        ~Tick() { if (n) s.add(n); }
    };
    Tick tick() { return Tick{*this}; }
private:
    static constexpr int64_t BATCH = 2048;
    Progress& p_; int phase_; double lo_, hi_; int64_t total_;
    std::atomic<int64_t> done_{0};
};

// ---------------------------------------------------------------- surface + BVH (surface.cpp)

struct Surface {
    std::vector<V3> p;                       // sommets soudés
    std::vector<std::array<int, 3>> tri;     // triangles non dégénérés
    std::vector<V3> fn;                      // normales de face (unitaires)
    std::vector<V3> vn;                      // pseudonormales de sommet
    std::unordered_map<uint64_t, V3> en;     // pseudonormales d'arête
    double volume = 0.0;
    int64_t reoriented = 0, holes = 0, open_edges = 0, nonmanifold_edges = 0, flipped_edges = 0, edges = 0, degenerate = 0, duplicates = 0;
};

inline uint64_t mix64(uint64_t x) {   // finaliseur splitmix64
    x += 0x9E3779B97F4A7C15ull; x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ull; x = (x ^ (x >> 27)) * 0x94D049BB133111EBull;
    return x ^ (x >> 31);
}
inline uint64_t ekey(int a, int b) { if (a > b) std::swap(a, b); return (uint64_t(a) << 32) | uint32_t(b); }

// Soude, retire les faces dégénérées, mesure la topologie et (option) remet les faces dans un sens cohérent.
bool build_surface(const SurfaceMesh& in, Surface& s, bool fix_orientation);

struct Node { double lo[3], hi[3]; int left = -1, right = -1, first = 0, count = 0; };

struct Bvh {
    const Surface* s = nullptr;
    std::vector<Node> nodes;
    std::vector<int> order;

    void build(const Surface& surf);
    int rec(int b, int e);
    static double box2(const Node& n, const V3& q) {
        double d = 0, a[3] = {q.x, q.y, q.z};
        for (int k = 0; k < 3; ++k) { double v = std::max({n.lo[k] - a[k], 0.0, a[k] - n.hi[k]}); d += v * v; }
        return d;
    }
};

// ---------------------------------------------------------------- nombre d'enroulement généralisé
// w(p) = (1/4π) Σ angle solide des triangles vus de p : 1 dedans, 0 dehors, entre les deux près des défauts.
// Les nœuds lointains sont remplacés par un dipôle (Barill et al., 2018) : normale pondérée par l'aire.
struct WData { V3 c, n; double area = 0, r = 0; };

struct Winding {
    const Bvh* bvh = nullptr;
    std::vector<WData> d;
    double beta = 3.0;

    void build(const Bvh& b);
    static double solid(const V3& a, const V3& b, const V3& c);
    double query(const V3& p, double beta_override = 0.0) const;
};

// point le plus proche d'un triangle (Ericson) ; region : 0 face, 1..3 arêtes ab/bc/ca, 4..6 sommets
V3 closest_on_tri(const V3& p, const V3& a, const V3& b, const V3& c, int& region);

struct Hit { bool found = false; double d2 = 0; V3 q; int tri = -1, region = 0; };

Hit nearest(const Bvh& bvh, const V3& p, double maxd);

// +1 dehors, -1 dedans (pseudonormale), 0 si p est sur la surface
int sign_at(const Surface& s, const Hit& h, const V3& p);


// ---------------------------------------------------------------- arêtes vives (piste C)
std::unique_ptr<FeatureEdges> make_features(const Surface& s, const Params& prm, double radius);
V3 snap_to_feature(const FeatureEdges* features, const V3& q, double radius);

}  // namespace tetranet::detail
