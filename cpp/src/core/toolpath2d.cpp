// Voir apex/toolpath2d.hpp : transcription de step6_toolpath2d._offset_walls
// et _generate_infill, offsets par Clipper2 (tolerance declaree, pas le bit).

#include "apex/toolpath2d.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "clipper2/clipper.h"

namespace apex {

namespace {

constexpr int kMaxInfillLines = 5000;           // _MAX_INFILL_LINES
constexpr double kZigzagMaxConnectorSpacings = 2.5;

struct Piece {
    Ring2 outer;               // ouvert (sans point de fermeture)
    std::vector<Ring2> holes;  // ouverts
};

double signed_area(const Ring2& r) {
    double a = 0.0;
    const std::size_t n = r.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Point2& p = r[i];
        const Point2& q = r[(i + 1) % n];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

double piece_area(const Piece& p) {
    double a = std::fabs(signed_area(p.outer));
    for (const auto& h : p.holes) a -= std::fabs(signed_area(h));
    return a;
}

// Retire un point de fermeture eventuel (premier == dernier).
Ring2 open_ring(const Ring2& r) {
    if (r.size() >= 2 && r.front()[0] == r.back()[0] && r.front()[1] == r.back()[1]) {
        return Ring2(r.begin(), r.end() - 1);
    }
    return r;
}

Clipper2Lib::Path64 to_path(const Ring2& r, double scale, bool positive) {
    Clipper2Lib::Path64 path;
    path.reserve(r.size());
    for (const auto& p : r) {
        path.emplace_back(static_cast<int64_t>(std::llround(p[0] * scale)),
                          static_cast<int64_t>(std::llround(p[1] * scale)));
    }
    if ((Clipper2Lib::Area(path) > 0) != positive) std::reverse(path.begin(), path.end());
    return path;
}

Ring2 from_path(const Clipper2Lib::Path64& path, double scale) {
    Ring2 r;
    r.reserve(path.size());
    for (const auto& p : path) r.push_back({static_cast<double>(p.x) / scale, static_cast<double>(p.y) / scale});
    return r;
}

void collect_pieces(const Clipper2Lib::PolyPath64& outer_node, double scale, std::vector<Piece>& out) {
    Piece piece;
    piece.outer = from_path(outer_node.Polygon(), scale);
    for (std::size_t i = 0; i < outer_node.Count(); ++i) {
        const auto* hole = outer_node.Child(i);
        piece.holes.push_back(from_path(hole->Polygon(), scale));
        for (std::size_t j = 0; j < hole->Count(); ++j) collect_pieces(*hole->Child(j), scale, out);
    }
    out.push_back(std::move(piece));
}

// shapely `buffer(delta, join_style="mitre")` d'un ensemble de morceaux.
std::vector<Piece> offset_pieces(const std::vector<Piece>& pieces, double delta, const FillRegionParams& P) {
    Clipper2Lib::Paths64 paths;
    for (const auto& p : pieces) {
        paths.push_back(to_path(p.outer, P.scale, true));
        for (const auto& h : p.holes) paths.push_back(to_path(h, P.scale, false));
    }
    Clipper2Lib::ClipperOffset co(P.miter_limit);
    co.AddPaths(paths, Clipper2Lib::JoinType::Miter, Clipper2Lib::EndType::Polygon);
    Clipper2Lib::PolyTree64 tree;
    co.Execute(delta * P.scale, tree);
    std::vector<Piece> out;
    for (std::size_t i = 0; i < tree.Count(); ++i) collect_pieces(*tree.Child(i), P.scale, out);
    return out;
}

// Anneau ferme a l'orientation de GEOS (exterieur horaire, trou anti-horaire).
Ring2 closed_oriented(Ring2 r, bool exterior) {
    const double a = signed_area(r);
    if ((exterior && a > 0) || (!exterior && a < 0)) std::reverse(r.begin(), r.end());
    if (!r.empty()) r.push_back(r.front());
    return r;
}

struct Bounds {
    double minx = std::numeric_limits<double>::infinity(), miny = minx;
    double maxx = -minx, maxy = -minx;
};

Bounds bounds_of(const std::vector<Piece>& pieces) {
    Bounds b;
    for (const auto& p : pieces) {
        for (const auto& q : p.outer) {
            b.minx = std::min(b.minx, q[0]); b.maxx = std::max(b.maxx, q[0]);
            b.miny = std::min(b.miny, q[1]); b.maxy = std::max(b.maxy, q[1]);
        }
    }
    return b;
}

std::vector<Ring2> all_rings(const std::vector<Piece>& pieces) {
    std::vector<Ring2> rings;
    for (const auto& p : pieces) {
        rings.push_back(p.outer);
        for (const auto& h : p.holes) rings.push_back(h);
    }
    return rings;
}

using Scan = std::vector<Ring2>;   // segments d'une ligne de balayage, dans le sens du balayage

// _clipped_scanlines : lignes paralleles coupees par la zone (pair-impair sur
// toutes les aretes, zone valide), regroupees par ligne, triees le long du
// balayage.
std::vector<Scan> clipped_scanlines(const std::vector<Piece>& zone, double density, double nozzle_width,
                                    double angle, double& spacing) {
    std::vector<Scan> scans;
    spacing = 0.0;
    const Bounds b = bounds_of(zone);
    if (!(b.maxx >= b.minx)) return scans;
    const double diagonal = std::hypot(b.maxx - b.minx, b.maxy - b.miny);
    if (diagonal <= 0) return scans;
    spacing = nozzle_width / std::min(density, 1.0);
    spacing = std::max({spacing, 2.0 * diagonal / kMaxInfillLines, 1e-6});
    const double cx = (b.minx + b.maxx) / 2.0, cy = (b.miny + b.maxy) / 2.0;
    const double half = diagonal;
    const double dx = std::cos(angle), dy = std::sin(angle);
    const double px = -dy, py = dx;
    const auto rings = all_rings(zone);
    const auto count = static_cast<std::int64_t>(std::ceil((2.0 * half) / spacing));
    std::vector<double> ts;
    for (std::int64_t i = 0; i < count; ++i) {
        const double off = -half + static_cast<double>(i) * spacing;
        const double bx = cx + px * off, by = cy + py * off;
        ts.clear();
        for (const auto& r : rings) {
            const std::size_t n = r.size();
            for (std::size_t k = 0; k < n; ++k) {
                const Point2& p = r[k];
                const Point2& q = r[(k + 1) % n];
                const double sp = (p[0] - bx) * px + (p[1] - by) * py;
                const double sq = (q[0] - bx) * px + (q[1] - by) * py;
                if ((sp >= 0) == (sq >= 0)) continue;
                const double tp = (p[0] - bx) * dx + (p[1] - by) * dy;
                const double tq = (q[0] - bx) * dx + (q[1] - by) * dy;
                ts.push_back(tp + (tq - tp) * (sp / (sp - sq)));
            }
        }
        if (ts.size() < 2) continue;
        std::sort(ts.begin(), ts.end());
        Scan scan;
        for (std::size_t k = 0; k + 1 < ts.size(); k += 2) {
            const double t0 = std::max(ts[k], -half), t1 = std::min(ts[k + 1], half);
            if (t1 - t0 <= 1e-12) continue;
            scan.push_back({{bx + dx * t0, by + dy * t0}, {bx + dx * t1, by + dy * t1}});
        }
        if (!scan.empty()) scans.push_back(std::move(scan));
    }
    return scans;
}

double dist(const Point2& a, const Point2& b) { return std::hypot(a[0] - b[0], a[1] - b[1]); }

// Distance d'un point a un anneau ferme et abscisse curviligne du point le
// plus proche (shapely `distance` / `project`).
void project_on_ring(const Ring2& ring, const Point2& p, double& distance, double& s, double& length) {
    distance = std::numeric_limits<double>::infinity();
    s = 0.0;
    double acc = 0.0;
    for (std::size_t k = 0; k + 1 < ring.size(); ++k) {
        const Point2& a = ring[k];
        const Point2& b = ring[k + 1];
        const double ex = b[0] - a[0], ey = b[1] - a[1];
        const double l2 = ex * ex + ey * ey;
        const double seg = std::sqrt(l2);
        double t = l2 > 0 ? ((p[0] - a[0]) * ex + (p[1] - a[1]) * ey) / l2 : 0.0;
        t = std::clamp(t, 0.0, 1.0);
        const double d = std::hypot(a[0] + ex * t - p[0], a[1] + ey * t - p[1]);
        if (d < distance) {
            distance = d;
            s = acc + t * seg;
        }
        acc += seg;
    }
    length = acc;
}

// shapely.ops.substring(ring, s0, s1) avec s0 <= s1 : point a s0, sommets
// strictement entre, point a s1.
Ring2 substring(const Ring2& ring, double s0, double s1) {
    Ring2 out;
    double acc = 0.0;
    auto at = [&](std::size_t k, double local) {
        const Point2& a = ring[k];
        const Point2& b = ring[k + 1];
        const double seg = dist(a, b);
        const double t = seg > 0 ? local / seg : 0.0;
        return Point2{a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t};
    };
    bool started = false;
    for (std::size_t k = 0; k + 1 < ring.size(); ++k) {
        const double seg = dist(ring[k], ring[k + 1]);
        if (!started && s0 <= acc + seg) {
            out.push_back(at(k, s0 - acc));
            started = true;
        }
        if (started) {
            if (s1 <= acc + seg) {
                out.push_back(at(k, s1 - acc));
                return out;
            }
            if (acc + seg > s0) out.push_back(ring[k + 1]);
        }
        acc += seg;
    }
    if (!started && !ring.empty()) out.push_back(ring.back());
    if (out.size() < 2 && !ring.empty()) out.push_back(ring.back());
    return out;
}

// _boundary_connector : points intermediaires du plus court chemin le long
// d'un bord portant a et b. `found` faux = aucun anneau ne les porte, ou
// chemin plus long que max_length.
bool boundary_connector(const Point2& a, const Point2& b, const std::vector<Ring2>& rings, double max_length,
                        Ring2& path_out) {
    const double tolerance = 1e-6 * std::max(1.0, max_length);
    for (const auto& ring : rings) {
        double da, sa, db, sb, length;
        project_on_ring(ring, a, da, sa, length);
        if (da > tolerance) continue;
        project_on_ring(ring, b, db, sb, length);
        if (db > tolerance) continue;
        auto pymod = [](double x, double m) {
            double r = std::fmod(x, m);
            if (r < 0) r += m;
            return r;
        };
        const double forward = pymod(sb - sa, length);
        const double backward = pymod(sa - sb, length);
        if (std::min(forward, backward) > max_length) return false;
        auto along = [&](double s0, double s1) {
            if (s1 >= s0) return substring(ring, s0, s1);
            Ring2 first = substring(ring, s0, length);
            Ring2 second = substring(ring, 0.0, s1);
            first.insert(first.end(), second.begin(), second.end());
            return first;
        };
        Ring2 path;
        if (forward <= backward) {
            path = along(sa, sb);
        } else {
            path = along(sb, sa);
            std::reverse(path.begin(), path.end());
        }
        path_out.clear();
        if (path.size() > 2) path_out.assign(path.begin() + 1, path.end() - 1);
        return true;
    }
    return false;
}

double orient(const Point2& a, const Point2& b, const Point2& c) {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

bool inside_rings(const std::vector<Ring2>& rings, const Point2& p) {
    bool inside = false;
    for (const auto& r : rings) {
        for (std::size_t k = 0; k + 1 < r.size(); ++k) {
            const Point2& a = r[k];
            const Point2& b = r[k + 1];
            if ((a[1] > p[1]) != (b[1] > p[1])) {
                const double x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
                if (p[0] < x) inside = !inside;
            }
        }
    }
    return inside;
}

double point_ring_distance(const std::vector<Ring2>& rings, const Point2& p) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& r : rings) {
        double d, s, l;
        project_on_ring(r, p, d, s, l);
        best = std::min(best, d);
    }
    return best;
}

// `inside.covers(LineString([a, b]))`, inside = zone dilatee de eps :
// aucune traversee franche d'un bord loin des extremites, et milieu dans la
// zone (ou a moins de eps de son bord).
bool segment_covered(const std::vector<Ring2>& rings, const Point2& a, const Point2& b, double eps) {
    for (const auto& r : rings) {
        for (std::size_t k = 0; k + 1 < r.size(); ++k) {
            const Point2& c = r[k];
            const Point2& d = r[k + 1];
            const double o1 = orient(a, b, c), o2 = orient(a, b, d);
            const double o3 = orient(c, d, a), o4 = orient(c, d, b);
            if (((o1 > 0) != (o2 > 0)) && ((o3 > 0) != (o4 > 0)) && o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0) {
                const double t = o3 / (o3 - o4);
                const Point2 x{a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t};
                if (dist(x, a) > eps && dist(x, b) > eps) return false;
            }
        }
    }
    const Point2 mid{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2};
    return inside_rings(rings, mid) || point_ring_distance(rings, mid) <= eps;
}

std::vector<Ring2> connect_zigzag(const std::vector<Scan>& scans, const std::vector<Piece>& zone, double spacing) {
    const double eps = std::max(1e-6, spacing * 1e-4);
    const double max_connector = kZigzagMaxConnectorSpacings * spacing;
    std::vector<Ring2> rings;
    for (const auto& p : zone) {
        rings.push_back(closed_oriented(p.outer, true));
        for (const auto& h : p.holes) rings.push_back(closed_oriented(h, false));
    }
    std::vector<std::vector<char>> used(scans.size());
    for (std::size_t k = 0; k < scans.size(); ++k) used[k].assign(scans[k].size(), 0);
    std::vector<Ring2> chains;
    Ring2 connector;
    for (std::size_t k = 0; k < scans.size(); ++k) {
        for (std::size_t j = 0; j < scans[k].size(); ++j) {
            if (used[k][j]) continue;
            used[k][j] = 1;
            Ring2 chain = scans[k][j];
            Point2 end = chain.back();
            std::size_t current = k;
            while (current + 1 < scans.size()) {
                bool have = false;
                double best_d = 0;
                std::size_t best_j = 0;
                bool best_rev = false;
                const auto& next = scans[current + 1];
                for (std::size_t jj = 0; jj < next.size(); ++jj) {
                    if (used[current + 1][jj]) continue;
                    for (int rev = 0; rev < 2; ++rev) {
                        const Point2& start = rev ? next[jj].back() : next[jj].front();
                        const double d = dist(start, end);
                        if (!have || d < best_d) {
                            have = true; best_d = d; best_j = jj; best_rev = rev != 0;
                        }
                    }
                }
                if (!have || best_d > max_connector) break;
                const Ring2& cand = next[best_j];
                const Point2 start = best_rev ? cand.back() : cand.front();
                if (!boundary_connector(end, start, rings, max_connector, connector)) {
                    if (best_d > 0 && !segment_covered(rings, end, start, eps)) break;
                    connector.clear();
                }
                used[current + 1][best_j] = 1;
                chain.insert(chain.end(), connector.begin(), connector.end());
                if (best_rev) chain.insert(chain.end(), cand.rbegin(), cand.rend());
                else chain.insert(chain.end(), cand.begin(), cand.end());
                end = chain.back();
                current += 1;
            }
            chains.push_back(std::move(chain));
        }
    }
    return chains;
}

void generate_infill(const std::vector<Piece>& boundary, double density, InfillKind pattern, double angle,
                     const FillRegionParams& P, std::vector<Ring2>& out) {
    if (boundary.empty() || density <= 0) return;
    if (pattern == InfillKind::Grid) {
        const double half = std::min(density, 1.0) / 2.0;
        generate_infill(boundary, half, InfillKind::Lines, angle, P, out);
        generate_infill(boundary, half, InfillKind::Lines, angle + 3.14159265358979323846 / 2.0, P, out);
        return;
    }
    double spacing = 0.0;
    if (pattern == InfillKind::CrossedZigzag) {
        const auto zone = offset_pieces(boundary, -P.nozzle_width / 2.0, P);
        if (zone.empty()) return;
        const auto scans = clipped_scanlines(zone, density, P.nozzle_width, angle, spacing);
        auto chains = connect_zigzag(scans, zone, spacing);
        for (auto& c : chains) out.push_back(std::move(c));
        return;
    }
    const auto scans = clipped_scanlines(boundary, density, P.nozzle_width, angle, spacing);
    for (const auto& scan : scans)
        for (const auto& seg : scan) out.push_back(seg);
}

}  // namespace

FillRegionResult fill_region(const Ring2& outer, const std::vector<Ring2>& holes, const FillRegionParams& P) {
    FillRegionResult result;
    Piece base;
    base.outer = open_ring(outer);
    for (const auto& h : holes) base.holes.push_back(open_ring(h));
    if (base.outer.size() < 3) return result;
    const double reference_area = piece_area(base);
    if (!(reference_area > 0)) return result;
    result.base_valid = true;

    std::vector<Piece> innermost;
    bool have_innermost = false;
    for (int n = P.wall_count; n >= 1; --n) {
        auto offset = offset_pieces({base}, -P.nozzle_width * n, P);
        std::vector<Piece> kept;
        for (auto& p : offset) {
            if (piece_area(p) > P.degenerate_area_ratio * reference_area) kept.push_back(std::move(p));
        }
        if (kept.empty()) continue;
        for (const auto& p : kept) {
            result.walls.push_back(closed_oriented(p.outer, true));
            for (const auto& h : p.holes) result.walls.push_back(closed_oriented(h, false));
        }
        if (!have_innermost) {
            innermost = kept;
            have_innermost = true;
        }
    }
    if (have_innermost) generate_infill(innermost, P.infill_density, P.pattern, P.sweep_angle_rad, P, result.infill);
    return result;
}

}  // namespace apex
