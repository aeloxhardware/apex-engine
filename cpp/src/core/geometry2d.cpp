// Voir apex/geometry2d.hpp.

#include "apex/geometry2d.hpp"

namespace apex {

bool point_in_polygon_2d(double x, double y, const double* polygon, std::int64_t n) {
    bool inside = false;
    for (std::int64_t i = 0; i < n; ++i) {
        const double x1 = polygon[i * 2], y1 = polygon[i * 2 + 1];
        const std::int64_t j = (i + 1) % n;
        const double x2 = polygon[j * 2], y2 = polygon[j * 2 + 1];
        if ((y1 > y) != (y2 > y)) {
            const double x_intersect = (x2 - x1) * (y - y1) / (y2 - y1) + x1;
            if (x < x_intersect) {
                inside = !inside;
            }
        }
    }
    return inside;
}

}  // namespace apex
