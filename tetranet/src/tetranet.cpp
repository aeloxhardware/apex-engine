#include "tetranet/tetranet.hpp"

namespace tetranet {

const char* version() { return "0.0.0"; }

Report mesh(const SurfaceMesh&, const Params&, TetMesh&) {
    Report r;
    r.status = Status::InvalidInput;
    r.reason = "TetraNet n'est pas encore implémenté.";
    return r;
}

}  // namespace tetranet
