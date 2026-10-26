#include "tetranet/tetranet.hpp"
#include <cstdio>
#include <cstring>

int main() {
    if (std::strcmp(tetranet::version(), "0.0.0") != 0) return 1;
    tetranet::SurfaceMesh s;
    tetranet::TetMesh t;
    auto r = tetranet::mesh(s, tetranet::Params{}, t);
    std::puts(r.reason.c_str());
    return r.status == tetranet::Status::InvalidInput ? 0 : 1;
}
