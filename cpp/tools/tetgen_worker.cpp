// apex_tetgen_worker : executable SEPARE qui fait tourner UNE tetraedrisation
// TetGen (AGPL-3.0, cpp/third_party/tetgen). Jamais lie dans apex_core :
// isolation des crashs natifs (CLAUDE.md limites #35/#56) et frontiere de
// programme claire. Interface : apex_tetgen_worker <entree.bin> <sortie.bin>.
//
// Entree (little-endian) : int32 nv, nf, nh, nr, slen ; slen octets de
// switches ; nv*3 double ; nf*3 int32 ; nh*3 double ; nr*5 double
// (x, y, z, id, max_vol) - memes donnees que le binding Python `tetgen`.
// Sortie : int32 statut (0 ok, 1 erreur TetGen) ; ok : int32 nn, nt, nn*3
// double, nt*4 int32 ; erreur : int32 n, n octets de message.
// Un crash natif se manifeste par un code de sortie non nul, sans sortie.
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "tetgen.h"

namespace {

template <class T>
bool read_vec(FILE* f, std::vector<T>& v, size_t n) {
  v.resize(n);
  return n == 0 || std::fread(v.data(), sizeof(T), n, f) == n;
}

template <class T>
void write_raw(FILE* f, const T* p, size_t n) {
  if (n) std::fwrite(p, sizeof(T), n, f);
}

int fail(const char* path, const std::string& msg) {
  FILE* f = std::fopen(path, "wb");
  if (!f) return 3;
  int32_t st = 1, n = static_cast<int32_t>(msg.size());
  write_raw(f, &st, 1);
  write_raw(f, &n, 1);
  write_raw(f, msg.data(), msg.size());
  std::fclose(f);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) return 2;
  FILE* in = std::fopen(argv[1], "rb");
  if (!in) return 2;
  int32_t hdr[5];
  if (std::fread(hdr, sizeof(int32_t), 5, in) != 5) return 2;
  const int nv = hdr[0], nf = hdr[1], nh = hdr[2], nr = hdr[3], slen = hdr[4];
  std::vector<char> sw;
  std::vector<double> verts, holes, regions;
  std::vector<int32_t> faces;
  if (!read_vec(in, sw, slen) || !read_vec(in, verts, size_t(nv) * 3) ||
      !read_vec(in, faces, size_t(nf) * 3) || !read_vec(in, holes, size_t(nh) * 3) ||
      !read_vec(in, regions, size_t(nr) * 5))
    return 2;
  std::fclose(in);
  std::string switches(sw.begin(), sw.end());

  tetgenio io, io_out, io_bg, io_unused;
  io.initialize();
  io.numberofpoints = nv;
  io.pointlist = new double[size_t(nv) * 3];
  std::memcpy(io.pointlist, verts.data(), sizeof(double) * size_t(nv) * 3);
  io.numberoffacets = nf;
  io.facetlist = new tetgenio::facet[nf];
  for (int i = 0; i < nf; i++) {
    tetgenio::facet* f = &io.facetlist[i];
    io.init(f);
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[1];
    tetgenio::polygon* p = &f->polygonlist[0];
    io.init(p);
    p->numberofvertices = 3;
    p->vertexlist = new int[3];
    for (int j = 0; j < 3; j++) p->vertexlist[j] = faces[size_t(i) * 3 + j];
  }
  // Meme ordre que le binding : regions puis trous alloues avant l'appel.
  io.numberofregions = nr;
  io.regionlist = new double[size_t(nr) * 5];
  std::memcpy(io.regionlist, regions.data(), sizeof(double) * size_t(nr) * 5);
  io.numberofholes = nh;
  io.holelist = new double[size_t(nh) * 3];
  std::memcpy(io.holelist, holes.data(), sizeof(double) * size_t(nh) * 3);

  tetgenbehavior behavior;
  char buf[1024];
  std::strncpy(buf, switches.c_str(), sizeof(buf) - 1);
  buf[sizeof(buf) - 1] = '\0';
  if (!behavior.parse_commandline(buf)) return fail(argv[2], "Unable to parse switches");

  const auto t0 = std::chrono::steady_clock::now();
  try {
    tetrahedralize(&behavior, &io, &io_out, &io_unused, &io_bg);
  } catch (int code) {
    return fail(argv[2], "TetGen a termine avec le code " + std::to_string(code));
  } catch (const std::exception& e) {
    // Cette variante de TetGen leve des std::runtime_error portant le meme
    // texte que celui relaye par le binding Python (la cascade de step1 y
    // cherche "self-intersection", "recoversubfaces", "split_segment").
    return fail(argv[2], e.what());
  } catch (...) {
    return fail(argv[2], "TetGen a leve une exception inconnue");
  }
  const double secs =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::fprintf(stderr, "apex_tetgen_worker: tetrahedralize %.3f s\n", secs);

  if (io_out.numberofpoints <= 0)
    return fail(argv[2],
                "Failed to tetrahedralize. You may need to repair surface by making it manifold");
  const int nn = io_out.numberofpoints, nt = io_out.numberoftetrahedra;
  std::vector<int32_t> elem(size_t(nt) * 4);
  for (size_t i = 0; i < elem.size(); i++) {
    // corners = 4 ; indices ramenes a base 0 comme le binding
    elem[i] = io_out.tetrahedronlist[i] - (io_out.firstnumber == 1 ? 1 : 0);
  }
  FILE* out = std::fopen(argv[2], "wb");
  if (!out) return 3;
  int32_t st = 0;
  write_raw(out, &st, 1);
  int32_t dims[2] = {nn, nt};
  write_raw(out, dims, 2);
  write_raw(out, io_out.pointlist, size_t(nn) * 3);
  write_raw(out, elem.data(), elem.size());
  std::fclose(out);
  return 0;
}
