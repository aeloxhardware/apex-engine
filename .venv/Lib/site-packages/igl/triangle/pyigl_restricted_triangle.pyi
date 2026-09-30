from typing import Annotated

from numpy.typing import ArrayLike
import scipy.sparse


def cdt(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], E: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], flags: str = '-c') -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Constrained Delaunay triangulation of the convex hull of 2D points and segments.

    Preprocesses input to remove duplicates before calling triangle. Unlike
    triangulate(), this returns an adjusted edge list on the output.

    @param[in] V  #V by 2 list of 2D vertex positions
    @param[in] E  #E by 2 list of constraint edge indices into V
    @param[in] flags  triangle option string; should contain "-c" unless segments
        already enclose all points (default: "-c")
    @return tuple (WV, WF, WE, J):
      - WV  #WV by 2 background mesh vertex positions
      - WF  #WF by 3 background mesh triangle indices into WV
      - WE  #WE by 2 constraint edge indices into WV (may be smaller than E
            if degenerate constraints were removed)
      - J   #V list of indices into WF/WE for each original vertex in V
    """

def refine(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], E: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], flags: str = 'q') -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Refine an existing 2D triangulation.

    @param[in] V  #V by 2 list of 2D vertex positions
    @param[in] E  #E by 2 list of segment indices (constraint edges)
    @param[in] F  #F by 3 list of triangle indices
    @param[in] flags  triangle option string (default: "q" for quality refinement)
    @return tuple (V2, F2):
      - V2  #V2 by 2 refined vertex positions
      - F2  #F2 by 3 refined triangle indices
    """

def remesh_at_points(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], B: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FI: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Given a set of unique points on a mesh specified by barycentric coordinates
    and a triangle index, remesh the mesh to include these points as vertices.
    Barycentric coordinates should be non-negative and sum to 1. Vertex-points are
    not inserted (but tracked in K). Edge-points are inserted and the edge is split
    on all incident faces. Face-points are inserted. The output preserves the
    (non-)manifoldness of the input. The input mesh can be in 3D (or any dimension).

    @param[in] V  #V by dim list of mesh vertex positions
    @param[in] F  #F by 3 list of triangle indices into rows of V
    @param[in] B  #B by 3 list of barycentric coordinates, ith row are coordinates
                  of ith sampled point in face FI(i)
    @param[in] FI  #B list of indices into F
    @param[out] VV  #VV by dim list of mesh vertex positions (top #V rows is always
                    V; bottom rows are new vertices)
    @param[out] FF  #FF by 3 list of triangle indices into rows of VV
    @param[out] J  #FF list of indices into F
    @param[out] K  #B list of indices into VV
    """

class SCAFData:
    def __init__(self) -> None: ...

def scaf_precompute(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')], F: Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')], V_init: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')], scaf_energy: "igl::MappingEnergyType", b: Annotated[ArrayLike, dict(dtype='int32', shape=(None), order='C')], bc: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')], soft_p: float, data: SCAFData) -> None:
    """
    Compute necessary information to start using SCAF

    @param[in] V           #V by 3 list of mesh vertex positions
    @param[in] F           #F by 3/3 list of mesh faces (triangles/tets)
    @param[in] V_init      #V by 3 list of initial mesh vertex positions
    @param[in,out] data  resulting precomputed data
    @param[in] slim_energy Energy type to minimize
    @param[in] b           list of boundary indices into V (soft constraint)
    @param[in] bc          #b by dim list of boundary conditions (soft constraint)
    @param[in] soft_p      Soft penalty factor (can be zero)
    """

def scaf_solve(iter_num: int, data: SCAFData) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')]:
    """
    Run iter_num iterations of SCAF, with precomputed data
    @param[in] data  precomputed data
    @param[in] iter_num  number of iterations to run
    @returns resulting V_o (in SLIMData): #V by dim list of mesh vertex positions
    """

def scaf_system(s: SCAFData) -> tuple[scipy.sparse.csc_matrix[float], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]]:
    """
    Set up the SCAF system L * uv = rhs, without solving it.
    @param[in] s:   igl::SCAFData. Will be modified by energy and Jacobian computation.
    @param[out] L:   m by m matrix
    @param[out] rhs: m by 1 vector
    with m = dim * (#V_mesh + #V_scaf - #V_frame)
    """

def triangulate(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], E: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)] = ..., H: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., VM: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)] = ..., EM: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)] = ..., flags: str = '') -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Triangulate the interior of a polygon using the triangle library.

    @param[in] V #V by 2 list of 2D vertex positions
    @param[in] E #E by 2 list of vertex ids forming unoriented edges of the boundary of the polygon
    @param[in] H #H by 2 coordinates of points contained inside holes of the polygon
    @param[in] VM #V list of markers for input vertices
    @param[in] EM #E list of markers for input edges
    @param[in] flags  string of options pass to triangle (see triangle documentation)
    @param[out] V2  #V2 by 2  coordinates of the vertives of the generated triangulation
    @param[out] F2  #F2 by 3  list of indices forming the faces of the generated triangulation
    @param[out] VM2  #V2 list of markers for output vertices
    @param[out] E2  #E2 by 2 list of output edges
    @param[out] EM2  #E2 list of markers for output edges
    """
