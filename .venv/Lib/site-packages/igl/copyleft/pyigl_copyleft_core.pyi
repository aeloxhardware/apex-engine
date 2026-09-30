from typing import Annotated, overload

from numpy.typing import ArrayLike


@overload
def marching_cubes(values: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)], points: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], x_res: int, y_res: int, z_res: int, isovalue: float = 0.0) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Perform marching cubes to extract an isosurface from a scalar field on a regular grid.

    @param[in] values  #points vector of scalar values (<0 inside, >0 outside)
    @param[in] points  #points by 3 grid point positions (x-fastest order:
        index = x + y*x_res + z*x_res*y_res)
    @param[in] x_res  grid resolution in x
    @param[in] y_res  grid resolution in y
    @param[in] z_res  grid resolution in z
    @param[in] isovalue  isovalue to reconstruct (default: 0.0)
    @return tuple (vertices, faces):
      - vertices  #V by 3 mesh vertex positions
      - faces  #F by 3 mesh triangle indices
    """

@overload
def marching_cubes(values: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)], points: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], indices: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], isovalue: float = 0.0) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Perform marching cubes on a sparse grid defined by explicit cube corner indices.

    @param[in] values  #points vector of scalar values
    @param[in] points  #points by 3 grid point positions
    @param[in] indices  #cubes by 8 array of corner indices into points/values
    @param[in] isovalue  isovalue to reconstruct (default: 0.0)
    @return tuple (vertices, faces):
      - vertices  #V by 3 mesh vertex positions
      - faces  #F by 3 mesh triangle indices
    """

def progressive_hulls(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')], F: Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')], max_m: int = 0) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='F')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='F')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Performs progressive hull simplification on a mesh, collapsing edges until a target number of faces is reached.

        @param[in] V  #V by dim list of vertex positions
        @param[in] F  #F by 3 list of face indices into V
        @param[in] max_m Target number of output faces
        @param[out] U Output vertex positions
        @param[out] G Output face indices into U
        @param[out] J Indices into F indicating the birth face for each face in G
    """

def quadprog(G: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], g0: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)], CE: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., ce0: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)] = ..., CI: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., ci0: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)] = ...) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]:
    """
    Solve a dense convex quadratic program.

    Minimizes:  0.5 x' G x + g0' x
    Subject to: CE' x + ce0  = 0
                CI' x + ci0 >= 0

    @param[in] G    #x by #x positive semi-definite quadratic cost matrix
    @param[in] g0   #x linear cost vector
    @param[in] CE   #x by #ce equality constraint matrix (optional)
    @param[in] ce0  #ce equality constraint rhs (optional)
    @param[in] CI   #x by #ci inequality constraint matrix (optional)
    @param[in] ci0  #ci inequality constraint rhs (optional)
    @return x  #x solution vector
    @raises RuntimeError if the solver fails
    """
