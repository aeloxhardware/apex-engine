from typing import Annotated, overload

from numpy.typing import ArrayLike


def convex_hull(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]:
    """
    Compute the convex hull of a set of points, returning only the triangular faces of the hull.

        @param[in] V  #V by 3 matrix of input points
        @return F: #F by 3 matrix of triangle indices into V
    """

def coplanar(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> bool:
    """
    Test whether all points lie on the same plane.

    @param[in] V  #V by 3 list of 3D vertex positions
    @return true if all points are coplanar
    """

def delaunay_triangulation(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]:
    """
    Compute the Delaunay triangulation of a 2D point set using CGAL.

    @param[in] V  #V by 2 list of 2D vertex positions
    @return F  #F by 3 triangle indices into V
    """

def extract_cells(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[int, Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Extract connected 3D cells partitioned by a triangle mesh.

    @param[in] V  #V by 3 array of vertices
    @param[in] F  #F by 3 array of triangle indices
    @return tuple (num_cells, cells):
      - num_cells  number of cells (cell 0 is the infinite outer cell)
      - cells  #F by 2 array of cell indices; cells[i,0] is the cell on the
        positive side of face i, cells[i,1] is on the negative side
    """

def fast_winding_number(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], N: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], expansion_order: int = 2, beta: float = 2.0) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]:
    """
    Evaluate the fast winding number for point data with adjustable accuracy.

        @param[in] P  #P by 3 list of point locations
        @param[in] N  #P by 3 list of point normals
        @param[in] Q  #Q by 3 list of query points for the winding number
        @param[in] expansion_order  Order of the Taylor expansion (0, 1, or 2)
        @param[in] beta  Barnes-Hut style accuracy parameter (recommended: 2)
        @return Vector of winding number values for each query point
    """

def intersect_other(VA: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FA: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], VB: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FB: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], detect_only: bool = False, first_only: bool = False, stitch_all: bool = False, slow_and_more_precise_rounding: bool = False, cutoff: int = 1000) -> tuple[Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Detect intersecting faces between two triangle meshes, providing detailed output.

        @param[in] VA  #V by 3 list of vertices for first mesh
        @param[in] FA  #F by 3 list of faces for first mesh
        @param[in] VB  #V by 3 list of vertices for second mesh
        @param[in] FB  #F by 3 list of faces for second mesh
        @param[in] detect_only  only detect intersections, do not resolve
        @param[in] first_only  only return first intersection
        @param[in] stitch_all  stitch all intersections
        @param[in] slow_and_more_precise_rounding  use slow and more precise rounding
        @param[in] cutoff  maximum number of intersections to resolve
        @return Tuple containing:
          - success: bool indicating if the operation succeeded
          - IF: # intersecting face pairs
          - VVAB: list of intersection vertex positions
          - FFAB: list of triangle indices into VVAB
          - JAB: list of indices into [FA;FB] denoting the birth triangle
          - IMAB: indices stitching duplicates from intersections
    """

@overload
def intersect_with_half_space(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], p: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)], n: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Intersect a PWN mesh with a half-space using a point and normal.

        @param[in] V  #V by 3 list of mesh vertex positions
        @param[in] F  #F by 3 list of triangle indices
        @param[in] p  3D point on plane
        @param[in] n  3D normal vector
        @return Tuple containing:
          - success: bool, true if successful
          - VC: vertices of resulting mesh
          - FC: face indices of resulting mesh
          - J: birth facet indices
    """

@overload
def intersect_with_half_space(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], equ: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Intersect a PWN mesh with a half-space using the plane equation.

        @param[in] V  #V by 3 list of mesh vertex positions
        @param[in] F  #F by 3 list of triangle indices
        @param[in] equ  Plane equation coefficients (a, b, c, d)
        @return Tuple containing:
          - success: bool, true if successful
          - VC: vertices of resulting mesh
          - FC: face indices of resulting mesh
          - J: birth facet indices
    """

def mesh_boolean(VA: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FA: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], VB: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., FB: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)] = ..., type_str: str) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Compute the boolean operation (union, intersection, difference, etc.) between two meshes.

    @param[in] VA #VA by dim matrix of mesh A vertices
    @param[in] FA #FA by simplex_size matrix of mesh A faces
    @param[in] VB #VB by dim matrix of mesh B vertices
    @param[in] FB #FB by simplex_size matrix of mesh B faces
    @param[in] type_str Type of boolean operation: "union", "intersection", "difference", etc.
    @param[out] VC #VC by dim matrix of result vertices
    @param[out] FC #FC by simplex_size matrix of result faces
    @param[out] J #FC list of indices indicating which input face contributed to each result face
    @return Tuple containing:
      - VC: Result vertices
      - FC: Result faces
      - J: Face origin indices
    """

def oriented_bounding_box(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]:
    """
    Given a set of points compute the rotation transformation of them such
    that their axis-aligned bounding box is as small as possible.

    igl::oriented_bounding_box is often faster and better

      @param[in] P  #P by 3 list of point locations
      @param[out] R  rotation matrix
    """

def outer_hull(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='bool', shape=(None), order='C')]]:
    """
    Compute the outer hull of a piecewise-constant winding number mesh.

    @param[in] V   #V by 3 list of vertex positions
    @param[in] F   #F by 3 list of triangle indices into V
    @return tuple (HV, HF, J, flip):
      - HV    #HV by 3 output vertex positions
      - HF    #HF by 3 output triangle indices into HV
      - J     #HF list of indices into F (birth faces)
      - flip  #HF list of bools, true if facet was flipped
    """

def peel_outer_hull_layers(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='bool', shape=(None), order='C')], int]:
    """
    Successively peel outer hull layers from a mesh, useful for boolean operations.

    @param[in] V  #V by 3 list of vertex positions
    @param[in] F  #F by 3 list of triangle indices into V
    @return tuple (I, flip, num_peels):
      - I         #F list of which peel layer each facet belongs to
      - flip      #F list of bools, true if facet orientation was flipped during peeling
      - num_peels number of peel layers
    """

def peel_winding_number_layers(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], int]:
    """
    Peel winding number layers from a mesh.

    @param[in] V  #V by 3 list of vertex positions
    @param[in] F  #F by 3 list of triangle indices into V
    @return tuple (W, num_layers):
      - W          #V list of integer winding numbers per vertex
      - num_layers number of winding number layers
    """

def point_areas(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], I: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], N: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    r"""
    Given a 3D set of points P, each with a list of k-nearest-neighbours,
    estimate the geodesic voronoi area associated with each point.

    The k nearest neighbours may be known from running igl::knn_octree on
    the output data from igl::octree. We reccomend using a k value
    between 15 and 20 inclusive for accurate area estimation.

    N is used filter the neighbours, to ensure area estimation only occurs
    using neighbors that are on the same side of the surface (ie for thin
    sheets), as well as to solve the orientation ambiguity of the tangent
    plane normal.

    \note This function *should* be implemented by pre-filtering I, rather
    than filtering in this function using N. In this case, the function
    would only take P and I as input.

    @param[in] P  #P by 3 list of point locations
    @param[in] I  #P by k list of k-nearest-neighbor indices into P
    @param[in] N  #P by 3 list of point normals
    @param[out] A  #P list of estimated areas
    @param[out] T  #P by 3 list of tangent plane normals for each point

    \see igl::knn
    """

def remesh_self_intersections(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], detect_only: bool = False, first_only: bool = False, stitch_all: bool = False, slow_and_more_precise_rounding: bool = False, cutoff: int = 1000) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Resolve self-intersections in a mesh, without returning unique vertex indices (IM).

        @param[in] V  #V by 3 list of vertex positions
        @param[in] F  #F by 3 list of face indices
        @param[in] detect_only  only detect intersections, do not resolve
        @param[in] first_only  only return first intersection
        @param[in] stitch_all  stitch all intersections
        @param[in] slow_and_more_precise_rounding  use slow and more precise rounding
        @param[in] cutoff  maximum number of intersections to resolve
        @return Tuple containing:
          - VV: remeshed vertex positions
          - FF: remeshed face indices
          - IF: intersecting face pairs
          - J: birth triangle indices
          - IM if stitch_all = true   #VV list from 0 to #VV-1
               elseif stitch_all = false  #VV list of indices into VV of unique vertices.
    """

def trim_with_solid(VA: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FA: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], VB: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], FB: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='bool', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Trim a mesh with another solid mesh, determining which faces lie inside or outside.

        @param[in] VA Vertex positions of mesh A
        @param[in] FA Triangle indices of mesh A
        @param[in] VB Vertex positions of mesh B (solid)
        @param[in] FB Triangle indices of mesh B
        @param[out] V Output vertex positions
        @param[out] F Output triangle indices
        @param[out] D Boolean vector indicating if each face is inside B
        @param[out] J Indices into FA showing parent triangle
    """
