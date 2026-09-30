import enum
from typing import Annotated, overload

from numpy.typing import ArrayLike


class Orientation(enum.Enum):
    POSITIVE = 1

    INSIDE = 1

    NEGATIVE = -1

    OUTSIDE = -1

    COLLINEAR = 0

    COPLANAR = 0

    COCIRCULAR = 0

    COSPHERICAL = 0

    DEGENERATE = 0

POSITIVE: Orientation = Orientation.POSITIVE

NEGATIVE: Orientation = Orientation.NEGATIVE

COLLINEAR: Orientation = Orientation.COLLINEAR

def cubic_winding_number(C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], q: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> float:
    """
    Exact winding number of a 2D cubic Bézier curve about a query point.

    @param[in] C  4 by 2 matrix of control points for the cubic Bézier curve
    @param[in] q  2D query point
    @return the (fractional) winding number of the curve about q
    """

def delaunay_triangulation(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]) -> Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')]:
    """
    Compute a Delaunay triangulation of a 2D point set using exact predicates.

    @param[in] V  #V by 2 list of 2D vertex positions
    @return F  #F by 3 list of triangle indices into V.
    """

@overload
def ear_clipping(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]) -> tuple[bool, Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')]]:
    """
    Ear-clipping triangulation of a 2D polygon.

    Reverses P if necessary so output orientation matches input.

    @param[in] P  #P by 2 list of 2D polygon vertices (CCW for proper mesh)
    @return (ok, eF) where ok is true if mesh is proper (input is a simple polygon)
      and eF is #eF by 3 list of triangle indices into P.
    """

@overload
def ear_clipping(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], RT: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]) -> tuple[Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Ear-clipping triangulation with reserved (non-clippable) vertices.

    @param[in] P   #P by 2 list of 2D polygon vertices
    @param[in] RT  #P list, vertices marked 1 are preserved (not clipped)
    @return (eF, I) where eF is the clipped ears (original indices into P) and
      I maps indices of the remaining polygon to original P indices.
    """

def find_intersections(V1: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], F1: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], V2: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], F2: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], first_only: bool = False) -> tuple[bool, Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Find intersecting triangle pairs between two meshes.

    Uses AABB tree and exact predicates.

    @param[in] V1        #V1 by 3 vertex positions of first mesh
    @param[in] F1        #F1 by 3 triangle indices of first mesh
    @param[in] V2        #V2 by 3 vertex positions of second mesh
    @param[in] F2        #F2 by 3 triangle indices of second mesh
    @param[in] first_only  stop after finding the first intersection
    @return (found, IF, CP) where found is true if any intersections exist,
      IF is #IF by 2 list of intersecting face index pairs (into F1, F2),
      CP is #IF list of whether each intersection is coplanar.
    """

def find_self_intersections(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], first_only: bool = False) -> tuple[bool, Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Find self-intersecting triangle pairs within a mesh.

    Uses AABB tree and exact predicates.

    @param[in] V         #V by 3 vertex positions
    @param[in] F         #F by 3 triangle indices
    @param[in] first_only  stop after finding the first intersection
    @return (found, IF, CP) where found is true if any self-intersections exist,
      IF is #IF by 2 list of intersecting face index pairs,
      CP is #IF list of whether each intersection is coplanar.
    """

def incircle(pa: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], pb: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], pc: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], pd: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> Orientation:
    """
    Decide whether a 2D point is inside/outside/on a circle.

    Uses exact arithmetic (Shewchuk predicates).

    @param[in] pa  2D point on circle
    @param[in] pb  2D point on circle
    @param[in] pc  2D point on circle
    @param[in] pd  2D query point
    @return INSIDE if pd is inside the circle defined by pa,pb,pc,
      OUTSIDE if outside, COCIRCULAR if exactly on the circle.
    """

def insphere(pa: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pb: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pc: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pd: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pe: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')]) -> Orientation:
    """
    Decide whether a 3D point is inside/outside/on a sphere.

    Uses exact arithmetic (Shewchuk predicates).

    @param[in] pa  3D point on sphere
    @param[in] pb  3D point on sphere
    @param[in] pc  3D point on sphere
    @param[in] pd  3D point on sphere
    @param[in] pe  3D query point
    @return INSIDE if pe is inside the sphere defined by pa,pb,pc,pd,
      OUTSIDE if outside, COSPHERICAL if exactly on the sphere.
    """

def lexicographic_triangulation(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]) -> Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')]:
    """
    Compute a lexicographic triangulation of a 2D point set using exact predicates.

    @param[in] V  #V by 2 list of 2D vertex positions
    @return F  #F by 3 list of triangle indices into V.
    """

def orient2d(pa: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], pb: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], pc: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> Orientation:
    """
    Compute the orientation of the triangle formed by pa, pb, pc.

    Uses exact arithmetic (Shewchuk predicates).

    @param[in] pa  2D point
    @param[in] pb  2D point
    @param[in] pc  2D point
    @return POSITIVE if pa,pb,pc are counterclockwise, NEGATIVE if clockwise,
      COLLINEAR if collinear.
    """

@overload
def orient3d(pa: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pb: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pc: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], pd: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')]) -> Orientation:
    """
    Compute the orientation of the tetrahedron formed by pa, pb, pc, pd.

    Uses exact arithmetic (Shewchuk predicates).

    @param[in] pa  3D point
    @param[in] pb  3D point
    @param[in] pc  3D point
    @param[in] pd  3D query point
    @return POSITIVE if pd is below the plane oriented by pa,pb,pc,
      NEGATIVE if above, COPLANAR if on the plane.
    """

@overload
def orient3d(A: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], B: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], D: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]) -> Annotated[ArrayLike, dict(dtype='int32', shape=(None), order='C')]:
    """
    Vectorized orient3d: compute orientation for each row-tuple (A[i], B[i], C[i], D[i]).

    @param[in] A  #P by 3 matrix of 3D points
    @param[in] B  #P by 3 matrix of 3D points
    @param[in] C  #P by 3 matrix of 3D points
    @param[in] D  #P by 3 matrix of 3D points
    @return #P vector of orientation values (+1, -1, or 0)
    """

def point_in_convex_hull(q: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], a: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], b: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], c: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], d: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> Orientation:
    """
    Test whether a 2D point lies in the convex hull of four points using exact predicates.

    @param[in] q  2D query point
    @param[in] a  2D point
    @param[in] b  2D point
    @param[in] c  2D point
    @param[in] d  2D point
    @return POSITIVE if q is strictly inside the convex hull of {a,b,c,d},
      NEGATIVE if strictly outside, COLLINEAR if on the boundary.
    """

def point_inside_convex_polygon(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, 2), order='F')], q: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> bool:
    """
    Check whether a 2D point lies inside a 2D convex polygon.

    @param[in] P  #P by 2 list of polygon vertices (n >= 3)
    @param[in] q  2D query point
    @return true if q is inside P.
    """

def polygons_to_triangles(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], I: Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]) -> tuple[Annotated[ArrayLike, dict(dtype='int32', shape=(None, None), order='F')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Triangulate each polygon of a polygon mesh with a fan.

    @param[in] V  #V by dim vertex positions (used for orientation in 3D)
    @param[in] I  #I list of polygon corner indices into rows of V
    @param[in] C  #polygons+1 cumulative polygon sizes (C[i+1]-C[i] = size of polygon i)
    @return (F, J) where F is #F by 3 list of triangle indices and J is #F list of
      source polygon indices.
    """

def segment_segment_intersect(A: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], B: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], C: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')], D: Annotated[ArrayLike, dict(dtype='float64', shape=(2), order='C')]) -> bool:
    """
    Test whether two 2D segments intersect using exact orient2d predicates.

    @param[in] A  1st endpoint of segment 1
    @param[in] B  2nd endpoint of segment 1
    @param[in] C  1st endpoint of segment 2
    @param[in] D  2nd endpoint of segment 2
    @return true if the segments intersect.
    """

def spline_winding_number(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], B1: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], B2: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], leaf: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)], Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]:
    """
    Winding number of a closed spline of cubic Bézier curves about query points.

    Uses the Eytzinger-layout AABB tree produced by
    igl.cycodebase.spline_eytzinger_aabb for acceleration.

    @param[in] P     #P by 2 matrix of spline control points
    @param[in] C     #C by 4 matrix of indices into P defining the cubic Bézier curves
    @param[in] B1    #B by 2 matrix of AABB min box corners
    @param[in] B2    #B by 2 matrix of AABB max box corners
    @param[in] leaf  #B vector of AABB leaf node indices/flags
    @param[in] Q     #Q by 2 matrix of query points
    @return W  #Q vector of winding numbers about each query point
    """

def triangle_triangle_intersect(a1: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], a2: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], a3: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], b1: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], b2: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')], b3: Annotated[ArrayLike, dict(dtype='float64', shape=(3), order='C')]) -> tuple[bool, bool]:
    """
    Test whether two 3D triangles intersect using exact predicates.

    @param[in] a1  1st vertex of triangle A
    @param[in] a2  2nd vertex of triangle A
    @param[in] a3  3rd vertex of triangle A
    @param[in] b1  1st vertex of triangle B
    @param[in] b2  2nd vertex of triangle B
    @param[in] b3  3rd vertex of triangle B
    @return (intersects, coplanar) where intersects is true if the triangles
      intersect, and coplanar is true if they are coplanar.
    """
