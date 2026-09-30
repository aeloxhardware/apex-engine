from typing import Annotated, overload

from numpy.typing import ArrayLike


@overload
def box_cubic(C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]]:
    """
    Compute the min/max box corners tightly containing a cubic Bézier curve.

    @param[in] C  4 by dim matrix of control points defining the cubic Bézier curve
    @return Tuple (B1, B2) where
      B1  1 by dim min corner of the bounding box
      B2  1 by dim max corner of the bounding box
    """

@overload
def box_cubic(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    """
    Compute bounding boxes for a collection of indexed cubic Bézier curves.

    @param[in] P  #P by dim matrix of control point locations
    @param[in] C  #C by 4 matrix of indices into P defining the cubics
    @return Tuple (B1, B2) where
      B1  #C by dim matrix of min corners of the bounding boxes
      B2  #C by dim matrix of max corners of the bounding boxes
    """

@overload
def point_cubic_squared_distance(Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    """
    Squared distance from each query point to a cubic Bézier curve.

    @param[in] Q  #Q by dim matrix of query points
    @param[in] C  4 by dim matrix of control points for the cubic Bézier curve
    @return Tuple (sqrD, S, K) where
      sqrD  #Q vector of smallest squared distances
      S     #Q vector of parameters of the closest points on the curve
      K     #Q by dim matrix of closest points on the curve
    """

@overload
def point_cubic_squared_distance(Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], D: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], B: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    """
    Squared distance to a cubic Bézier curve, reusing precomputed monomial bases.

    Compute the bases once with igl.cubic_monomial_bases(C) -> (M, D, B) and pass
    (D, B) here so repeated queries against the same curve skip recomputing them.
    The result matches the two-argument overload.

    @param[in] Q  #Q by dim matrix of query points
    @param[in] C  4 by dim matrix of control points for the cubic Bézier curve
    @param[in] D  3 by dim matrix of monomial coefficients for dC/dt (from cubic_monomial_bases)
    @param[in] B  6-vector of inner products of the monomial bases (from cubic_monomial_bases)
    @return Tuple (sqrD, S, K) where
      sqrD  #Q vector of smallest squared distances
      S     #Q vector of parameters of the closest points on the curve
      K     #Q by dim matrix of closest points on the curve
    """

@overload
def point_spline_squared_distance(Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    """
    Squared distance from each query point to a spline of cubic Bézier curves.

    @param[in] Q  #Q by dim matrix of query points
    @param[in] P  #P by dim matrix of spline control points
    @param[in] C  #C by 4 matrix of indices into P defining the cubic Bézier curves
    @return Tuple (sqrD, I, S, K) where
      sqrD  #Q vector of smallest squared distances
      I     #Q vector of indices of the closest cubic (row of C)
      S     #Q vector of parameters of the closest points on that cubic
      K     #Q by dim matrix of closest points on the spline
    """

@overload
def point_spline_squared_distance(Q: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], B1: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], B2: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], leaf: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]]:
    """
    Squared distance to a spline, reusing a prebuilt Eytzinger AABB.

    Build the AABB once with igl.cycodebase.spline_eytzinger_aabb(P, C) and pass its
    (B1, B2, leaf) here so repeated queries against the same spline skip rebuilding
    the acceleration structure. The result matches the unaccelerated overload.

    @param[in] Q     #Q by dim matrix of query points
    @param[in] P     #P by dim matrix of spline control points
    @param[in] C     #C by 4 matrix of indices into P defining the cubic Bézier curves
    @param[in] B1    #B by dim matrix of AABB min box corners (from spline_eytzinger_aabb)
    @param[in] B2    #B by dim matrix of AABB max box corners (from spline_eytzinger_aabb)
    @param[in] leaf  #B vector of AABB leaf node indices/flags (from spline_eytzinger_aabb)
    @return Tuple (sqrD, I, S, K) where
      sqrD  #Q vector of smallest squared distances
      I     #Q vector of indices of the closest cubic (row of C)
      S     #Q vector of parameters of the closest points on that cubic
      K     #Q by dim matrix of closest points on the spline
    """

def roots(coef: Annotated[ArrayLike, dict(dtype='float64', shape=(None), writable=False)], xmin: float, xmax: float) -> tuple[int, Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]]:
    """
    Find the real roots of a polynomial within an interval [xmin, xmax].

    @param[in] coef  #coef list of monomial coefficients (low to high order); the
      polynomial degree is len(coef)-1
    @param[in] xmin  lower bound of the search interval
    @param[in] xmax  upper bound of the search interval
    @return Tuple (n, R) where
      n  number of roots found in [xmin, xmax]
      R  degree-length vector whose first n entries are the roots (ascending); the
         remaining entries are NaN
    """

def spline_eytzinger_aabb(P: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')]]:
    """
    Compute an Eytzinger-layout AABB tree for a spline of cubic Bézier curves.

    @param[in] P  #P by dim matrix of spline control points
    @param[in] C  #C by 4 matrix of indices into P defining the cubic Bézier curves
    @return Tuple (B1, B2, leaf) where
      B1    #B by dim matrix of AABB min box corners
      B2    #B by dim matrix of AABB max box corners
      leaf  #B vector of AABB leaf node indices/flags
    """
