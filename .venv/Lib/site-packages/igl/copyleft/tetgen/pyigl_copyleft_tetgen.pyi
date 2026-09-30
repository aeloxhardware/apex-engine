from typing import Annotated

from numpy.typing import ArrayLike


class CDTParam:
    def __init__(self) -> None: ...

    @property
    def use_bounding_box(self) -> bool: ...

    @use_bounding_box.setter
    def use_bounding_box(self, arg: bool, /) -> None: ...

    @property
    def bounding_box_scale(self) -> float: ...

    @bounding_box_scale.setter
    def bounding_box_scale(self, arg: float, /) -> None: ...

    @property
    def flags(self) -> str: ...

    @flags.setter
    def flags(self, arg: str, /) -> None: ...

def cdt(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], param: CDTParam = ...) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Create a constrained Delaunay tetrahedralization of the convex hull of a mesh.

    The input mesh must be non-self-intersecting.

    @param[in] V  #V by 3 list of input mesh vertices
    @param[in] F  #F by 3 list of input mesh facets
    @param[in] param  CDTParam struct with fields:
        use_bounding_box    force meshing of the bounding box (default: false)
        bounding_box_scale  scale factor for bounding box (default: 1.01)
        flags               tetgen flags, do not include "c" (default: "Y")
    @return tuple (TV, TT, TF):
      - TV  #TV by 3 output vertex positions (input vertices come first)
      - TT  #TT by 4 tetrahedra indices into TV
      - TF  #TF by 3 surface triangle indices (facets from F, possibly subdivided)
    @raises RuntimeError on failure
    """

def mesh_with_skeleton(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)], C: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], P: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)] = ..., BE: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)] = ..., CE: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)] = ..., samples_per_bone: int = 10, tetgen_flags: str = '') -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')]]:
    """
    Mesh the interior of a surface with a graded tetrahedral mesh conforming to a skeleton.

    @param[in] V  #V by 3 list of surface mesh vertex positions
    @param[in] F  #F by 3 list of triangle indices
    @param[in] C  #C by 3 list of skeleton joint positions
    @param[in] P  #P list of point handle indices into C (optional)
    @param[in] BE #BE by 2 list of bone edge indices into C (optional)
    @param[in] CE #CE by 2 list of cage edge indices into P (optional)
    @param[in] samples_per_bone  number of samples added per bone (default: 10)
    @param[in] tetgen_flags  tetgen option string (default: "" → "pq2Y")
    @return tuple (VV, TT, FF):
      - VV  #VV by 3 tet-mesh vertex positions
      - TT  #TT by 4 tetrahedra indices
      - FF  #FF by 3 surface triangle indices
    @raises RuntimeError on failure
    """

def tetrahedralize(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)] = ..., H: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., VM: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)] = ..., FM: Annotated[ArrayLike, dict(dtype='int64', shape=(None), writable=False)] = ..., R: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)] = ..., flags: str = '') -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None), order='C')], Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), order='C')], int]:
    """
    Mesh the interior of a surface mesh (V,F) using tetgen

    @param[in] V  #V by 3 vertex position list
    @param[in] F  #F list of polygon face indices into V (0-indexed)
    @param[in] H  #H by 3 list of seed points inside holes
    @param[in] VM  #VM list of vertex markers
    @param[in] FM  #FM list of face markers
    @param[in] R  #R by 5 list of region attributes            
    @param[in] flags string of tetgen options (See tetgen documentation) e.g.
        "pq1.414a0.01" tries to mesh the interior of a given surface with
          quality and area constraints
        "" will mesh the convex hull constrained to pass through V (ignores F)
    @param[out] TV  #TV by 3 vertex position list
    @param[out] TT  #TT by 4 list of tet face indices
    @param[out] TF  #TF by 3 list of triangle face indices ('f', else
      `boundary_facets` is called on TT)
    @param[out] TR  #TT list of region ID for each tetrahedron      
    @param[out] TN  #TT by 4 list of indices neighbors for each tetrahedron ('n')
    @param[out] PT  #TV list of incident tetrahedron for a vertex ('m')
    @param[out] FT  #TF by 2 list of tetrahedrons sharing a triface ('nn')
    @param[out] num_regions Number of regions in output mesh
    """
