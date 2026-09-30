import enum
from typing import Annotated, overload

from numpy.typing import ArrayLike
import scipy.sparse


class EigsType(enum.Enum):
    EIGS_TYPE_SM = 0

    EIGS_TYPE_LM = 1

EIGS_TYPE_SM: EigsType = EigsType.EIGS_TYPE_SM

EIGS_TYPE_LM: EigsType = EigsType.EIGS_TYPE_LM

@overload
def eigs(A: scipy.sparse.csc_matrix[float], B: scipy.sparse.csc_matrix[float], k: int, type: EigsType) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]]:
    """
    Compute the first/last k eigen pairs of the generalized eigenvalue problem:

        A u = s B u

    Solutions are approximate and sorted. Only EIGS_TYPE_SM (smallest magnitude)
    is well supported.

    @param[in] A  #A by #A symmetric sparse matrix
    @param[in] B  #A by #A symmetric positive-definite sparse matrix
    @param[in] k  number of eigen pairs to compute
    @param[in] type  EigsType.EIGS_TYPE_SM or EigsType.EIGS_TYPE_LM
    @param[out] U  #A by k matrix of sorted eigenvectors (descending)
    @param[out] S  k vector of sorted eigenvalues (descending)
    @return Tuple (U, S)
    """

@overload
def eigs(A: scipy.sparse.csc_matrix[float], B: scipy.sparse.csc_matrix[float], k: int, sigma: float) -> tuple[Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')], Annotated[ArrayLike, dict(dtype='float64', shape=(None), order='C')]]:
    """
    Compute k eigen pairs of A u = s B u with a shift: A ← A + sigma * B.

    @param[in] A  #A by #A symmetric sparse matrix
    @param[in] B  #A by #A symmetric positive-definite sparse matrix
    @param[in] k  number of eigen pairs to compute
    @param[in] sigma  shift to apply to A
    @param[out] U  #A by k matrix of sorted eigenvectors (descending)
    @param[out] S  k vector of sorted eigenvalues (descending)
    @return Tuple (U, S)
    """

def lscm(V: Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), writable=False)], F: Annotated[ArrayLike, dict(dtype='int64', shape=(None, None), writable=False)]) -> Annotated[ArrayLike, dict(dtype='float64', shape=(None, None), order='C')]:
    """
    Compute a free-boundary least-squares conformal map parametrization using
    spectral decomposition. Assumes the mesh is a single connected component
    topologically equivalent to a disk.

    Implements "Spectral Conformal Parameterization" [Mullen et al. 2008].

    @param[in] V   #V by 3 list of mesh vertex positions
    @param[in] F   #F by 3 list of mesh face indices (triangles)
    @param[out] UV #V by 2 list of UV coordinates
    @return UV coordinates
    """
