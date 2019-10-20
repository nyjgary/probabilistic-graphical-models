# relevant functions copied from current version (0.7.7) of pyldpc package, to supplement our usage of its older version

import math
import numpy as np
import scipy
from scipy.stats import norm
from scipy.sparse import csr_matrix
pi = math.pi


# def int2bitarray(N, k):
#     """
#     Changes array's base from int (base 10) to binary (base 2)
#     Parameters:
#     ===========
#     N: int N
#     k: Width of the binary array you would like to change N into.
#     N must not be greater than 2^k - 1.
#     >> Examples: int2bitarray(6,3) returns [1, 1, 0]
#                  int2bitarray(6,5) returns [0, 0, 1, 1,0]
#                  int2bitarray(255,8) returns [1, 1, 1, 1, 1, 1, 1, 1]
#                  int2bitarray(255,10) returns [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
#     """

#     binary_string = bin(N)
#     length = len(binary_string)
#     bitarray = np.zeros(k, 'int')
#     for i in range(length-2):
#         bitarray[k-i-1] = int(binary_string[length-i-1])

#     return bitarray


# def bitarray2int(bitarray):

#     """ Changes array's base from binary (base 2) to int (base 10).
#     Parameters:
#     ===========
#     bitarray: Binary Array.
#     >> Examples: bitarray2int([1, 1, 0]) returns 6
#                  bitarray2int([0, 0, 1, 1,0]) returns 6
#                  bitarray2int([1, 1, 1, 1, 1, 1, 1, 1]) returns 255
#     """

#     bitstring = "".join([str(i) for i in bitarray])

#     return int(bitstring, 2)


def binaryproduct(X, Y):

    """ Binary Matrices or Matrix-vector product in Z/2Z.
    Works with scipy.sparse.csr_matrix matrices X, Y too."""

    A = X.dot(Y)

    if type(A) != scipy.sparse.csr_matrix:
        return A % 2

    return A.toarray() % 2


def gaussjordan(X, change=0):

    """
    Description:
    Performs the row reduced echelon form of X and returns it.
    If change = 1, all changes in the X's rows are applied to
    identity matrix P:
    Let A be our parameter X. refA the reduced echelon form of A.
    P is the square invertible matrix:
    P.A = Aref.
    -------------------------------------------------
    Parameters:
    X: 2D-Array.
    change : boolean (default = 0)
    ------------------------------------------------
    change = 0  (default)
     >>> Returns 2D-Array Row Reduced Echelon form of Matrix
    change = 1
    >>> Returns Tuple of 2D-arrays (refX, P) where P is described above.
    """

    A = np.copy(X)
    m, n = A.shape

    if change:
        P = np.identity(m).astype(int)

    pivot_old = -1
    for j in range(n):
        filtre_down = A[pivot_old+1:m, j]
        pivot = np.argmax(filtre_down)+pivot_old+1

        if A[pivot, j]:
            pivot_old += 1
            if pivot_old != pivot:
                aux = np.copy(A[pivot, :])
                A[pivot, :] = A[pivot_old, :]
                A[pivot_old, :] = aux
                if change:
                    aux = np.copy(P[pivot, :])
                    P[pivot, :] = P[pivot_old, :]
                    P[pivot_old, :] = aux

            for i in range(m):
                if i != pivot_old and A[i, j]:
                    if change:
                        P[i, :] = abs(P[i, :]-P[pivot_old, :])
                    A[i, :] = abs(A[i, :]-A[pivot_old, :])

        if pivot_old == m-1:
            break

    if change:
        return A, P
    return A


# def binaryrank(X):
#     """ Computes rank of a binary Matrix using Gauss-Jordan algorithm"""
#     A = np.copy(X)
#     m, n = A.shape

#     A = gaussjordan(A)

#     return sum([a.any() for a in A])


# def f1(y, sigma):
#     """ Normal Density N(1,sigma) """
#     f = norm.pdf(y, loc=1, scale=sigma)
#     return f


# def fm1(y, sigma):
#     """ Normal Density N(-1,sigma) """

#     f = norm.pdf(y, loc=-1, scale=sigma)
#     return f


# def bits2i(H, i):
#     """
#     Computes list of elements of N(i)-j:
#     List of variables (bits) connected to Parity node i.
#     """
#     if type(H) != scipy.sparse.csr_matrix:
#         m, n = H.shape
#         return list(np.where(H[i])[0])

#     indj = H.indptr
#     indi = H.indices

#     return [indi[a] for a in range(indj[i], indj[i+1])]


# def nodes2j(H, j):

#     """
#     Computes list of elements of M(j):
#     List of nodes (PC equations) connecting variable j.
#     """

#     return bits2i(H.T, j)


# def bitsandnodes(H):

#     m, n = H.shape

#     bits = [bits2i(H, i) for i in range(m)]
#     nodes = [nodes2j(H, j)for j in range(n)]

#     return bits, nodes


# def incode(H, x):

#     """ Computes Binary Product of H and x. If product is null, x is in the code.
#         Returns appartenance boolean.
#     """

#     return (binaryproduct(H, x) == 0).all()


# def gausselimination(A, b):

#     """ Applies Gauss Elimination Algorithm to X in order to solve a
#     linear system X.X = B. X is transformed to row echelon form:
#          |1 * * * * * |
#          |0 1 * * * * |
#          |0 0 1 * * * |
#          |0 0 0 1 * * |
#          |0 0 0 0 1 * |
#          |0 0 0 0 0 1 |
#          |0 0 0 0 0 0 |
#          |0 0 0 0 0 0 |
#          |0 0 0 0 0 0 |
#     Same row operations are applied on 1-D Array vector B.
#     Both arguments are sent back.
#     --------------------------------------
#     Parameters:
#     X: 2D-array.
#     B:      1D-array. Size must equal number of rows of X.
#     -----------------------------------
#     Returns:
#     Modified arguments X, B as described above.
#          """
#     if type(A) == scipy.sparse.csr_matrix:
#         A = A.toarray().copy()
#     else:
#         A = A.copy()
#     b = b.copy()
#     n, k = A.shape

#     for j in range(min(k, n)):
#         listedepivots = [i for i in range(j, n) if A[i, j]]
#         if len(listedepivots):
#             pivot = np.min(listedepivots)
#         else:
#             continue
#         if pivot != j:
#             aux = (A[j, :]).copy()
#             A[j, :] = A[pivot, :]
#             A[pivot, :] = aux

#             aux = b[j].copy()
#             b[j] = b[pivot]
#             b[pivot] = aux

#         for i in range(j+1, n):
#             if A[i, j]:
#                 A[i, :] = abs(A[i, :]-A[j, :])
#                 b[i] = abs(b[i]-b[j])

#     return A, b

def coding_matrix_systematic(X, sparse=True):

    """
    Description:
    Solves H.G' = 0 and finds the coding matrix G in the systematic form :
    [I_k  A] by applying permutations on X.
    CAUTION: RETURNS TUPLE (Hp,tGS) WHERE Hp IS A MODIFIED VERSION OF THE
    GIVEN PARITY CHECK X, tGS THE TRANSPOSED
    SYSTEMATIC CODING X ASSOCIATED TO Hp. YOU MUST USE THE RETURNED TUPLE
    IN CODING AND DECODING, RATHER THAN THE UNCHANGED
    PARITY-CHECK X H.
    -------------------------------------------------
    Parameters:
    X: 2D-Array. Parity-check matrix.
    sparse: (optional, default True): use scipy.sparse matrices
    to speed up computation if n>100.
    ------------------------------------------------
    >>> Returns Tuple of 2D-arrays (Hp,GS):
        Hp: Modified H: permutation of columns (The code doesn't change)
        tGS: Transposed Systematic Coding matrix associated to Hp.
    """

    H = X.copy()
    m, n = H.shape

    if n > 100 and sparse:
        sparse = True
    else:
        sparse = False

    P1 = np.identity(n, dtype=int)

    Hrowreduced = gaussjordan(H)

    k = n - sum([a.any() for a in Hrowreduced])

    # After this loop, Hrowreduced will have the form H_ss : | I_(n-k)  A |

    while(True):
        zeros = [i for i in range(min(m, n)) if not Hrowreduced[i, i]]
        indice_colonne_a = min(zeros)
        list_ones = [j for j in range(indice_colonne_a+1, n)
                     if Hrowreduced[indice_colonne_a, j]]
        if not len(list_ones):
            break

        indice_colonne_b = min(list_ones)

        aux = Hrowreduced[:, indice_colonne_a].copy()
        Hrowreduced[:, indice_colonne_a] = Hrowreduced[:, indice_colonne_b]
        Hrowreduced[:, indice_colonne_b] = aux

        aux = P1[:, indice_colonne_a].copy()
        P1[:, indice_colonne_a] = P1[:, indice_colonne_b]
        P1[:, indice_colonne_b] = aux

    # NOW, Hrowreduced has the form: | I_(n-k)  A | ,
    # the permutation above makes it look like :
    # |A  I_(n-k)|

    P1 = P1.T
    identity = list(range(n))
    sigma = identity[n-k:] + identity[:n-k]

    P2 = np.zeros(shape=(n, n), dtype=int)
    P2[identity, sigma] = np.ones(n)

    if sparse:
        P1 = csr_matrix(P1)
        P2 = csr_matrix(P2)
        H = csr_matrix(H)

    P = binaryproduct(P2, P1)

    if sparse:
        P = csr_matrix(P)

    Hp = binaryproduct(H, np.transpose(P))

    GS = np.zeros((k, n), dtype=int)
    GS[:, :k] = np.identity(k)
    GS[:, k:] = (Hrowreduced[:n-k, n-k:]).T

    return Hp, GS.T