#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT prop:ppt3-moduli (in rank_five_surface.tex): near the rational state
sigma of cor:ppt3-extreme5, the fully-PPT states of rank profile (5,5,5,5) form a smooth manifold of real dimension 27
(trace one), on which local filters have 18-dimensional orbits, leaving nine filter-invariant moduli.

All ranks are exact, over Q(sqrt(-3)).  A real-linear map Z -> A(Z) + B(conj Z) has real rank equal to the complex rank
of its complexification (Z, Y) -> A(Z) + B(Y), whose matrix is assembled below column by column.
  MO-1  rank-five tangent: (Z, Y) -> Z M^+ + M Y (M = [u, w, v, e, o], 8 x 5) has complex rank 55 = dim of the
        rank-five Hermitian matrices.
  MO-2  each partial transpose sigma^{T_c} has an exact three-dimensional kernel K_c.
  MO-3  the 27 differentials of the Schur-complement equations, K_c^+ (dS)^{T_c} K_c (c = A, B, C), restricted to the
        rank-five tangent, have rank 27: the profile-(5,5,5,5) locus is smooth of dimension 55 - 27 = 28 near sigma.
  MO-4  local filters: (X, Y) -> X sigma + sigma Y over local gl(2)^3 generators has complex rank 19 (the orbit of
        F sigma F^+ including scale), contained in the 28-dimensional tangent.
"""
import sys
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import verify_ppt_gme_rank5_flat_profile as F  # noqa: E402


def rank(cols):
    m = sp.Matrix.hstack(*[sp.Matrix(c) for c in cols])
    return DomainMatrix.from_Matrix(m, extension=True).rank()


def flat(m):
    return list(m)


def main():
    p = F.weights(8, 3, 2)
    N = F.numerator(p)
    u, w, v, e, o = F.parts(p['s'])
    M = sp.Matrix.hstack(u, w, v, e, o)
    units85 = []
    for i in range(8):
        for j in range(5):
            Z = sp.zeros(8, 5)
            Z[i, j] = 1
            units85.append(Z)

    tang = [Z * M.H for Z in units85] + [M * Z.T for Z in units85]     # (Z, 0) and (0, Y), Y = Z^T as a 5 x 8 unit
    r1 = rank([flat(t) for t in tang])
    print(f"MO-1 complex rank of the rank-five tangent map: {r1} (expected 55)")
    assert r1 == 55

    kernels = []
    for c in range(3):
        dm = DomainMatrix.from_Matrix(F.pt(N, c), extension=True)
        ker = dm.nullspace().to_Matrix().T
        assert ker.cols == 3
        kernels.append(sp.Matrix(ker))
        assert all(sp.simplify(x) == 0 for x in F.pt(N, c) * ker)
        print(f"MO-2 cut {c}: kernel of sigma^T_c has dimension 3")

    def constraints(dS):
        out = []
        for c, K in enumerate(kernels):
            out += flat((K.H * F.pt(dS, c) * K).applyfunc(sp.expand))
        return out
    cols = [constraints(t) for t in tang]
    r3 = rank(cols)
    print(f"MO-3 rank of the 27 Schur-complement differentials on the tangent: {r3} (expected 27); "
          f"locus dimension {r1 - r3}, trace one {r1 - r3 - 1}")
    assert r3 == 27

    I2 = sp.eye(2)
    orbit = []
    for q in range(3):
        for k in range(4):
            x = sp.zeros(2, 2)
            x[k // 2, k % 2] = 1
            ops = [I2, I2, I2]
            ops[q] = x
            X = sp.kronecker_product(*ops)
            orbit += [X * N, N * X.T]                                     # (X, 0) and (0, Y), Y = X^T
    r4 = rank([flat(t) for t in orbit])
    both = rank([flat(t) for t in tang] + [flat(t) for t in orbit])
    print(f"MO-4 complex rank of the local-filter orbit map: {r4} (expected 19); with the tangent: {both} (expected 55)")
    assert r4 == 19 and both == 55
    # orbit directions satisfy the constraints (filters preserve all four ranks)
    assert all(sp.simplify(x) == 0 for t in orbit for x in constraints(t))
    print(f"moduli: {r1 - r3} - {r4} = {r1 - r3 - r4} filter-invariant directions")


if __name__ == '__main__':
    main()
