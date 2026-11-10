#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT cor:ppt3-extreme5 and prop:ppt3-moduli on the WHOLE profile-(5,5,5,5) surface
(in rank_five_surface.tex, appendix_surface.tex).

rho_tau = [u]+V[v]+tau[w]+G(vw+wv)+E[e]+O[o] (alpha=A=1), tau=(Q^2-1)/(2Q), Q=tau+sqrt(1+tau^2) in (sqrt3, 1+sqrt2).
The Z_3 (cyclic), parity and reflection symmetries of rho_tau split both face-direction systems into small blocks that
are evaluated at the single cut A|BC.  Each block has a square minor that is an explicit polynomial in Q; we compute it,
factor it, and check that no factor has a root inside the window.
  SB-1  rho_tau^{T_A} K_A(Q) = 0 for an explicit degree-3 polynomial K_A(Q) with independent columns (kernel dimension 3
        since rank rho^{T_A} = 5 on the window by thm:ppt3-orbit).
  SB-2  extremality: the 25 real parameters of Hermitian X_D split into six block types; the block minors are
        nonzero on the window, so the face-direction system has rank 24 (nullspace = rho_tau).
  SB-3  moduli: the 27 Schur-complement differentials restricted to the 55-dim rank-5 tangent split into sectors
        (charge q=0,1 and parity); each sector map is onto, via an explicit 5x5 or 4x4 minor, so rank 27 for q=0,1,2.
  SB-4  independent exact cross-check at Q=2 (tau=3/4) WITHOUT the block reduction: complexified face system over all
        three cuts has rank 24, the 27 constraints rank 27 on the 55-dim tangent, the filter orbit rank 19.
"""
import sympy as sp
from sympy.polys.matrices import DomainMatrix

Q = sp.Symbol('Q', positive=True)
s3 = sp.sqrt(3)
om = (-1 + sp.I * s3) / 2
omb = sp.conjugate(om)
I = sp.I


def ket(b):
    v = sp.zeros(8, 1)
    v[int(b, 2)] = 1
    return v


u = ket('000') - ket('011') - ket('101') - ket('110')
v = ket('111')
w = ket('001') + ket('010') + ket('100')
e = ket('011') + om * ket('101') + omb * ket('110')
o = ket('100') + om * ket('010') + omb * ket('001')
# orthogonal complement of R: rows of K in appendix_family.tex (proof of prop:ppt3-witness)
up = 3 * ket('000') + ket('011') + ket('101') + ket('110')
eb, ob = e.applyfunc(sp.conjugate), o.applyfunc(sp.conjugate)
VEC = [u, v, w, e, o, up, eb, ob]
CHARGE = [0, 0, 0, 2, 2, 0, 1, 1]          # cyclic shift acts on e,o,... by omega^charge (e, o: omega^2)
PAR = [1, -1, -1, 1, -1, 1, 1, -1]         # global parity Z^{x3}


def pt(M, c):
    out = sp.zeros(8, 8)
    for a in range(8):
        for b in range(8):
            bits_a = [(a >> k) & 1 for k in (2, 1, 0)]
            bits_b = [(b >> k) & 1 for k in (2, 1, 0)]
            ta, tb = bits_a[:], bits_b[:]
            ta[c], tb[c] = bits_b[c], bits_a[c]
            out[int(''.join(map(str, ta)), 2), int(''.join(map(str, tb)), 2)] = M[a, b]
    return out


def pair(i, j):
    return VEC[i] * VEC[j].H


def PT(i, j, c=0):
    return pt(pair(i, j), c)


def weights(Qv):
    tau = (Qv**2 - 1) / (2 * Qv)
    kap = 1 - tau**2
    E = kap * Qv**2
    O = kap / tau
    G = (E - 1) / tau
    V = (E**2 - 2 * kap * E + kap) / (tau * kap)
    return tau, V, G, E, O


def Xmat(Qv):
    tau, V, G, E, O = weights(Qv)
    return sp.Matrix([[1, 0, 0, 0, 0], [0, V, G, 0, 0], [0, G, tau, 0, 0], [0, 0, 0, E, 0], [0, 0, 0, 0, O]])


def rho_of(Qv):
    X = Xmat(Qv)
    return sum((X[i, j] * pair(i, j) for i in range(5) for j in range(5)), sp.zeros(8, 8))


def rk(m):
    return DomainMatrix.from_Matrix(sp.Matrix(m), extension=True).rank()


def pivots(m):
    return DomainMatrix.from_Matrix(sp.Matrix(m), extension=True).rref()[1]


def K_poly():
    r3 = s3 * I
    k0 = sp.zeros(8, 1)
    k0[1], k0[2], k0[4] = (3 - r3) * Q, (3 + r3) * Q, 3 * (Q**2 - 1)
    k1 = sp.zeros(8, 1)
    k1[1] = Q * (-(3 + r3) * Q**2 + 15 - 3 * r3)
    k1[2] = Q * (-(3 - r3) * Q**2 + 15 + 3 * r3)
    k1[7] = 12
    k2 = sp.zeros(8, 1)
    k2[0], k2[3] = -3 * Q * (Q**2 - 1), -6 * Q
    k2[5] = (1 - omb) * Q**2 + 3 * (1 + omb)
    k2[6] = (1 - om) * Q**2 + 3 * (1 + om)
    return sp.Matrix.hstack(k0, k1, k2)


def check_window_free(expr, name):
    """expr: polynomial in Q with coefficients in Q(sqrt-3).  No root of |expr|^2 in the open window (sqrt3, 1+sqrt2)."""
    expr = sp.expand(sp.simplify(expr))
    P = sp.expand(sp.simplify(expr * sp.conjugate(expr)))
    P = sp.Poly(P, Q)
    assert all(c.is_rational for c in P.all_coeffs()), P
    lo, hi = sp.sqrt(3).evalf(60), (1 + sp.sqrt(2)).evalf(60)
    bad, boundary = [], []
    for f, _ in sp.factor_list(P.as_expr())[1]:
        if sp.Poly(f, Q).degree() == 0:
            continue
        for r in sp.Poly(f, Q).real_roots():
            rv = sp.N(r, 60)
            if abs(rv - lo) < 1e-40 or abs(rv - hi) < 1e-40:
                boundary.append(f)
            elif lo < rv < hi:
                bad.append((f, rv))
    assert not bad, (name, bad)
    fac = sp.factor(sp.simplify(expr))
    print(f"   {name}: minor = {fac}")
    return fac


def sector_minor(cols_mats, K, nrows_sel=None):
    """Columns: images M*K (or K^+ M K entries); return (rank at Q=2, minor polynomial of a maximal square submatrix)."""
    M = sp.Matrix.hstack(*cols_mats)
    M2 = M.subs(Q, 2)
    r = rk(M2)
    cpiv = list(pivots(M2))
    sub = M.extract(list(range(M.rows)), cpiv)
    rpiv = list(pivots(sub.subs(Q, 2).T))
    return r, sub.extract(rpiv, list(range(len(cpiv)))).det(method='berkowitz'), cpiv


def main():
    # ---------- SB-1
    K = K_poly()
    rho = rho_of(Q)
    rhoTA = pt(rho, 0).applyfunc(sp.cancel)
    assert (rhoTA * K).applyfunc(sp.cancel) == sp.zeros(8, 3)
    assert sp.expand(K[4, 0]) != 0 and K[7, 1] == 12 and K[7, 0] == 0 and K[4, 1] == 0 and K[7, 2] == 0 and K[4, 2] == 0
    assert K[3, 2] == -6 * Q and K[3, 0] == 0 and K[3, 1] == 0
    print("SB-1 PASS: rho_tau^{T_A} K_A(Q) = 0 identically in Q; columns independent for Q>1 (pivots 3(Q^2-1), 12, -6Q)")
    assert rk(rhoTA.subs(Q, 2)) == 5

    # ---------- SB-2: extremality blocks
    def Ei(i, j):
        return PT(i, j) * K

    groups = [
        ('q=1, parity +   [ue, vo, wo]', [Ei(0, 3), Ei(1, 4), Ei(2, 4)], 3),
        ('q=1, parity -   [uo, ve, we]', [Ei(0, 4), Ei(1, 3), Ei(2, 3)], 3),
        ('q=0, parity -, sym  [uv, uw, eo]', [Ei(0, 1) + Ei(1, 0), Ei(0, 2) + Ei(2, 0), Ei(3, 4) + Ei(4, 3)], 3),
        ('q=0, parity -, anti [uv, uw, eo]', [I * (Ei(0, 1) - Ei(1, 0)), I * (Ei(0, 2) - Ei(2, 0)), I * (Ei(3, 4) - Ei(4, 3))], 3),
        ('q=0, parity +, anti [vw]', [I * (Ei(1, 2) - Ei(2, 1))], 1),
        ('q=0, parity +, sym  [uu, vv, vw, ww, ee, oo]', [Ei(0, 0), Ei(1, 1), Ei(1, 2) + Ei(2, 1), Ei(2, 2), Ei(3, 3), Ei(4, 4)], 5),
    ]
    total = 2 * (3 + 3)                                   # q=2 is the Hermitian conjugate of q=1
    nunk = 2 * (3 + 3) + 3 + 3 + 1 + 6
    assert nunk == 25
    for name, mats, expect in groups:
        cols = [sp.expand(m).reshape(24, 1) for m in mats]
        r, minor, _ = sector_minor(cols, K)
        assert r == expect, (name, r, expect)
        check_window_free(minor, f"SB-2 [{name}] rank {r}")
        if not name.startswith('q=1'):
            total += r
        else:
            pass
    # the rho_tau direction lies in the last group's kernel, so the total rank is 24 and the face is {rho_tau}
    assert total == 24, total
    # direction rho_tau itself is annihilated (it is a face direction)
    rho_dir = sum((Xmat(Q)[i, j] * Ei(i, j) for i in range(5) for j in range(5)), sp.zeros(8, 3))
    assert rho_dir.applyfunc(sp.cancel) == sp.zeros(8, 3)
    print(f"SB-2 PASS: face-direction system has rank {total} = 25-1 for every Q in the window: every rho_tau,a is extreme")

    # ---------- SB-3: moduli blocks (tangent: pairs (i,j), min(i,j)<5 in the basis u,v,w,e,o,u',eb,ob)
    KH = K.H
    rows_plus = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 2)]
    rows_minus = [(0, 2), (1, 2), (2, 0), (2, 1)]
    ntan = 0
    total_rank = 0
    for q in (0, 1, 2):
        for p, rows in ((1, rows_plus), (-1, rows_minus)):
            idx = [(i, j) for i in range(8) for j in range(8)
                   if not (i >= 5 and j >= 5) and (CHARGE[i] - CHARGE[j]) % 3 == q and PAR[i] * PAR[j] == p]
            ntan += len(idx)
            if q == 2:                                    # q=2 is the Hermitian conjugate of q=1
                total_rank += len(rows)
                continue
            cols = [sp.Matrix([sp.expand((KH * PT(i, j) * K)[r, c]) for (r, c) in rows]) for (i, j) in idx]
            r, minor, _ = sector_minor(cols, K)
            assert r == len(rows), (q, p, r)
            check_window_free(minor, f"SB-3 [q={q}, parity {p:+d}] {len(idx)} tangent directions, onto rank {r}")
            total_rank += r
    assert ntan == 55 and total_rank == 27, (ntan, total_rank)
    print("SB-3 PASS: 55-dim tangent, 27 independent Schur-complement differentials for every Q in the window "
          "(manifold dimension 28, trace one 27, minus 18 filter dimensions = 9 moduli)")

    # ---------- SB-4: exact cross-check at Q=2 without the block reduction
    Q0 = sp.Integer(2)
    rho0 = rho_of(Q0)
    ker = []
    for c in range(3):
        kc = DomainMatrix.from_Matrix(pt(rho0, c), extension=True).nullspace().to_Matrix().T
        assert kc.cols == 3
        ker.append(kc)
    cols = []
    for i in range(5):                                    # complexified face system: D=B X B^+, both D^T K=0 and K^+ D^T=0
        for j in range(5):
            D = pair(i, j)
            col = []
            for c in range(3):
                col += list(pt(D, c) * ker[c]) + list(ker[c].H * pt(D, c))
            cols.append(col)
    assert rk(sp.Matrix(cols).T) == 24
    mats = [pair(i, j) for i in range(8) for j in range(8) if not (i >= 5 and j >= 5)]
    assert len(mats) == 55
    cons = [sum([list((ker[c].H * pt(D, c) * ker[c])) for c in range(3)], []) for D in mats]
    assert rk(sp.Matrix(cons).T) == 27
    I2 = sp.eye(2)
    orb = []
    for q in range(3):
        for k in range(4):
            x = sp.zeros(2, 2)
            x[k // 2, k % 2] = 1
            ops = [I2, I2, I2]
            ops[q] = x
            Xl = sp.kronecker_product(*ops)
            orb += [list(Xl * rho0), list(rho0 * Xl.T)]
    assert rk(sp.Matrix(orb).T) == 19
    print("SB-4 PASS: at Q=2 (tau=3/4, no block reduction) face-system rank 24, Schur-differential rank 27 on the 55-dim "
          "tangent, filter-orbit rank 19")
    print("ALL PASS")


if __name__ == '__main__':
    main()
