#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT rem:ppt3-edge-overlap (in rank_five_surface.tex, appendix_surface.tex):
for the profile-(5,5,5,5) states rho_tau (A = 1, Q = tau + sqrt(1 + tau^2)), the A|BC product vectors psi(z) of the range R_1 have
partial conjugates psi^{T_A}(z) = (|0> + conj(z)|1>) x phi(z) whose overlap with the kernel K_A(Q) of rho_tau^{T_A} vanishes
only at the endpoint product rays.

  EO-1  B q_1(z) = (|0>+z|1>) x phi(z) with phi(z) read off, and rho_tau^{T_A} K_A(Q) = 0 (K_A from verify_ppt3_surface_blocks.py).
  EO-2  k_1^+ psi^{T_A}(z) = 6((Q^2-1)|z|^2 + Q(z^2-1)) for z = x + i y; its imaginary part is 12 Q x y.
  EO-3  on the branch y = 0 (x^2 = Q/(Q^2+Q-1)) and on x = 0 (y^2 = Q/(Q^2-Q-1)), k_3^+ psi^{T_A} equals
        -6 Q x (Q^2-3)(Q^2+2Q-1)/(Q^2+Q-1), resp. -6 i Q y (Q^2-3)(Q^2-2Q-1)/(Q^2-Q-1); k_2^+ |111> = 12.
  EO-4  on sqrt3 < Q < 1+sqrt2 the factors Q^2-3, Q^2+2Q-1, Q^2+Q-1, Q^2-Q-1 are positive and Q^2-2Q-1 negative, so no ray
        has all three overlaps zero (checked by the sign of each factor at its real roots).
  EO-5  at Q = sqrt3 the four rays z = +-sqrt(2 sqrt3 - 3), +-i sqrt(3 + 2 sqrt3) have all three overlaps zero; at Q = 1+sqrt2
        the rays z = +-i do.
"""
import sympy as sp
from verify_ppt3_surface_blocks import K_poly, rho_of, pt, Q

x, y = sp.symbols('x y', real=True)
zz = sp.Symbol('zz')
s3 = sp.sqrt(3)


def ket(b):
    v = sp.zeros(8, 1)
    v[int(b, 2)] = 1
    return v


def main():
    K = K_poly()
    k1, k2, k3 = K[:, 0], K[:, 1], K[:, 2]
    # EO-1: psi(z) = (|0> + z|1>) x phi(z); u = |000> - h, s = 1, as in verify_ppt3_slocc_chirality.py (SC-1)
    om = (-1 + sp.I * s3) / 2
    omb = sp.conjugate(om)
    h = ket('011') + ket('101') + ket('110')
    bu, bv, bw = ket('000') - h, ket('111'), ket('001') + ket('010') + ket('100')
    be = ket('011') + om * ket('101') + omb * ket('110')
    bo = ket('100') + om * ket('010') + omb * ket('001')
    U, O = 2 * zz, 1 + zz**2
    W, E, V = zz * U - O, -zz * O, -zz * U - zz**2 * O
    p = (U * bu + V * bv + W * bw + E * be + O * bo).applyfunc(sp.expand)
    M = sp.Matrix(2, 4, list(p))
    assert (M[1, :] - zz * M[0, :]).applyfunc(sp.expand) == sp.zeros(1, 4)
    phi = M[0, :]
    # the 8-vector psi^{T_A} has components (1 on A=0, conj z on A=1) x phi(z)
    zc = sp.Symbol('zc')
    psiT = sp.Matrix([(1 if i == 0 else zc) * phi[j] for i in range(2) for j in range(4)])
    assert (pt(rho_of(Q), 0).applyfunc(sp.cancel) * K).applyfunc(sp.cancel) == sp.zeros(8, 3)
    print("EO-1 PASS: psi(z) = (|0>+z|1>) x phi(z); rho^{T_A} K_A = 0 identically in Q")

    def ov(k, zv, zcv):
        return sp.expand((k.H * psiT.subs({zz: zv, zc: zcv}))[0])

    zxy, zcxy = x + sp.I * y, x - sp.I * y
    # EO-2
    o1 = ov(k1, zxy, zcxy)
    claim = 6 * ((Q**2 - 1) * (x**2 + y**2) + Q * (sp.expand(zxy**2) - 1))
    assert sp.expand(o1 - claim) == 0
    assert sp.expand(sp.im(sp.expand(o1)) - 12 * Q * x * y) == 0
    print("EO-2 PASS: k_1^+ psi^T = 6((Q^2-1)|z|^2 + Q(z^2-1)), imaginary part 12 Q x y")
    # EO-3
    xc, yc = Q / (Q**2 + Q - 1), Q / (Q**2 - Q - 1)
    X, Y = sp.Symbol('X', positive=True), sp.Symbol('Y', positive=True)
    o3x = ov(k3, X, X)                                # z = X real
    o3y = ov(k3, sp.I * Y, -sp.I * Y)                 # z = i Y
    o1x, o1y = ov(k1, X, X), ov(k1, sp.I * Y, -sp.I * Y)
    assert sp.simplify(o1x.subs(X, sp.sqrt(xc))) == 0 and sp.simplify(o1y.subs(Y, sp.sqrt(yc))) == 0
    want_x = -6 * Q * X * (Q**2 - 3) * (Q**2 + 2 * Q - 1) / (Q**2 + Q - 1)
    want_y = -6 * sp.I * Q * Y * (Q**2 - 3) * (Q**2 - 2 * Q - 1) / (Q**2 - Q - 1)
    assert sp.simplify(o3x.subs(X, sp.sqrt(xc)) - want_x.subs(X, sp.sqrt(xc))) == 0
    assert sp.simplify(o3y.subs(Y, sp.sqrt(yc)) - want_y.subs(Y, sp.sqrt(yc))) == 0
    e111 = ket('111')
    assert sp.simplify((k2.H * e111)[0] - 12) == 0
    print("EO-3 PASS: k_3^+ psi^T on the two branches, and k_2^+ |111> = 12")
    # EO-4: roots of every factor, and its sign on the window
    lo, hi = s3, 1 + sp.sqrt(2)
    for name, f, sign in [("Q^2-3", Q**2 - 3, 1), ("Q^2+2Q-1", Q**2 + 2 * Q - 1, 1), ("Q^2+Q-1", Q**2 + Q - 1, 1),
                          ("Q^2-Q-1", Q**2 - Q - 1, 1), ("Q^2-2Q-1", Q**2 - 2 * Q - 1, -1)]:
        roots = [r for r in sp.solve(f, Q) if r.is_real]
        assert all(not (lo < r < hi) for r in roots), name
        mid = (lo + hi) / 2
        assert sp.sign(sp.N(f.subs(Q, mid), 30)) == sign, name
    print("EO-4 PASS: no factor has a root in (sqrt3, 1+sqrt2); signs +,+,+,+,- on the window")
    # EO-5
    for Qv, rays in [(s3, [sp.sqrt(2 * s3 - 3), -sp.sqrt(2 * s3 - 3), sp.I * sp.sqrt(3 + 2 * s3), -sp.I * sp.sqrt(3 + 2 * s3)]),
                     (1 + sp.sqrt(2), [sp.I, -sp.I])]:
        for zr in rays:
            for k in (k1, k2, k3):
                val = ov(k.subs(Q, Qv), zr, sp.conjugate(zr))
                assert sp.simplify(sp.radsimp(val)) == 0, (Qv, zr)
    print("EO-5 PASS: at Q = sqrt3 the rays z = +-sqrt(2 sqrt3 - 3), +-i sqrt(3 + 2 sqrt3), and at Q = 1+sqrt2 the rays z = +-i, are in the overlap kernel")
    print("ALL PASS")


if __name__ == '__main__':
    main()
