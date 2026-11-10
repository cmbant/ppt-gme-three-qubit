#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT thm:ppt3-orbit, lem:ppt3-h0 and cor:ppt3-extreme5 (analytic proof of GME; in fully_ppt_three_qubit_gme.tex):
an analytic supported witness for the two-parameter surface rho(m,n,a) of thm:ppt3-orbit, with no solver data.

Notation: B_s=[u_s,v,w,e,o], u_s=|000>-s(|011>+|101>+|110>), q_s(z)=(2z,-z^2(z^2+3s),z^2-s,-z(z^2+s),z^2+s) (coordinates
of the A|BC product rays), r=|z|^2, z=x+iy.  The supported witness is the diagonal form H_s=diag(-12 s^3, 1, 11 s^2, -6 s, 14 s^2) on
(u_s,v,w,e,o), with H_5:=H_1.  Parts: (A) completeness of the product rays; (B) the identity q_s^+ H_s q_s=(r-s)^2[(r-5s)^2+12 s x^2],
the gauge reduction to s=1 (D_s=diag(1,s^(3/2),s^(1/2),s,s^(1/2))), the expansion used in the paper's proof and cyclic invariance of H_s;
(C) sigma of Cor 4b: H_(1/2) and tr(H X_sigma)=-493577/32646438; (D) the surface is one orbit of local filters diag(1,lam)^(x3), so everything depends on tau=W/alpha only;
(E) exact PSD/PPT window 1/sqrt3<tau<1 (Gram factor and leading minors in Q=tau+sqrt(1+tau^2)) and (F) GME window
3/4<=tau<1 (tr(H_5 X_t)=-(Q^2-2Q-1)S_6(Q)/(2Q(Q^2-1)), shifted S_6 and its sign bound); (G) the tau-independent
form H0 with coefficients in Q(sqrt3) whose two axis restrictions are the squares (Y-1)^2(Y-y1)^2, (X-x0)^2(X^2+pX+q), which is
negative on rho_tau for 1/sqrt3<tau<=3/4 and closes the window; (H) the upper endpoint is an equal mixture of two
fully product states, with rank profile (2,2,2,2), at every finite local-filter parameter a>0.
"""
import sys
from pathlib import Path

import sympy as sp
from sympy import Matrix, Rational as R, symbols, sqrt, I

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import certify_ppt_gme_dps2_sym as C  # noqa: E402
import verify_ppt_gme_rank5_flat_profile as V4  # noqa: E402

OM = V4.OMEGA
x, y = symbols('x y', real=True)
s = symbols('s', positive=True)
z = symbols('z')


def qvec(zz, ss):
    return Matrix([2 * zz, -zz**2 * (zz**2 + 3 * ss), zz**2 - ss, -zz * (zz**2 + ss), zz**2 + ss])


def form(ss):
    """Supported witness H_s = diag(-12 s^3, 1, 11 s^2, -6 s, 14 s^2) in the basis (u_s,v,w,e,o)."""
    ss = sp.sympify(ss)
    return sp.diag(-12 * ss**3, 1, 11 * ss**2, -6 * ss, 14 * ss**2)


def bmat(ss):
    u, w, v, e, o = V4.parts(ss)
    return u.row_join(v).row_join(w).row_join(e).row_join(o)


def section_a():
    b = bmat(s)
    ray = (b * qvec(z, s)).applyfunc(sp.expand)
    assert all(sp.expand(ray[4 + j] - z * ray[j]) == 0 for j in range(4))        # second row = z * first row
    # the 4x4 minors of the pencil row1 - z row0 in the unknown coordinates (u,v,w,e,o)
    cu, cv, cw, ce, co = symbols('cu cv cw ce co')
    r0 = [cu, cw + OM**2 * co, cw + OM * co, -s * cu + ce]
    r1 = [cw + co, -s * cu + OM * ce, -s * cu + OM**2 * ce, cv]
    cs = (cu, cv, cw, ce, co)
    pencil = Matrix([[sp.diff(r1[j] - z * r0[j], c_) for c_ in cs] for j in range(4)])
    expected = [s + z**2, z * (s + z**2), z**2 - s, z**2 * (3 * s + z**2), 2 * z]
    for k, ex in enumerate(expected):
        minor = pencil.extract(range(4), [i for i in range(5) if i != 4 - k]).det()
        ratio = sp.simplify(minor / (2 * OM + 1))
        assert sp.simplify(ratio - ex) == 0 or sp.simplify(ratio + ex) == 0
    # no common finite root (s>0) and the only degree-4 minor has leading coefficient 1 (point at infinity is v)
    assert sp.resultant(z**2 + s, z**2 - s, z) == 4 * s**2
    assert sp.degree(expected[3], z) == 4 and all(sp.degree(e_, z) <= 3 for k, e_ in enumerate(expected) if k != 3)
    print('A: PASS cut-A product rays of span(u_s,v,w,e,o) are exactly B_s q_s(z), z in C, and v (all s>0)')


def section_b():
    H = form(s)
    zz = x + sp.I * y
    r = x**2 + y**2
    val = sp.expand((qvec(zz, s).H * H * qvec(zz, s))[0])
    assert sp.expand(val - (r - s)**2 * ((r - 5 * s)**2 + 12 * s * x**2)) == 0
    h5 = form(1)
    assert h5 == sp.diag(-12, 1, 11, -6, 14)
    val1 = sp.expand((qvec(zz, 1).H * h5 * qvec(zz, 1))[0])
    assert sp.expand(val1 - (r - 1)**2 * ((r - 5)**2 + 12 * x**2)) == 0
    # the expansion in the paper's proof: d = Re z^2, squared moduli of the entries of q_1, then collection
    d, rr = symbols('d rr', real=True)
    mods = [4 * rr, rr**2 * (rr**2 + 9 + 6 * d), rr**2 + 1 - 2 * d, rr * (rr**2 + 1 + 2 * d), rr**2 + 1 + 2 * d]
    for m, c_ in zip(mods, list(qvec(zz, 1))):
        assert sp.expand(m.subs({rr: r, d: x**2 - y**2}) - sp.expand(c_ * sp.conjugate(c_))) == 0
    total = -12 * mods[0] + 1 * mods[1] + 11 * mods[2] - 6 * mods[3] + 14 * mods[4]
    assert sp.expand(total - (rr**4 - 6 * rr**3 + 34 * rr**2 - 54 * rr + 25 + 6 * d * (rr - 1)**2)) == 0
    assert sp.expand((rr - 1)**2 * (rr**2 - 4 * rr + 25) - (rr**4 - 6 * rr**3 + 34 * rr**2 - 54 * rr + 25)) == 0
    assert sp.expand((rr - 5)**2 + 12 * ((rr + d) / 2) - (rr**2 - 4 * rr + 25 + 6 * d)) == 0   # x^2 = (r+d)/2
    # gauge: q_s(sqrt(s) z) = sqrt(s) D_s q_1(z), H_s = s^3 D_s^-1 H_5 D_s^-1
    D = sp.diag(1, s**R(3, 2), sp.sqrt(s), s, sp.sqrt(s))
    zq = symbols('zq')
    assert (qvec(sp.sqrt(s) * zq, s) - sp.sqrt(s) * D * qvec(zq, 1)).applyfunc(sp.simplify) == sp.zeros(5, 1)
    assert (H - s**3 * D.inv() * h5 * D.inv()).applyfunc(sp.simplify) == sp.zeros(5, 5)
    # cyclic permutation i -> CYCLE[i] acts on B_s by a diagonal matrix of phases; H_s is diagonal, hence invariant
    perm = sp.zeros(8, 8)
    for i in range(8):
        perm[C.CYCLE[i], i] = 1
    b = bmat(s)
    dd = sp.diag(*[sp.simplify((b.col(j).H * perm * b.col(j))[0] / (b.col(j).H * b.col(j))[0]) for j in range(5)])
    assert all(sp.simplify(t) == 0 for t in (perm * b - b * dd).applyfunc(sp.simplify))
    assert all(sp.simplify(sp.Abs(dd[j, j]) - 1) == 0 for j in range(5))
    assert all(sp.simplify(t) == 0 for t in (dd.H * H * dd - H))
    print('B: PASS q_s^+ H_s q_s = (r-s)^2[(r-5s)^2 + 12 s x^2] >= 0 (gauge reduction to s=1 checked; H_5 = diag(-12,1,11,-6,14), '
          'expansion r^4-6r^3+34r^2-54r+25+6d(r-1)^2); H_s diagonal, invariant under the cyclic phases')


def weights_x(p):
    return Matrix([[p['A'], 0, 0, 0, 0], [0, p['V'], p['G'], 0, 0], [0, p['G'], p['W'], 0, 0],
                   [0, 0, 0, p['E'], 0], [0, 0, 0, 0, p['O']]])


def section_c():
    p = V4.weights(8, 3, 2)
    h = form(R(1, 2))
    assert h == sp.diag(R(-3, 2), 1, R(11, 4), -3, R(7, 2))
    num = V4.numerator(p)
    b = bmat(R(1, 2))
    xs = weights_x(p)
    assert (num - b * xs * b.H).applyfunc(sp.simplify) == sp.zeros(8, 8)
    # X_sigma in the basis (u_{1/2}, v, w, e, o): 16323219 X_sigma = 30000 X, diagonal (4*825000, 2021184, 1440000, 951720, 450625), vw 145200
    nx = (30000 * xs).applyfunc(sp.nsimplify)
    expect = sp.diag(4 * 825000, 2021184, 1440000, 951720, 450625).as_mutable()
    expect[1, 2] = expect[2, 1] = 145200
    assert nx == expect and 30000 * sp.expand(num.trace()) == 16323219
    value = sp.trace(h * xs) / sp.expand(num.trace())
    assert value == R(-493577, 32646438) and value < 0
    print(f'C: PASS Cor 4b: s=1/2, H_(1/2)=diag(-3/2,1,11/4,-3,7/2), 16323219 X_sigma = diag(4*825000,2021184,1440000,951720,450625)+145200(vw), '
          f'tr(H sigma)={value} = {sp.N(value, 8)} < 0; q^+Hq=(r-1/2)^2[(r-5/2)^2+6x^2]')


def section_d():
    m, n, a = symbols('m n a', positive=True)
    p1, pa = V4.weights(m, n, 1), V4.weights(m, n, a)
    lam2 = 1 / a                                   # F = diag(1,lam)^(x3): u_s->u_{s lam^2}, w->lam w, v->lam^3 v, e->lam^2 e, o->lam o
    scal = dict(A=1, W=lam2, V=lam2**3, G=lam2**2, E=lam2**2, O=lam2, s=lam2)
    assert all(sp.simplify(pa[k] - a * scal[k] * p1[k]) == 0 for k in ('A', 'W', 'V', 'G', 'E', 'O'))
    assert sp.simplify(pa['s'] - lam2 * p1['s']) == 0
    print('D: PASS rho(m,n,a) = a F_lam rho(m,n,1) F_lam^+ with lam^2=1/a: the surface is one local-filter orbit, tau=W/alpha only')


Qv = symbols('Q', positive=True)


def surface_q():
    """Weights of rho_t at alpha = A = 1 in the coordinate Q = t + sqrt(1+t^2), t = (Q^2-1)/(2Q) (paper, proof of thm:ppt3-orbit)."""
    t = (Qv**2 - 1) / (2 * Qv)
    K = 1 - t**2
    E = K * Qv**2
    p = dict(A=1, W=t, O=K / t, E=E, G=(E - 1) / t, V=(E**2 - 2 * K * E + K) / (t * K), s=1)
    return t, K, {k: sp.cancel(v) for k, v in p.items()}


def vanishes(expr):
    return sp.cancel(sp.expand(sp.together(sp.expand(expr)))) == 0


def section_e():
    """PPT window: tV-G^2 = 2Q^2K(Q^2-3); odd PT block = Gram [[L,M],[M^+,M^+L^-1M]]; even PT block det 0, minors 1, E, m3."""
    t, K, p = surface_q()
    assert vanishes(t * p['V'] - p['G']**2 - 2 * Qv**2 * K * (Qv**2 - 3))
    u, w, v, e, o = V4.parts(R(1))
    rho = (u * u.H + p['W'] * w * w.H + p['V'] * v * v.H + p['G'] * (w * v.H + v * w.H)
           + p['E'] * e * e.H + p['O'] * o * o.H)
    pta = V4.pt(rho, 0)
    odd, even = pta.extract(V4.ODD, V4.ODD), pta.extract(V4.EVEN, V4.EVEN)
    L = Matrix([[t + p['O'], t + OM * p['O']], [t + OM**2 * p['O'], t + p['O']]])
    M = Matrix([[-1, 1 + OM * p['E']], [-1, 1 + OM**2 * p['E']]])
    assert vanishes(L.det() - 3 * K)
    gram = Matrix(sp.BlockMatrix([[L, M], [M.H, M.H * L.inv() * M]]))
    assert all(vanishes(c_) for c_ in odd - gram)
    assert all(vanishes(c_) for c_ in (M.H * L.inv() * M - Matrix([[t + p['O'], p['G']], [p['G'], p['V']]])))
    b_ = t + OM * p['O']
    shown = Matrix([[1, -1, b_, sp.conjugate(b_)], [-1, 1 + p['E'], p['G'], p['G']],
                    [sp.conjugate(b_), p['G'], 1 + p['E'], 1 + OM**2 * p['E']], [b_, p['G'], 1 + OM * p['E'], 1 + p['E']]])
    assert all(vanishes(c_) for c_ in even - shown)                       # the displayed even block
    assert vanishes(even.det())
    assert vanishes(p['E']**2 - 2 * K * (1 + 2 * t**2) * p['E'] + K**2)    # the bracket of the displayed det factorisation
    minors = [even.extract(range(k), range(k)).det() for k in (1, 2, 3)]
    expected = [1, p['E'], (Qv**2 - 3) * (Qv**4 + 3) * K**2 / (4 * t**2)]
    assert all(vanishes(a_ - b_) for a_, b_ in zip(minors, expected))
    assert sp.solve(Qv**2 - 3, Qv) == [sqrt(3)] and vanishes(t.subs(Qv, sqrt(3)) - 1 / sqrt(3))
    print('E: PASS exact in Q = tau+sqrt(1+tau^2): tV-G^2 = 2Q^2K(Q^2-3); odd PT_A block is the Gram matrix of (L, M) with '
          'det L = 3K > 0 (rank 2); even block has det 0 and leading minors 1, E, (Q^2-3)(Q^4+3)K^2/(4t^2) (rank 3); '
          'so the PSD/PPT window is Q^2 >= 3, i.e. tau >= 1/sqrt3, profile (5,5,5,5) inside')


L5 = Qv**5 + 4 * Qv**4 + 2 * Qv**3 - 8 * Qv**2 - 27 * Qv + 4
L10 = (Qv**10 + 8 * Qv**9 + 20 * Qv**8 - 16 * Qv**7 - 98 * Qv**6 - 96 * Qv**5 - 44 * Qv**4 + 272 * Qv**3
       + 233 * Qv**2 - 168 * Qv + 16)


def section_f():
    """GME for 3/4 <= tau < 1 (Q in [2, 1+sqrt2)): tr(H_5 X_t) = -(Q^2-2Q-1) S_6(Q) / (2Q(Q^2-1)) and S_6(2+x) < 0 on [0, sqrt2-1]."""
    t, K, p = surface_q()
    f0 = sp.cancel((form(1) * weights_x(p)).trace())
    r2 = Qv**2 - 2 * Qv - 1
    S6 = Qv**6 - Qv**5 - 9 * Qv**4 + 2 * Qv**3 + 11 * Qv**2 + 27 * Qv - 3
    assert vanishes(f0 + r2 * S6 / (2 * Qv * (Qv**2 - 1)))
    xx = symbols('xx')
    shifted = sp.Poly(sp.expand(S6.subs(Qv, 2 + xx)), xx).all_coeffs()[::-1]
    assert shifted == [-1, -81, -33, 50, 41, 11, 1]
    # sign argument on [0, 1/2] (which contains [0, sqrt2-1] since sqrt2 < 3/2): 50x^3 <= 25x^2 < 33x^2 for x > 0,
    # 41x^4 + 11x^5 + x^6 <= (41/8 + 11/16 + 1/32) x < 6x < 81x
    half = R(1, 2)
    assert 50 * half == 25 < 33 and R(41, 8) + R(11, 16) + R(1, 32) == R(187, 32) < 6 < 81
    assert sqrt(2) < R(3, 2) and sp.expand(r2.subs(Qv, 1 + sqrt(2))) == 0
    # the bounds as exact polynomial inequalities: x^2(25 - 50x) >= 0 and x(187/32 - 41x^3 - 11x^4 - x^5) >= 0 on [0, 1/2]
    assert sp.expand(25 * xx**2 - 50 * xx**3 - xx**2 * 25 * (1 - 2 * xx)) == 0
    for k, bound in ((4, R(41, 8)), (5, R(11, 16)), (6, R(1, 32))):            # x^k <= x * (1/2)^(k-1) on [0, 1/2]
        assert sp.expand(half**(k - 1) * xx - xx**k - xx * (half**(k - 1) - xx**(k - 1))) == 0 and half**(k - 1) > 0
    assert shifted[4] * half**3 + shifted[5] * half**4 + shifted[6] * half**5 == R(187, 32)
    # numerical sanity of S_6(2+x) < 0 on a grid of [0, 1/2], exact rationals
    assert all(sp.Poly(sp.expand(S6.subs(Qv, 2 + xx)), xx).eval(R(i, 200)) < 0 for i in range(0, 101))
    assert t.subs(Qv, 2) == R(3, 4) and sp.simplify(t.subs(Qv, 1 + sqrt(2)) - 1) == 0      # Q=2 <-> t=3/4, Q -> 1+sqrt2 <-> t -> 1
    # control: S_6 changes sign (the argument needs the upper end of x): S_6(4) > 0, so the sign argument needs the bounded x
    assert S6.subs(Qv, 4) > 0
    print('F: PASS exact: at a=1 (every a by D), tr(H_5 X_t) = -(Q^2-2Q-1) S_6(Q)/(2Q(Q^2-1)), S_6(2+x) = -1-81x-33x^2+50x^3+41x^4+11x^5+x^6 < 0 '
          'on [0, 1/2] (hence on [0, sqrt2-1]) by 50x^3 <= 25x^2, 41x^4+11x^5+x^6 <= (187/32)x < 6x, and Q^2-2Q-1 < 0 for Q < 1+sqrt2: '
          'tr(H_5 X_t) < 0 for 3/4 <= tau < 1 (control: S_6(4) > 0)')


def section_g():
    """H0 (s=1): r P = Y (r-1)^2 (r-y1)^2 + X (r-x0)^2 (r^2+p r+q) + 16 nu X Y r with nu=7+4sqrt3, y1=3+2sqrt3, x0=2sqrt3-3,
    p=26+16sqrt3, q=nu^2.  P is concave in X on [0,r] (the X^2 coefficient is -16 nu), so P >= the chord through its two axis values."""
    r3 = sqrt(3)
    nu, y1, x0, pp = 7 + 4 * r3, 3 + 2 * r3, 2 * r3 - 3, 26 + 16 * r3
    q0 = nu**2
    H = sp.diag(-24 - 8 * r3, 1, 9, 12 + 4 * r3, 12 + 12 * r3)
    H[1, 2] = H[2, 1] = -nu
    zz = x + sp.I * y
    q = qvec(zz, 1)
    val = sp.expand((q.H * H * q)[0])
    r = x**2 + y**2
    target = (y**2 * (r - 1)**2 * (r - y1)**2 + x**2 * (r - x0)**2 * (r**2 + pp * r + q0) + 16 * nu * x**2 * y**2 * r)
    assert sp.expand(sp.expand(r * val) - sp.expand(target)) == 0
    assert sp.expand(x0**2 * q0 - y1**2) == 0 and all(c_.is_positive for c_ in (nu, y1, x0, pp, q0))
    t, K, p = surface_q()
    P = Qv**4 + (8 + 2 * r3) * Qv**3 - (2 + 4 * r3) * Qv**2 - (20 + 14 * r3) * Qv + 1 + 4 * r3
    f0 = (H * weights_x(p)).trace()
    assert vanishes(f0 + (Qv**2 - 3) * (Qv**2 - 2 * Qv - 1) * P / (2 * Qv * (Qv**2 - 1)))
    d1, d2, d3 = (sp.diff(P, Qv, k) for k in (1, 2, 3))
    assert sp.expand(d1.subs(Qv, r3) - (28 + 12 * r3)) == 0 and sp.expand(P.subs(Qv, 2) - (33 - 24 * r3)) == 0
    assert all(c_ > 0 for c_ in sp.Poly(d3, Qv).all_coeffs()) and d2.subs(Qv, r3) > 0 and 33 - 24 * r3 < 0
    assert P.subs(Qv, R(5, 2)) > 0                       # control: H0 stops detecting before tau = 21/20 (Q=5/2)
    assert sp.expand(P.subs(Qv, R(11, 5)) - (34956 - 15540 * r3) / 625) == 0 and P.subs(Qv, R(11, 5)) > 0
    # Lower endpoint tau = 1/sqrt(3) (Q = sqrt(3)): explicit 4-term single-cut product decomposition on A|BC
    x_low = weights_x({k: v.subs(Qv, r3) for k, v in p.items()})
    w_im, w_re = (3 * r3 - 5) / 48, (3 * r3 + 5) / 48
    assert w_im > 0 and w_re > 0
    q_im_p, q_im_m = qvec(I * sqrt(y1), 1), qvec(-I * sqrt(y1), 1)
    q_re_p, q_re_m = qvec(sqrt(x0), 1), qvec(-sqrt(x0), 1)
    low_decomp = w_im * (q_im_p * q_im_p.H + q_im_m * q_im_m.H) + w_re * (q_re_p * q_re_p.H + q_re_m * q_re_m.H)
    assert (x_low - low_decomp).applyfunc(sp.simplify) == sp.zeros(5)
    # Uniqueness of H0 from double roots at z = i, i*sqrt(y1), sqrt(x0) on symmetric forms
    hu, hv, hw, he, ho, hvw = symbols('hu hv hw he ho hvw')
    Hsym = sp.diag(hu, hv, hw, he, ho)
    Hsym[1, 2] = Hsym[2, 1] = hvw
    rr = symbols('rr', positive=True)
    p_im = sp.expand(4 * rr * hu + rr**2 * (rr - 3)**2 * hv + (rr + 1)**2 * hw + rr * (rr - 1)**2 * he
                     + (rr - 1)**2 * ho + 2 * rr * (rr - 3) * (rr + 1) * hvw)
    p_re = sp.expand(4 * rr * hu + rr**2 * (rr + 3)**2 * hv + (rr - 1)**2 * hw + rr * (rr + 1)**2 * he
                     + (rr + 1)**2 * ho - 2 * rr * (rr + 3) * (rr - 1) * hvw)
    eqs = [p_im.subs(rr, 1), sp.diff(p_im, rr).subs(rr, 1),
           p_im.subs(rr, y1), sp.diff(p_im, rr).subs(rr, y1),
           p_re.subs(rr, x0), sp.diff(p_re, rr).subs(rr, x0)]
    mat_sys = Matrix([[sp.expand(eq).coeff(v_) for v_ in (hu, hv, hw, he, ho, hvw)] for eq in eqs])
    ns = mat_sys.nullspace()
    assert len(ns) == 1 and (ns[0] / ns[0][1] - Matrix([-24 - 8 * r3, 1, 9, 12 + 4 * r3, 12 + 12 * r3, -nu])).applyfunc(sp.simplify) == sp.zeros(6, 1)
    print('G: PASS H0 = diag(-24-8s3, 1, 9, 12+4s3, 12+12s3), H0_vw = -(7+4s3): r q^+H0q = Y(r-1)^2(r-y1)^2 + X(r-x0)^2(r^2+pr+q) '
          '+ 16 nu X Y r >= 0 (y1=3+2s3, x0=2s3-3, p=26+16s3, q=(7+4s3)^2); tr(H0 X_tau) = -(Q^2-3)(Q^2-2Q-1)P(Q)/(2Q(Q^2-1)) '
          'with P convex on Q >= sqrt3, P\'(sqrt3) = 28+12s3 > 0, P(2) = 33-24s3 < 0: negative for 1/sqrt3 < tau <= 3/4, zero at 1/sqrt3; '
          'explicit 4-term single-cut product decomposition of X_{1/sqrt3} and uniqueness of H0 from the 6 endpoint double roots checked')


def window_polys():
    t, h = symbols('t h', positive=True)                    # h = sqrt(1+t^2)
    Kk = 1 - t**2
    E = Kk * (h + t)**2
    p = dict(A=1, W=t, V=(E**2 - 2 * E * Kk + Kk) / (t * Kk), G=(E - 1) / t, E=E, O=Kk / t, s=1)
    return t, h, p


def section_h():
    t, h, p = window_polys()
    limits = {key: sp.limit(sp.cancel(value).subs(h, sqrt(1+t**2)), t, 1, dir='-')
              for key, value in p.items()}
    assert limits == dict(A=1, W=1, V=1, G=-1, E=0, O=0, s=1)
    # lam^2=1/a, so the formula covers every finite positive filter parameter.
    lam = symbols('lam', positive=True)
    u, w, v, _, _ = V4.parts(lam**2)
    endpoint = u*u.H/lam**2 + (w-lam**2*v)*(w-lam**2*v).H
    rays = [sp.kronecker_product(*([Matrix([1, sign*I*lam])]*3)) for sign in (1, -1)]
    mixture = sum((q*q.H for q in rays), sp.zeros(8))/(2*lam**2)
    assert (endpoint-mixture).applyfunc(sp.expand) == sp.zeros(8)
    assert sp.factor(endpoint.trace()) == (1+lam**2)**3/lam**2
    # Each partial transpose conjugates one local factor of each product ray.
    # Their Gram determinant is strictly positive at every lam>0, so each rank is exactly two.
    determinant = sp.expand((1+lam**2)**6-(1-lam**2)**6)
    assert sp.factor(determinant) == 4*lam**2*(lam**4+3)*(3*lam**4+1)
    assert sp.expand(determinant-4*lam**2*(3+10*lam**4+3*lam**8)) == 0
    for cut in range(3):
        transposed_rays = []
        for sign in (1, -1):
            factors = [Matrix([1, sign*I*lam]) for _ in range(3)]
            factors[cut] = factors[cut].conjugate()
            transposed_rays.append(sp.kronecker_product(*factors))
        target = sum((q*q.H for q in transposed_rays), sp.zeros(8))/(2*lam**2)
        assert (V4.pt(endpoint, cut)-target).applyfunc(sp.expand) == sp.zeros(8)
        frame = Matrix.hstack(*transposed_rays)
        assert sp.expand((frame.H*frame).det()-determinant) == 0
    frame = Matrix.hstack(*rays)
    assert sp.expand((frame.H*frame).det()-determinant) == 0
    print('H: PASS tau -> 1: (A,W,V,G,E,O) -> (1,1,1,-1,0,0); for every a>0, the normalized limit is the equal mixture '
          'of eta_+^(tensor 3) and eta_-^(tensor 3), eta_+/-=(sqrt(a)|0> +/- i|1>)/sqrt(a+1); fully separable, profile (2,2,2,2)')


def main():
    print('Proposition 4d: analytic supported witness on the surface of Proposition 4c; exact verification')
    section_a()
    section_b()
    section_c()
    section_d()
    section_e()
    section_f()
    section_g()
    section_h()
    print('ALL PASS')


if __name__ == '__main__':
    main()
