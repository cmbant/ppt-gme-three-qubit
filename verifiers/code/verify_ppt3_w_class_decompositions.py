#!/usr/bin/env python3
"""Exact regressions for w_class_decompositions.tex, WC-1--WC-4.

The theorem is analytic. Check the actual eight-component pure vectors and
their orbit-projector sums, the symbolic surface budgets, and exact rational
members of both branches. Zero hyperdeterminant certifies membership of the
mixed W class, not biseparability or the general PPT-to-W conjecture.
"""
import sympy as s

I, R = s.I, s.Rational
omega = (-1 + I * s.sqrt(3)) / 2
u = s.Matrix([1, 0, 0, -1, 0, -1, -1, 0])
v = s.eye(8)[:, 7]
w = s.Matrix([0, 1, 1, 0, 1, 0, 0, 0])
e = s.Matrix([0, 0, 0, 1, 0, omega, omega**2, 0])
o = s.Matrix([0, omega**2, omega, 0, 1, 0, 0, 0])
parity = s.diag(*[(-1)**k.bit_count() for k in range(8)])
cycle = s.zeros(8)
for k in range(8):
    cycle[((k & 1) << 2) | (k >> 1), k] = 1


def zero(x):
    if isinstance(x, s.MatrixBase):
        assert all(s.simplify(y) == 0 for y in x)
    else:
        assert s.simplify(x) == 0


def outer(x):
    return x * x.conjugate().T


def det3(x):
    a, b, c, d, f, g, h, j = x
    return s.expand(
        a*a*j*j + b*b*h*h + c*c*g*g + d*d*f*f
        - 2*(a*j*(b*h + c*g + d*f) + b*h*c*g + b*h*d*f + c*g*d*f)
        + 4*(a*d*g*h + j*b*c*f))


def cyclic(x):
    return sum((cycle**k * x * (cycle.T)**k for k in range(3)), s.zeros(8)) / 3


def pinch(x):
    return (x + parity * x * parity) / 2


def state(A, V, W, E, O, G=0):
    return (A*outer(u) + V*outer(v) + W*outer(w) + E*outer(e)
            + O*outer(o) + G*(w*v.T + v*w.T))


def diagonal_decomposition(A, V, W, E, O):
    r = min(A, E)
    q = A-r
    assert min(r, q, E-r, V-4*q, W, O) >= 0
    return (r*cyclic(outer(u+e)) + q*pinch(outer(u+2*v))
            + (E-r)*outer(e) + (V-4*q)*outer(v) + W*outer(w) + O*outer(o))


def main():
    z = s.symbols('z')
    b = s.symbols('b', real=True)
    x_sym, y_sym, c_sym = s.symbols('x_sym y_sym c_sym', real=True)
    g_ket = s.Matrix([1, 0, 0, 0, 0, 0, 0, 0])
    h_ket = s.Matrix([0, 0, 0, 1, 0, 1, 1, 0])
    zero(det3(u+z*e) - 4*(z**3-1))
    zero(det3(u+z*v) - (z*z-4))
    zero(det3(v+b*(w-o)))
    zero(det3((g_ket - (x_sym + y_sym)*h_ket) + x_sym*e + y_sym*e.conjugate()))
    zero(det3(x_sym*e + I*(v - b*w + c_sym*o)) - 4*(c_sym**3 - b**3 + 3*b*c_sym*x_sym**2))
    product = s.kronecker_product(*[s.Matrix([1, I])]*3)
    zero(u+I*(w-v)-product)
    assert det3(u) == -4
    zero(det3(v+b*w) - 4*b**3)  # Omitting o destroys the W identity.
    for x in (g_ket, h_ket, v, w, e, e.conjugate(), o, o.conjugate(), u+e, u+2*v, product):
        zero(det3(x))
    print('WC-1 PASS: Cayley hyperdeterminants, parity-tetrahedron & cross-sector atoms, product atom, GHZ and missing-o controls')

    zero(cyclic(outer(u+e)) - outer(u) - outer(e))
    zero(pinch(outer(u+2*v)) - outer(u) - 4*outer(v))
    zero(pinch(outer(product)) - outer(u) - outer(w-v))
    zero(cyclic(outer(v+b*(w-o))) - outer(v)
         - b*(w*v.T+v*w.T) - b*b*(outer(w)+outer(o)))
    for t in (R(327, 100), R(10, 3), R(169, 50), R(12)):
        zero(state(t, 3*t, 3, 3, 1) - diagonal_decomposition(t, 3*t, 3, 3, 1))
    a = (15+9*s.sqrt(17))/16
    v6 = 8*(9-a)/5
    zero(state(a, v6, 3, 3, 1) - diagonal_decomposition(a, v6, 3, 3, 1))
    # An exact member outside the sufficient diagonal budget must be rejected.
    try:
        diagonal_decomposition(R(1), R(0), R(0), R(0), R(0))
    except AssertionError:
        pass
    else:
        raise AssertionError('GHZ-only coefficient data accepted')
    print('WC-2 PASS: actual cyclic/parity projector sums and diagonal-family decompositions')

    t, E = s.symbols('tau E', positive=True)
    k = 1-t*t
    O = k/t
    G = (E-1)/t
    V = (E*E-2*k*E+k)/(t*k)
    zero(O*V-G*G-(2*E-1))
    h = (1-E)/t
    zero(t-h-(E-k)/t)
    zero(V-h-E*(E-k)/(t*k))
    zero(E-1+h-(1-E)*(1/t-1))
    Q = s.symbols('Q', positive=True)
    tq = (Q-1/Q)/2
    kq = 1-tq*tq
    Eq = kq*Q*Q
    zero((t*V-G*G).subs({t: tq, E: Eq}) - 2*Q*Q*kq*(Q*Q-3))
    # The analytic sign proof uses Q >= sqrt(3), Q < 1+sqrt(2),
    # E >= kappa >= 1-tau, and the indicated E >= 1 / E < 1 split.
    zero(Eq-kq-kq*(Q*Q-1))
    zero(k-(1-t)-t*(1-t))
    print('WC-3 PASS: symbolic surface Gram determinant and every residual-budget identity')

    branches = set()
    for q in (R(7, 4), R(2), R(7, 3), R(12, 5)):
        tq = (q-1/q)/2
        kq = 1-tq*tq
        Eq, Oq = kq*q*q, kq/tq
        Gq = (Eq-1)/tq
        Vq = (Eq*Eq-2*kq*Eq+kq)/(tq*kq)
        target = state(1, Vq, tq, Eq, Oq, Gq)
        assert q*q >= 3 and kq > 0 and Vq*tq-Gq*Gq >= 0
        if Eq >= 1:
            branches.add('E>=1')
            bq = Gq/Vq
            coefficients = (1, Vq, tq-Gq*Gq/Vq, Eq-1, Oq-Gq*Gq/Vq)
            rebuilt = (cyclic(outer(u+e)) + Vq*cyclic(outer(v+bq*(w-o)))
                       + coefficients[2]*outer(w) + coefficients[3]*outer(e)
                       + coefficients[4]*outer(o))
        else:
            branches.add('E<1')
            hq = -Gq
            coefficients = (hq, 1-hq, tq-hq, Vq-hq, Eq-1+hq, Oq)
            rebuilt = (hq*pinch(outer(product)) + (1-hq)*cyclic(outer(u+e))
                       + (tq-hq)*outer(w) + (Vq-hq)*outer(v)
                       + (Eq-1+hq)*outer(e) + Oq*outer(o))
        assert min(coefficients) >= 0
        zero(rebuilt-target)
        # The filter covariance is checked on all eight computational amplitudes.
        for lam in (R(2, 3), R(3, 2)):
            f = s.diag(*[lam**i.bit_count() for i in range(8)])
            zero(f*rebuilt*f-f*target*f)
            for x in (u, u+e, u+2*v, product, v+w-o,
                      s.Matrix([1, 2+I, 3, 4-I, 5, 6+I, 7, 8])):
                zero(det3(f*x)-lam**6*det3(x))
    assert branches == {'E>=1', 'E<1'}
    # Phase-coexistence line at t = -1/200:
    m_bs = s.Matrix([1232, 30493, 1693, 3700, 925, 32]) / 54375
    m_0 = s.Matrix([293, 577, 273, 188, 76, -105]) / 3360
    A_l, V_l, W_l, E_l, O_l, G_l = (1 - R(-1, 200)) * m_0 + R(-1, 200) * m_bs
    mu_l = -G_l
    assert mu_l == R(10930399, 348000000) and (A_l - mu_l) - E_l == R(7, 31250)
    zero(state(A_l, V_l, W_l, E_l, O_l, G_l)
         - (mu_l * pinch(outer(product)) + diagonal_decomposition(A_l - mu_l, V_l - mu_l, W_l - mu_l, E_l, O_l)))
    # Full-rank robust state tilde_rho (10^7 * tilde_rho):
    c0, c1, c2, ce, ceb, cv, cg, cw, co, cob = (
        R(955592), R(-930700), R(906857), R(519407), R(208476),
        R(390647), R(-520181), R(693325), R(556418), R(104))
    assert c0*c2 - c1**2 == 382804344 and cv*cw - cg**2 == 257058514
    s_r, b_r = -c1 / c0, -cg / cv
    x1, y1 = R(5, 9) * s_r, R(4, 9) * s_r
    x2_sq = (b_r**3 - 1) / (3 * b_r)
    x2 = s.sqrt(x2_sq)
    us_r = g_ket - s_r * h_ket
    atom_even = us_r + x1 * e + y1 * e.conjugate()
    atom_cross = x2 * e + I * (v - b_r * w + o)
    zero(det3(atom_even))
    zero(det3(atom_cross))
    r_h = c2 - c0 * s_r**2
    r_e = ce - c0 * x1**2 - cv * x2_sq
    r_eb = ceb - c0 * y1**2
    r_w = cw - cv * b_r**2
    r_o = co - cv
    assert min(c0, cv, r_h, r_e, r_eb, r_w, r_o, cob) > 0
    tilde_num = (c0 * outer(g_ket) + c1 * (g_ket * h_ket.T + h_ket * g_ket.T) + c2 * outer(h_ket)
                 + ce * outer(e) + ceb * outer(e.conjugate())
                 + cv * outer(v) + cg * (v * w.T + w * v.T) + cw * outer(w)
                 + co * outer(o) + cob * outer(o.conjugate()))
    tilde_rebuilt = (c0 * cyclic(outer(atom_even)) + cv * cyclic(pinch(outer(atom_cross)))
                     + r_h * outer(h_ket) + r_e * outer(e) + r_eb * outer(e.conjugate())
                     + r_w * outer(w) + r_o * outer(o) + cob * outer(o.conjugate()))
    zero(tilde_num - tilde_rebuilt)
    print('WC-4 PASS: exact surface members of both branches, phase-coexistence line, full-rank state tilde_rho, and local-filter covariance')
    print('ALL PASS: constructed-state mixed-W classification; general PPT-to-W and H/S/R/I remain open')


if __name__ == '__main__':
    main()
