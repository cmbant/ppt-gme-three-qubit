#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT cor:ppt3-extreme5 (in fully_ppt_three_qubit_gme.tex): a fully-PPT, genuinely multipartite entangled three-qubit state
of rank five whose three partial transposes all have rank five (profile (5,5,5,5)).

State: sigma = N / tr N with N = A u u^+ + W w w^+ + V v v^+ + G (w v^+ + v w^+) + E e e^+ + O o o^+,
u = |000> - s(|011>+|101>+|110>), s = alpha/A, v = |111>, w, e, o as in the note, and weights from the Pythagorean
parametrisation alpha = m^2-n^2, W = 2mn, h = m^2+n^2, A = a*alpha,
   K = alpha^2-W^2, O = K/W, E = K (h+W)^2/(alpha^2 A), G = alpha (E - alpha^2/A)/W,
   V = (alpha^2 E^2 - 2 E (alpha^2/A) K + (alpha^2/A)^2 K)/(W K).
Section 1 proves symbolically, for all (m, n, a), that the even PT block has determinant 0 and that two explicit vectors
span a kernel of the odd PT block (so the PT of the whole two-parameter family has rank <= 5 at cut A|BC; the other two
cuts follow by cyclic invariance, checked exactly for the state).  Section 2 takes (m, n, a) = (8, 3, 2) (integer table in the note) and checks in
exact arithmetic over Q(sqrt(-3)): rank-5 PSD for rho and every partial transpose (characteristic polynomials), and
tr(Y sigma) < 0 for the level-two DPS dual witness Y_* of ppt_gme_sym_certificate.json (key 'surface'), whose own certificate
(Z1, Z2, P >= 0 and the 144-coefficient product identity) is re-checked here with certify_ppt_gme_dps2_sym.py.  Section 4 checks extremality exactly.  Positivity of the family away from this point is NOT claimed.
"""
import itertools
import json
import sys
from pathlib import Path

import sympy as sp
from sympy import Rational as R, Matrix, sqrt, I
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import certify_ppt_gme_dps2_sym as C  # noqa: E402

OMEGA = R(-1, 2) + sqrt(3) * I / 2
EVEN, ODD = [0, 3, 5, 6], [1, 2, 4, 7]


def ket(bits):
    v = sp.zeros(8, 1)
    v[int(bits, 2)] = 1
    return v


def parts(s):
    u = ket('000') - s * (ket('011') + ket('101') + ket('110'))
    w = ket('001') + ket('010') + ket('100')
    e = ket('011') + OMEGA * ket('101') + OMEGA**2 * ket('110')
    o = ket('100') + OMEGA * ket('010') + OMEGA**2 * ket('001')
    return u, w, ket('111'), e, o


def weights(m, n, a):
    m, n, a = sp.sympify(m), sp.sympify(n), sp.sympify(a)
    al, W, h = m**2 - n**2, 2 * m * n, m**2 + n**2
    A = a * al
    K = al**2 - W**2
    O = K / W
    E = K * (h + W)**2 / (al**2 * A)
    be = al**2 / A
    G = al * (E - be) / W
    V = (al**2 * E**2 - 2 * E * be * K + be**2 * K) / (W * K)
    return dict(al=al, W=W, A=A, s=al / A, O=O, E=E, G=G, V=V, be=be)


def numerator(p):
    u, w, v, e, o = parts(p['s'])
    d = lambda x: x.H
    return (p['A'] * u * d(u) + p['W'] * w * d(w) + p['V'] * v * d(v) + p['G'] * (w * d(v) + v * d(w))
            + p['E'] * e * d(e) + p['O'] * o * d(o))


def psd_rank(matrix):
    """Rank r if the Hermitian matrix is PSD, else None (det(xI-M)=x^(n-r) q(x), q with strictly alternating signs)."""
    assert all(C.zero(x) for x in matrix - matrix.conjugate().T)
    dm = DomainMatrix.from_Matrix(matrix, extension=True)
    cs = [sp.expand(dm.domain.to_sympy(c)) for c in dm.charpoly()]
    n = matrix.rows
    sgn = [sp.sign((-1)**k * c) for k, c in enumerate(cs)]
    r = max(k for k in range(n + 1) if sgn[k] != 0)
    ok = all(sgn[k] > 0 for k in range(r + 1)) and all(sgn[k] == 0 for k in range(r + 1, n + 1))
    return r if ok else None


def pt(matrix, k):
    return C.partial_transpose(matrix, (2, 2, 2), {k})


def section1_family():
    """Symbolic in (m, n, a): PT_A even block singular, odd block has two explicit kernel vectors."""
    m, n, a = sp.symbols('m n a', positive=True)
    z = sp.Symbol('z')
    p = weights(m, n, a)
    pa = pt(numerator(p), 0)
    be_, bo = pa.extract(EVEN, EVEN), pa.extract(ODD, ODD)
    to_z = lambda x: sp.expand(x).subs(sqrt(3) * I, 2 * z + 1)      # omega -> z, z^2+z+1 = 0
    num, _ = sp.fraction(sp.cancel(sp.together(be_.applyfunc(to_z).det())))
    assert sp.rem(sp.Poly(sp.expand(num), z), sp.Poly(z**2 + z + 1, z)).is_zero
    al, W, O, E, be = (p[k] for k in ('al', 'W', 'O', 'E', 'be'))
    gram = Matrix([[2 * W + O / 2, sqrt(3) * O / 2], [sqrt(3) * O / 2, 3 * O / 2]])
    x = gram.inv() * Matrix([-sqrt(2) * al, 0])
    y = gram.inv() * Matrix([sqrt(2) * (be - E / 2), sqrt(6) * E / 2])
    for (c1, c2), last in ((x, 2), (y, 3)):
        k = sp.zeros(4, 1)
        k[0], k[1], k[last] = (c1 + I * c2) / sqrt(2), (c1 - I * c2) / sqrt(2), -1
        assert all(sp.simplify(t) == 0 for t in bo * k)
    print('family: PASS symbolic in (m,n,a): det(even PT block)=0 and two kernel vectors of the odd PT block')


def section2_state(m, n, a):
    p = weights(m, n, a)
    table = dict(A=R(110), W=R(48), O=R(721, 48), E=R(7931, 250), G=R(121, 25), V=R(42108, 625), s=R(1, 2))
    assert all(p[k] == v for k, v in table.items())
    num = numerator(p).applyfunc(sp.expand)
    trace = sp.expand(num.trace())
    assert trace == R(5441073, 10000)
    # integer form quoted in the note: 30000 N = 825000 ut ut^+ + 1440000 w w^+ + 2021184 v v^+ + 145200 (w v^+ + v w^+)
    #                                           + 951720 e e^+ + 450625 o o^+,  ut = 2|000> - |011> - |101> - |110>,  trace 16323219
    ut, w, v, e, o = parts(R(1))
    ut = 2 * ket('000') - ket('011') - ket('101') - ket('110')
    d = lambda x: x.H
    table_n = (825000 * ut * d(ut) + 1440000 * w * d(w) + 2021184 * v * d(v) + 145200 * (w * d(v) + v * d(w))
               + 951720 * e * d(e) + 450625 * o * d(o)).applyfunc(sp.expand)
    assert all(C.zero(x) for x in table_n - 30000 * num) and sp.expand(table_n.trace()) == 16323219
    sigma = (num / trace).applyfunc(sp.expand)
    assert psd_rank(sigma) == 5
    cyc = C.CYCLE
    assert all(C.zero(sigma[cyc[i], cyc[j]] - sigma[i, j]) for i, j in itertools.product(range(8), repeat=2))
    parity = [bin(i).count('1') % 2 for i in range(8)]
    assert all(C.zero(sigma[i, j]) for i, j in itertools.product(range(8), repeat=2) if parity[i] != parity[j])
    prof = [psd_rank(pt(sigma, k)) for k in range(3)]
    assert prof == [5, 5, 5], prof
    ea, oa = psd_rank(pt(sigma, 0).extract(EVEN, EVEN)), psd_rank(pt(sigma, 0).extract(ODD, ODD))
    assert (ea, oa) == (3, 2)
    print(f'state: PASS exact PSD rank 5, cyclic and parity invariant; PT ranks {prof} (cut A: even block {ea}, odd block {oa})')
    return sigma, num


def section3_witness(sigma):
    data = json.loads((HERE / 'ppt_gme_sym_certificate.json').read_text())
    y = C.check_certificate(data['surface'], data['denominator'])        # level-two certificate of Y_*: exact, all checks inside
    g = -C.expectation(y, sigma)
    h = sp.trace(y) / 8
    assert h == 1 and g > R(407, 100000)
    white = g / (g + h)
    assert R(407, 100407) < white and R(81, 20000) < R(407, 100407)                        # the 4.05e-3 of cor:ppt3-extreme5
    print(f'GME: PASS level-two DPS dual certificate valid; tr(Y sigma) = {sp.N(-g, 8)} < 0 (rational); '
          f'white-noise tolerance p < {sp.N(white, 5)} (exactly g > 407/100000, tr Y/8 = 1, g/(g+1) > 407/100407 > 81/20000)')
    return g


def negative_controls(m, n, a, sigma):
    p = weights(m, n, a)
    q = dict(p)
    q['E'] = p['E'] + R(1, 1000)                                          # leaves the surface
    nq = numerator(q)
    bad = (nq / sp.expand(nq.trace())).applyfunc(sp.expand)
    prof = [psd_rank(pt(bad, k)) for k in range(3)]
    assert prof != [5, 5, 5]
    print(f'control: PASS E -> E+1/1000 leaves the surface; PT profile test returns {prof} (None = not PSD)')
    # Another point of the same surface (a = 1/2): fully PPT of profile (5,5,5,5), but the witness is silent there.
    p2 = weights(m, n, R(1, 2))
    n2 = numerator(p2)
    s2 = (n2 / sp.expand(n2.trace())).applyfunc(sp.expand)
    assert psd_rank(s2) == 5 and [psd_rank(pt(s2, k)) for k in range(3)] == [5, 5, 5]
    data = json.loads((HERE / 'ppt_gme_sym_certificate.json').read_text())
    y = C.load_witness(data['surface'], data['denominator'])[0]
    assert C.expectation(y, s2) > 0
    print('control: PASS surface point a=1/2 is fully PPT of profile (5,5,5,5) but tr(Y sigma) > 0 (no detection)')


def kernel_columns(matrix):
    """Exact right-kernel basis of a matrix over Q(sqrt(-3)), as columns."""
    dm = DomainMatrix.from_Matrix(matrix, extension=True)
    ns = dm.nullspace()
    return ns.to_Matrix().T


def section4_extremality(sigma, p):
    """sigma is an extreme point of the set of fully-PPT states: any Hermitian face direction D = B X B^+ must have
    D^(T_c) killing ker sigma^(T_c) for all three cuts; the real-linear system for X (25 unknowns) has rank 24."""
    u, w, v, e, o = parts(p['s'])
    b = u.row_join(w).row_join(v).row_join(e).row_join(o)
    ker = [kernel_columns(pt(sigma, k)) for k in range(3)]
    assert all(kc.cols == 3 for kc in ker)
    basis = []
    for i in range(5):
        for j in range(i, 5):
            x = sp.zeros(5, 5)
            x[i, j] = x[j, i] = 1
            basis.append(x)
            if i != j:
                x = sp.zeros(5, 5)
                x[i, j], x[j, i] = I, -I
                basis.append(x)
    assert len(basis) == 25
    cols = []
    for x in basis:
        d = (b * x * b.H).applyfunc(sp.expand)
        col = []
        for k in range(3):
            m = (pt(d, k) * ker[k]).applyfunc(sp.expand)
            col += [t.coeff(I, 0) for t in m] + [t.coeff(I, 1) for t in m]      # real and imaginary parts in Q(sqrt 3)
        cols.append(col)
    system = Matrix(cols).T
    rank = DomainMatrix.from_Matrix(system, extension=True).rank()      # rank over Q(sqrt 3) = rank over R
    assert rank == 24, rank
    # the single null direction is sigma itself
    xs = (b.H * b).inv() * b.H * sigma * b * (b.H * b).inv()
    d = (b * xs * b.H).applyfunc(sp.expand)
    assert all(C.zero(t) for k in range(3) for t in pt(d, k) * ker[k])
    print('extremality: PASS exact rank 24 of 25 for the face-direction system; sigma is an extreme fully-PPT state')


def main():
    print('Corollary 4b: exact fully-PPT GME three-qubit state, rank 5, all partial transposes of rank 5')
    section1_family()
    sigma, num = section2_state(8, 3, 2)
    section3_witness(sigma)
    section4_extremality(sigma, weights(8, 3, 2))
    negative_controls(8, 3, 2, sigma)
    print('ALL PASS')


if __name__ == '__main__':
    main()
