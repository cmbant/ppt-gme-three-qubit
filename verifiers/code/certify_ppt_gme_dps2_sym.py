#!/usr/bin/env python3
"""Exact certificates for PPT THREE QUBIT cor:ppt3-noise, cor:ppt3-robust and the witness Y_* of cor:ppt3-noise-reused /
cor:ppt3-extreme5 (in noise_certificates.tex), from the symmetry-reduced data of ppt_gme_sym_certificate.json.

Data (untrusted numerical output of explore_ppt_gme_sym_dps2.py), integers in units 1/denominator:
  Y   = sum_k y_k E_k, 10 parameters, the invariant Hermitian basis of the cyclic / parity / T = K.S_BC symmetries
        (E_k built from a1=|000>, a2=|011>+|101>+|110>, e, ebar, v=|111>, w=|001>+|010>+|100>, o, obar);
  Zi  = F_even R_i,even F_even^+ + F_odd R_i,odd F_odd^+ with two real symmetric 8x8 blocks R, F the T-real basis
        |a1 a2> x {|00>, |11>, |01>+|10>, sqrt(-3)(|01>-|10>)} of A1 A2 B1 B2 (parity blocks).
Everything is exact (SymPy over Q(sqrt(-3)) and Fractions); no solver, sampling or floating point.  Zi >= 0 iff both blocks are
>= 0 (F has full column rank), tested by exact LDL.  P = V^+ (Y_{A1B} x I - Z1^{T_B} - Z2^{T_A2}) V is rational in the T-real basis
of Sym^2(A) x B and is tested by exact LDL; the product-vector identity is checked on all 144 coefficients; Y, partial transposes of
states are tested by characteristic polynomials.  Checked: the certificate for every witness, every published tolerance of
cor:ppt3-noise (white p, local q, GME ball, PPT ball, uniform bounds over the PPT interval) and of cor:ppt3-robust, with controls.
"""
import itertools
import json
from fractions import Fraction
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
R = sp.Rational
S3 = sp.sqrt(3) * sp.I                                   # sqrt(-3)
OMEGA = (-1 + S3) / 2
BRACKET = (R(32567, 10000), R(33886, 10000))
CYCLE = [(n >> 1) + 4 * (n & 1) for n in range(8)]
TRW = [1, 0, 3, 3, 3, 1, 0, 3, 3, 3]                     # tr E_k
PRINTED_RHO = [955592, -930700, 906857, 519407, 208476, 390647, -520181, 693325, 556418, 104]   # 10^7 tilde rho, paper


def zero(x): return sp.expand(x) == 0
def vec(s):
    v = sp.zeros(8, 1); v[int(s, 2)] = 1; return v


def y_basis():
    a1, a2 = vec('000'), vec('011') + vec('101') + vec('110')
    e = vec('011') + OMEGA * vec('101') + OMEGA**2 * vec('110')
    v, w = vec('111'), vec('001') + vec('010') + vec('100')
    o = vec('100') + OMEGA * vec('010') + OMEGA**2 * vec('001')
    d = lambda x, y=None: (x * (x if y is None else y).H).applyfunc(sp.expand)
    mix = lambda x, y: d(x, y) + d(y, x)
    return [d(a1), mix(a1, a2), d(a2), d(e), d(e.conjugate()), d(v), mix(v, w), d(w), d(o), d(o.conjugate())]
YB = y_basis()


def t_real_basis(nfree, nvec):
    """Columns |free index> x {|00>,|11>,|01>+|10>,sqrt(-3)(|01>-|10>)}: T = K.S_BC fixes each; ordered by (free, kind)."""
    cols = []
    for f in range(nfree):
        for kind in range(4):
            c = sp.zeros(nvec, 1); b = 4 * f
            if kind == 0: c[b] = 1
            elif kind == 1: c[b + 3] = 1
            elif kind == 2: c[b + 1] = c[b + 2] = 1
            else: c[b + 1], c[b + 2] = S3, -S3
            cols.append((c, f, kind))
    return cols
FZ = t_real_basis(4, 16)           # free index = (a1, a2)
PARITY_Z = [(f // 2 + f % 2 + (k >= 2)) % 2 for _, f, k in FZ]
FQ = sp.Matrix.hstack(*[c for c, _, _ in t_real_basis(3, 12)])            # T-real basis of Sym^2(A) x B (columns of V^+ ...V)
GQ = sp.diag(*[1, 1, 2, 6] * 3)


def ldl_pd(m):
    """Exact LDL of a rational symmetric matrix: True iff positive definite."""
    a = [[Fraction(int(sp.numer(x)), int(sp.denom(x))) for x in row] for row in m.tolist()]
    n = len(a)
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    for k in range(n):
        if a[k][k] <= 0: return False
        for i in range(k + 1, n):
            f = a[i][k] / a[k][k]
            for j in range(k + 1, n): a[i][j] -= f * a[k][j]
    return True


def real_form(m, cols, gram):
    """G^{-1} F^+ M F G^{-1}: rational symmetric because M is T-invariant and F is T-real over Q(sqrt(-3))."""
    r = (gram.inv() * cols.H * m * cols * gram.inv()).applyfunc(sp.expand)
    assert all(sp.im(x) == 0 for x in r) and all(sp.expand(sp.re(x)).is_Rational for x in r)
    return r.applyfunc(sp.re)


def digits(index, dims):
    out = []
    for d in reversed(dims):
        index, dg = divmod(index, d); out.append(dg)
    return tuple(reversed(out))


def partial_transpose(m, dims, sites):
    lab = [digits(i, dims) for i in range(m.rows)]; idx = {l: i for i, l in enumerate(lab)}
    return sp.Matrix(m.rows, m.cols, lambda i, j: m[idx[tuple(lab[j][k] if k in sites else x for k, x in enumerate(lab[i]))],
                                                    idx[tuple(lab[i][k] if k in sites else x for k, x in enumerate(lab[j]))]])


def positive_definite(m):
    """Hermitian M > 0 iff all coefficients of the characteristic polynomial alternate in sign (exact, over Q(sqrt(-3)))."""
    if not all(zero(x) for x in m - m.conjugate().T): return False
    dm = DomainMatrix.from_Matrix(m, extension=True)
    return all(sp.expand((-1)**k * dm.domain.to_sympy(c)) > 0 for k, c in enumerate(dm.charpoly()))


def load_witness(rec, den):
    y = sum((R(c, den) * M for c, M in zip(rec['y'], YB)), sp.zeros(8))
    assert sum(c * t for c, t in zip(rec['y'], TRW)) == 8 * den        # tr Y = 8 exactly
    zs = []
    for key in ('Z1', 'Z2'):
        z = sp.zeros(16)
        for par, tri in enumerate(rec[key]):
            blk = sp.zeros(8); it = iter(tri)
            for i in range(8):
                for j in range(i, 8): blk[i, j] = blk[j, i] = R(next(it), den)
            assert ldl_pd(blk)                                          # exact LDL of the 8x8 block
            f = sp.Matrix.hstack(*[c for (c, _, _), p in zip(FZ, PARITY_Z) if p == par])
            z += f * blk * f.H
        zs.append(z.applyfunc(sp.expand))
    return (y.applyfunc(sp.expand), *zs)


def extension_residual(y, z1, z2):
    v = sp.zeros(16, 12)
    for a1, a2, b in itertools.product(range(2), range(2), range(4)): v[8 * a1 + 4 * a2 + b, 4 * (a1 + a2) + b] = 1
    lifted = sp.Matrix(16, 16, lambda i, j: y[4 * (i // 8) + i % 4, 4 * (j // 8) + j % 4] if (i // 4) % 2 == (j // 4) % 2 else 0)
    return v.T * (lifted - partial_transpose(z1, (2, 2, 4), {2}) - partial_transpose(z2, (2, 2, 4), {1})) * v, v.T * v


def product_identity(y, z1, z2, p):
    """All 144 coefficients of |a|^2 <ab|Y|ab> - c^+ P c - <aab*|Z1|aab*> - <aa*b|Z2|aa*b> (key (r,s,i,j): a1*^r a1^s b_i* b_j...)."""
    co = {}
    def add(key, value): co[key] = co.get(key, 0) + value
    for i, j, nb in itertools.product(range(8), range(8), range(2)): add((i // 4 + nb, j // 4 + nb, i % 4, j % 4), y[i, j])
    for i, j in itertools.product(range(12), repeat=2): add((i // 4, j // 4, i % 4, j % 4), -p[i, j])
    for i, j in itertools.product(range(16), repeat=2):
        a1, a2, b = digits(i, (2, 2, 4)); c1, c2, d = digits(j, (2, 2, 4))
        add((a1 + a2, c1 + c2, d, b), -z1[i, j]); add((a1 + c2, a2 + c1, b, d), -z2[i, j])
    return len(co) == 144 and all(zero(x) for x in co.values())


def check_certificate(rec, den):
    """Exact level-two certificate of a witness record; returns Y (8x8 over Q(sqrt(-3)))."""
    y, z1, z2 = load_witness(rec, den)
    assert all(zero(y[CYCLE[i], CYCLE[j]] - y[i, j]) for i, j in itertools.product(range(8), repeat=2)) and all(zero(x) for x in y - y.H)
    p, gram = extension_residual(y, z1, z2)
    pr, gr = real_form(p, FQ, GQ), real_form(gram, FQ, GQ)
    assert ldl_pd(pr) and product_identity(y, z1, z2, p) and not ldl_pd(pr - gr / 100)   # last: damaged margin rejected
    return y


def expectation(y, state):
    value = sp.expand(sp.trace(y * state)); assert zero(value - sp.conjugate(value)) and value.is_Rational
    return value


def replace_qubit(state, site):
    out = sp.zeros(8); mask = 1 << (2 - site)
    for i, j in itertools.product(range(8), repeat=2):
        if (i & mask) == (j & mask): out[i, j] = (state[i & ~mask, j & ~mask] + state[i | mask, j | mask]) / 2
    return out


def depolarised(state, q):
    out = sp.zeros(8)
    for mask in range(8):
        m = state
        for k in range(3):
            if mask & (1 << k): m = replace_qubit(m, k)
        out += (1 - q)**(3 - mask.bit_count()) * q**mask.bit_count() * m
    return out.applyfunc(sp.expand)


BASIS = sp.Matrix([[1, 0, 0, 0, 0], [0, 0, 1, 0, OMEGA**2], [0, 0, 1, 0, OMEGA], [-1, 0, 0, 1, 0], [0, 0, 1, 0, 1],
                   [-1, 0, 0, OMEGA, 0], [-1, 0, 0, OMEGA**2, 0], [0, 1, 0, 0, 0]])
def family_numerator(t): return (BASIS * sp.diag(t, 3 * t, 3, 3, 1) * BASIS.H).applyfunc(sp.expand)
def family(t): return family_numerator(t) / (7 * (t + 3))


# Published constants: white p, local q, GME ball, PPT lower bound c of the state, outward spectral bounds lo < Y < hi, and the
# control white-noise level just beyond the witness' own threshold.
CONSTANTS = {
    R(327, 100): ('0.00033', '0.00014', '0.00010', '0.0001226', '-0.0827', '3.1389', '0.00034'),
    R(10, 3): ('0.00087', '0.00037', '0.00027', '0.0006906', '-0.1380', '3.0809', '0.000887'),
    R(169, 50): ('0.00140', '0.00061', '0.00044', '0.0002497', '-0.1757', '3.0290', '0.00145'),
    'state': ('0.0233', '0.0100', '0.0073', '0.0000249', '-0.4147', '2.8568', '0.0234'),
}


def check_witness(tag, y, state, den):
    white, local, radius, c, lo, hi, control = map(R, CONSTANTS[tag])
    assert zero(sp.trace(state) - 1)
    h, value = sp.trace(y) / 8, expectation(y, state)
    assert h == 1 and value < 0 and (1 - white) * value + white * h < 0 and not (1 - control) * value + control * h < 0
    assert positive_definite(y - lo * sp.eye(8)) and positive_definite(hi * sp.eye(8) - y) and sp.expand(value + radius * (hi - lo)) < 0
    assert sp.expand(expectation(y, depolarised(state, local))) < 0
    assert all(positive_definite(partial_transpose(state, (2, 2, 2), {k}) - c * sp.eye(8)) for k in range(3))
    print(f'{tag}: PASS white p={white} (control {control} fails), local q<={local}, GME D-ball={radius}, PPT eigenvalues > {c}; '
          f'tr(Y rho) = {value} ~ {float(value):.5e}, {lo} < Y < {hi}')
    return value, h


def main():
    print('Exact symmetry-reduced level-two DPS certificates (Y: 10 parameters, Z1, Z2: 72 each)')
    assert BRACKET[0] < (15 + 9 * sp.sqrt(17)) / 16 < 2 * (sp.sqrt(37) - 1) / 3 < BRACKET[1]
    data = json.loads((HERE / 'ppt_gme_sym_certificate.json').read_text()); den = data['denominator']
    for rec in data['family']:
        t = R(*rec['t'])
        y = check_certificate(rec, den)
        print(f'{t}: PASS extension identity (all 144 coefficients), block LDL of Z1, Z2 and P > 0, cyclic symmetry; negative control')
        value, h = check_witness(t, y, family(t), den)
        b = expectation(y, family_numerator(0)); a = expectation(y, family_numerator(1)) - b
        assert zero((a * t + b) / (7 * (t + 3)) - value)               # 7(s+3) tr(Y rho_s) = a s + b
        if t in (R(327, 100), R(10, 3)):
            assert all(sp.expand(a * s + b) < 0 for s in BRACKET)
            print(f'{t}: PASS negative on the full rational bracket {BRACKET}')
        else:
            assert a < 0 and sp.expand(a * R('3.27475') + b) < 0
            print(f'{t}: PASS negative for every s >= 3.27475')
        if t == R(327, 100):
            lo, hi = R(CONSTANTS[t][4]), R(CONSTANTS[t][5])
            for endpoint in BRACKET:
                g = expectation(y, family(endpoint))
                assert sp.expand((1 - R('0.00024')) * g + R('0.00024')) < 0 and sp.expand(g + R('0.000075') * (hi - lo)) < 0
            print('family: PASS uniform white p=0.00024 and GME D-ball=0.000075 throughout the PPT interval')
    rob = data['robust']
    state = sum((R(c, den) * M for c, M in zip(rob['rho'], YB)), sp.zeros(8)).applyfunc(sp.expand)
    assert rob['rho'] == PRINTED_RHO and den == 10**7      # the coordinates displayed in appendix_noise.tex
    assert sum(c * t for c, t in zip(rob['rho'], TRW)) == den
    assert all(zero(state[CYCLE[i], CYCLE[j]] - state[i, j]) for i, j in itertools.product(range(8), repeat=2))
    parity = [bin(i).count('1') % 2 for i in range(8)]
    assert all(zero(state[i, j]) for i, j in itertools.product(range(8), repeat=2) if parity[i] != parity[j])
    assert positive_definite(state)
    print('state: PASS trace one, full rank, exact cyclic and global-parity symmetry (10 parameters)')
    y = check_certificate(rob, den)
    print('state: PASS extension identity (all 144 coefficients), block LDL of Z1, Z2 and P > 0, cyclic symmetry; negative control')
    check_witness('state', y, state, den)
    check_certificate(data['surface'], den)
    print('surface: PASS certificate of Y_* (used by verify_ppt_gme_reused_witness_noise.py and verify_ppt_gme_rank5_flat_profile.py)')
    print('ALL PASS')


if __name__ == '__main__':
    main()
