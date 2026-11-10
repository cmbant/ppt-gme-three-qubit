#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT prop:ppt3-noise-sharp (in noise_certificates.tex): explicit
biseparable decompositions showing that white noise p_up destroys GME, for rho_* = rho(312,266,603)/tr
(cor:ppt3-noise-reused), rho_(10/3) (cor:ppt3-noise) and the rational full-rank state of ppt_gme_sym_certificate.json
(cor:ppt3-robust); p_up is read from the data file.

Data (untrusted numerical output of explore_ppt_gme_sym_upper.py): a few orbit representatives (z_r, B_r) with a_r = (1, z_r) a
Gaussian-rational qubit vector and B_r a 4x4 Gaussian-rational matrix.  The states are invariant under the cyclic permutation C,
the parity Z(x)Z(x)Z and T = K.S_BC, so each representative is expanded into its four images
g(|a><a| (x) B), g in {1, Z, T, ZT} with Z: (a, B) -> (Z a, ZZ B ZZ), T: (a, B) -> (conj a, S conj(B) S) (S swaps the qubits of BC).
Checked in exact arithmetic over Q(i, sqrt(-3)):
  1. every B_r is positive definite, so tau = sum_{r,g} g(|a_r><a_r| (x) B_r) is separable across A|BC;
  2. X = (1-p_up) rho + p_up I/8 - Cyc(tau), with Cyc the average over the three cyclic qubit permutations, is Hermitian
     and cyclic invariant, and ||X - c I||_2^2 <= c^2 with c = tr X / 8 > 0.  By Gurvits-Barnum (an operator within
     Frobenius distance 1 of the identity is separable) X/c, hence X, is separable across A|BC;
  3. therefore (1-p_up) rho + p_up I/8 = Cyc(tau + X) is a mixture of states separable across A|BC, B|CA and C|AB, i.e.
     biseparable, and so is every (1-p) rho + p I/8 with p >= p_up (a further mixture with the separable I/8);
  4. negative control: at the certified GME noise level (2153/500000, 87/100000, 233/10000) the same test fails, as it must.
Together with the certified lower bounds the GME white-noise thresholds are pinned to [lower, p_up] for each state.
"""
import itertools
import json
import sys
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import certify_ppt_gme_dps2_sym as C  # noqa: E402
import verify_ppt_gme_rank5_flat_profile as F  # noqa: E402

R = sp.Rational
CONTROL = {'rho_star': R(2153, 500000), 'rho_10/3': R(87, 100000), 'robust_state': R(233, 10000)}
ZZ = sp.diag(1, -1, -1, 1)                                  # Z (x) Z on the qubits of BC
SWAP = sp.Matrix(4, 4, lambda i, j: 1 if i == [0, 2, 1, 3][j] else 0)


def gauss(pair, den): return R(pair[0], den) + sp.I * R(pair[1], den)


def cyc(m):
    perm = sp.zeros(8)
    for i in range(8): perm[C.CYCLE[i], i] = 1
    return ((m + perm * m * perm.T + perm.T * m * perm) / 3).applyfunc(sp.expand)


def states():
    num = F.numerator(F.weights(19, 7, R(201, 104))).applyfunc(sp.expand)
    data = json.loads((HERE / 'ppt_gme_sym_certificate.json').read_text())
    rob = sum((R(c, data['denominator']) * M for c, M in zip(data['robust']['rho'], C.YB)), sp.zeros(8))
    return {'rho_star': (num / sp.expand(num.trace())).applyfunc(sp.expand), 'rho_10/3': C.family(R(10, 3)),
            'robust_state': rob.applyfunc(sp.expand)}


def images(a, b):
    """The four images g(|a><a| (x) B) of one representative."""
    za = sp.Matrix([a[0], -a[1]])
    cj = lambda m: m.conjugate()
    pairs = [(a, b), (za, ZZ * b * ZZ), (cj(a), SWAP * cj(b) * SWAP), (cj(za), SWAP * cj(ZZ * b * ZZ) * SWAP)]
    return [sp.kronecker_product(x * x.H, m) for x, m in pairs]


def gb_separable(x):
    """x Hermitian with ||x - cI||_2^2 <= c^2, c = tr x / 8 > 0 (Gurvits-Barnum ball)."""
    c = sp.expand(x.trace()) / 8
    d = x - c * sp.eye(8)
    norm2 = sp.expand(sum(sp.expand(e * sp.conjugate(e)) for e in d))
    return bool(c > 0) and bool(sp.expand(c**2 - norm2) >= 0), c, norm2


def main():
    data = json.loads((HERE / 'ppt_gme_noise_upper_certificate.json').read_text())
    rhos = states()
    print('prop:ppt3-noise-sharp: explicit biseparable decompositions of white-noise mixtures (orbit representatives)')
    for case in data['cases']:
        name, p_up = case['name'], R(*case['p_up'])
        rho = rhos[name]
        assert sp.expand(rho.trace()) == 1 and all(C.zero(e) for e in rho - cyc(rho))
        tau = sp.zeros(8)
        for rep in case['reps']:
            a = sp.Matrix([1, gauss(rep['z'], case['den_z'])])
            b = sp.Matrix([[gauss(z, case['den_b']) for z in row] for row in rep['B']])
            assert C.positive_definite(b)
            for image in images(a, b): tau += image
        cyc_tau = cyc(tau.applyfunc(sp.expand))

        def residual(p):
            return ((1 - p) * rho + p * sp.eye(8) / 8 - cyc_tau).applyfunc(sp.expand)
        x = residual(p_up)
        assert all(C.zero(e) for e in x - x.H) and all(C.zero(e) for e in x - cyc(x))
        ok, c, norm2 = gb_separable(x)
        assert ok and p_up > CONTROL[name]
        print(f'{name}: PASS {len(case["reps"])} orbit representatives (4 images each), positive definite B_r; residual X has tr X/8 = {sp.N(c, 4)}, '
              f'||X - cI||_2 = {sp.N(sp.sqrt(norm2), 4)} <= c (Gurvits-Barnum): (1-p) rho + p I/8 is biseparable '
              f'for every p >= {p_up}')
        bad, _, _ = gb_separable(residual(CONTROL[name]))
        assert not bad
        print(f'{name}: control PASS the same test fails at the certified GME level p = {CONTROL[name]}')
    print('ALL PASS')


if __name__ == '__main__':
    main()
