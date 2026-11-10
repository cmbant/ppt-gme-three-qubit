#!/usr/bin/env python3
"""Exact check of PPT THREE QUBIT cor:ppt3-noise-reused (DS-6; in noise_certificates.tex).

The level-two DPS dual witness Y_* of ppt_gme_sym_certificate.json (key 'surface'; optimal for rho_* among symmetric witnesses with
spectrum in (-1/4, 31/10) and tr(Y sigma) <= -0.004076 for the state sigma of cor:ppt3-extreme5) is used on
the rational-weight surface state rho_* = rho(312,266,603)/tr, i.e. Pythagorean data (m,n) = (19,7), a = 201/104,
tau = 133/156, built with the surface helpers of verify_ppt_gme_rank5_flat_profile.py. Checked exactly over Q(sqrt(-3)):
  1. rho_* is PSD of rank 5 with every partial transpose PSD of rank 5 (profile (5,5,5,5));
  2. Y_* Hermitian and cyclic-invariant, and its certificate (extension identity, Z1, Z2, P > 0) by certify_ppt_gme_dps2_sym.py;
  3. tr(Y rho_*) < -173/40000, tr Y / 8 = 1, white noise p = 2153/500000 still detected;
  4. uniform local depolarisation q = 81/50000 still detected (exact endpoint value; rates q_j <= q follow by
     composing local depolarising channels, as in the proof);
  5. -I/4 < Y < 31/10 I; every state within trace distance 129/100000 of rho_* is detected; around
     rho_o = (1-32/10000) rho_* + (32/10000) I/8 the ball of radius 331/1000000 is detected by Y and stays fully PPT
     (rho_o^T >= (p/8) I and the PPT stability radius 2c/(1+sqrt2) of lem:ppt3-dps, kappa=(1+sqrt2)/2 < 169/140);
  6. the PPT radii of the table of cor:ppt3-noise and of cor:ppt3-robust against the constants of
     certify_ppt_gme_dps2_sym.py;
  7. negative controls: p = 432/100000, q = 82/50000 and radii just beyond each stated bound are NOT certified.
"""
import sys
from pathlib import Path
import json

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import certify_ppt_gme_dps2_sym as C  # noqa: E402
import verify_ppt_gme_rank5_flat_profile as F  # noqa: E402

R = sp.Rational


def main():
    print('cor:ppt3-noise-reused: reuse the surface witness Y on the surface state rho(312,266,603)')
    p = F.weights(19, 7, R(201, 104))
    assert (p['al'], p['W'], p['A']) == (312, 266, 603) and 3 * p['W']**2 > p['al']**2 > p['W']**2
    num = F.numerator(p).applyfunc(sp.expand)
    rho = (num / sp.expand(num.trace())).applyfunc(sp.expand)
    assert F.psd_rank(rho) == 5 and [F.psd_rank(F.pt(rho, k)) for k in range(3)] == [5, 5, 5]
    print('state: PASS exact PSD rank 5, every partial transpose PSD of rank 5 (tau=133/156, a=201/104)')

    data = json.loads((HERE / 'ppt_gme_sym_certificate.json').read_text())
    y = C.check_certificate(data['surface'], data['denominator'])      # exact certificate incl. cyclic symmetry and Hermiticity
    print('witness: PASS Y_* Hermitian and cyclic; level-two certificate (144-coefficient identity, Z1, Z2, P > 0)')

    g = -C.expectation(y, rho)
    h = sp.trace(y) / 8
    assert h == 1 and sp.expand(g - R(173, 40000)) > 0
    white = lambda q: sp.expand(-(1 - q) * g + q * h)
    assert white(R(2153, 500000)) < 0 and not white(R(432, 100000)) < 0
    print(f'white: PASS tr(Y rho_*) = {sp.N(-g, 12)} < -173/40000; p=2153/500000 detected; p=432/100000 not (control)')

    depolarised = lambda q: C.depolarised(rho, q)           # D_q on every qubit
    assert sp.expand(C.expectation(y, depolarised(R(81, 50000)))) < 0
    assert not sp.expand(C.expectation(y, depolarised(R(82, 50000)))) < 0
    print('local: PASS uniform depolarisation q=81/50000 detected (smaller q_j by the semigroup argument); '
          'q=82/50000 not (control)')

    lo, hi = R(-1, 4), R(31, 10)
    assert C.positive_definite(y - lo * sp.eye(8)) and C.positive_definite(hi * sp.eye(8) - y)
    width = hi - lo                                       # 67/20
    kappa = (1 + R(99, 70)) / 2                           # > (1+sqrt2)/2, since 99^2 - 2*70^2 = 1 (lem:ppt3-dps)
    assert 99**2 - 2 * 70**2 == 1

    # GME-only radius around rho_*: g - width * r > 0
    assert sp.expand(g - width * R(129, 100000)) > 0 and not sp.expand(g - width * R(130, 100000)) > 0
    # joint ball around rho_o = (1-p1) rho_* + p1 I/8: detected by Y, and rho_o^T >= (p1/8) I (rho_*^T >= 0)
    p1, r0 = R(32, 10000), R(331, 1000000)
    centre = ((1 - p1) * rho + p1 * sp.eye(8) / 8).applyfunc(sp.expand)
    g1 = C.expectation(y, centre)
    assert sp.expand(g1 + width * r0) < 0 and not sp.expand(g1 + width * R(34, 100000)) < 0
    assert kappa * r0 < p1 / 8 and not kappa * R(332, 1000000) < p1 / 8
    # rho_o^T - (p1/8) I = (1-p1) rho_*^T, PSD by the rank check above, so rho_o^T >= (p1/8) I at every cut
    assert all(all(C.zero(e) for e in F.pt(centre, k) - (p1 / 8) * sp.eye(8) - (1 - p1) * F.pt(rho, k)) for k in range(3))
    print('ball: PASS -I/4 < Y < 31I/10; GME-only radius 129/100000 around rho_*; radius 331/1000000 around rho_o '
          'detected and fully PPT (kappa < 169/140); 34/100000 (GME) and 332/1000000 (PPT) not (controls)')

    # design-point table of cor:ppt3-noise and the full-rank state: PPT radius < c / kappa, c from the DS-1..DS-4 constants
    for tag, radius in ((R(327, 100), R('0.000101')), (R(10, 3), R('0.00057')), (R(169, 50), R('0.0002')),
                        ('state', R('0.00002'))):
        c = R(C.CONSTANTS[tag][3])
        assert kappa * radius < c
    print('table: PASS PPT radii 1.01e-4, 5.7e-4, 2.0e-4 (family) and 2.0e-5 (full-rank state) satisfy kappa*r < c')
    print('ALL PASS')


if __name__ == '__main__':
    main()
