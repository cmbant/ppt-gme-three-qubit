#!/usr/bin/env python3
"""Independent exact re-derivation of the GME and PPT claims of PPT THREE QUBIT thm:ppt3-family (in fully_ppt_three_qubit_gme.tex).

Class: exact symbolic (sympy, arithmetic in Q(omega), omega^2+omega+1 = 0).  Does NOT reuse verify_chiral_ppt_gme_family.py.

What is recomputed from scratch, for each of the three cuts separately (no appeal to cyclic transport):
  E1  K B = 0 and rank K = 3, so the row space of K is the left annihilator of R = ran(B), dim R = 5.
  E2  The product vectors a (x) b of the cut that lie in R: M(z) = K((1,z)(x)I_4) has no common zero of its four 3x3 minors
      (Groebner basis = [1]) and M(infinity) has rank 3, so the product kernel is one-dimensional for every point of the
      projective line.  The kernel is written down by cofactors and expanded in the (orthogonal) basis B.
  E3  The witness form c(z)^dag H_5 c(z) on that curve equals 3 (r-1)^2 [(r-5)^2 + 12 x^2], r = |z|^2, x = Re z, and the
      point at infinity gives a positive value; hence it is >= 0 on all product vectors of R at every cut.
  E4  The H_s identity of prop:ppt3-witness(b) for every s > 0: with q_s(z) as printed in the paper and z = x + i y,
      sum_i h_i |q_s,i(z)|^2 = (r-s)^2 [(r-5s)^2 + 12 s x^2], H_s = diag(-12 s^3, 1, 11 s^2, -6 s, 14 s^2); at s = 1 this
      is the E3 form up to the factor 3 of the B-normalisation.
  E5  tr(H_5 X) = (29-9t)/(7(t+3)) for rho_t = B X B^dag, so it is negative for t > 29/9, and t_- > 29/9.
  E6  PPT interval: for every cut the coefficients e_k(t) of the characteristic polynomial of the partial transpose
      (a Hermitian matrix) are >= 0 exactly on [t_-, t_+] (exact polynomial positivity sets, determinant zero at both ends).
  E7  Negative control: H_5 with its w-weight lowered from 11 to 10 is negative at z = i, so the witness check can fail.
Logical chain (same as the note, only the computations are independent): a biseparable decomposition of rho_t has all pure
components in ran(rho_t) = R, each a product vector across some cut, so tr(H_5 X) >= 0; E5 gives tr(H_5 X) < 0.
"""
import itertools
import sympy as sp

w, z, zb, t, lam = sp.symbols('w z zb t lam')
MIN = w**2 + w + 1


def red(e):
    e = sp.expand(e)
    return sp.rem(sp.Poly(e, w), sp.Poly(MIN, w)).as_expr() if e.has(w) else e


def conj(e):
    return red(sp.expand(e).subs({w: w**2, z: zb, zb: z}, simultaneous=True))


def ket(s):
    v = [0] * 8
    v[int(s, 2)] = 1
    return v


def vadd(*vs):
    return [sum(c) for c in zip(*vs)]


def sc(a, v):
    return [a * x for x in v]


u = vadd(ket('000'), sc(-1, ket('011')), sc(-1, ket('101')), sc(-1, ket('110')))
v = ket('111')
ww = vadd(ket('001'), ket('010'), ket('100'))
e = vadd(ket('011'), sc(w, ket('101')), sc(w**2, ket('110')))
o = vadd(ket('100'), sc(w, ket('010')), sc(w**2, ket('001')))
B = sp.Matrix([u, v, ww, e, o]).T.applyfunc(red)
K = sp.Matrix([[0, 0, 0, 1, 0, w, w**2, 0], [0, w**2, w, 0, 1, 0, 0, 0], [3, 0, 0, 1, 0, 1, 1, 0]])
NORMS = [4, 1, 3, 3, 3]
H5 = sp.diag(-12, 1, 11, -6, 14)
results = []


def report(name, ok, msg):
    results.append(ok)
    print(f"{name} {'PASS' if ok else 'FAIL'} exact symbolic: {msg}")


def nonzero(expr):
    assert expr.free_symbols <= {w}  # Only an exact element of Q(omega), never a symbolic nonvanishing claim.
    return sp.rem(sp.Poly(sp.expand(expr), w), sp.Poly(MIN, w)).as_expr() != 0


def perm_idx(c):
    order = [c] + [i for i in range(3) if i != c]
    mp = []
    for n in range(8):
        bits = [(n >> (2 - j)) & 1 for j in range(3)]
        old = [0] * 3
        for j, qq in enumerate(order):
            old[qq] = bits[j]
        mp.append(old[0] * 4 + old[1] * 2 + old[2])
    return mp


# E1
gram = sp.Matrix(5, 5, lambda i, j: red(sum(conj(B[n, i]) * B[n, j] for n in range(8))))
rankK = any(nonzero(K[:, list(cs)].det()) for cs in itertools.combinations(range(8), 3))
report('E1', (K * B).applyfunc(red) == sp.zeros(3, 5) and rankK and gram == sp.diag(*NORMS),
       'K B = 0, rank K = 3, B^dag B = diag(4,1,3,3,3) (so dim R = 5 and R^perp = row space of K)')

# E2, E3 per cut
for c in range(3):
    mp = perm_idx(c)
    Bc = sp.Matrix(8, 5, lambda n, j: B[mp[n], j])
    Kc = sp.Matrix(3, 8, lambda i, n: K[i, mp[n]])
    K0, K1 = Kc[:, :4], Kc[:, 4:]
    M = K0 + z * K1
    mins = [red(M[:, [q for q in range(4) if q != j]].det()) for j in range(4)]
    G = sp.groebner([*mins, MIN], z, w, order='lex', domain=sp.QQ)
    inf_minors = [red(K1[:, [q for q in range(4) if q != j]].det()) for j in range(4)]
    inf_ok = any(nonzero(m) for m in inf_minors)
    kb = [((-1)**j) * mins[j] for j in range(4)]
    x = sp.Matrix(kb + [z * q for q in kb])
    cvec = sp.Matrix([red(sum(conj(Bc[n, i]) * x[n] for n in range(8))) / NORMS[i] for i in range(5)]).applyfunc(sp.expand).applyfunc(red)
    recon = (Bc * cvec - x).applyfunc(red) == sp.zeros(8, 1)
    # point at infinity: b in ker K1, x = (0, b)
    kinf = [((-1)**j) * inf_minors[j] for j in range(4)]
    xinf = sp.Matrix([0, 0, 0, 0] + kinf)
    cinf = sp.Matrix([red(sum(conj(Bc[n, i]) * xinf[n] for n in range(8))) / NORMS[i] for i in range(5)]).applyfunc(red)
    recon_inf = (Bc * cinf - xinf).applyfunc(red) == sp.zeros(8, 1)
    finf = red(sum(conj(cinf[i]) * H5[i, i] * cinf[i] for i in range(5)))
    report(f'E2[cut {c}]', list(G.exprs) == [1] and inf_ok and recon and recon_inf,
           'no common zero of the 3x3 minors (Groebner = [1]); rank 3 at infinity; kernel vector expanded in B and verified')
    form = red(sum(conj(cvec[i]) * H5[i, i] * cvec[i] for i in range(5)))
    r = z * zb
    xr = (z + zb) / 2
    target = 3 * (r - 1)**2 * ((r - 5)**2 + 12 * xr**2)
    finf_val = sp.simplify(finf.subs(w, sp.Rational(-1, 2) + sp.sqrt(3) * sp.I / 2))
    report(f'E3[cut {c}]', sp.expand(form - target) == 0 and finf_val > 0,
           'c^dag H_5 c = 3 (|z|^2-1)^2 [(|z|^2-5)^2 + 12 (Re z)^2] >= 0; value at infinity = %s > 0' % finf_val)
    if c == 0:
        cv0 = cvec

# E4: the H_s identity of prop:ppt3-witness(b), every s > 0 (polynomial identity in real x, y, s)
sx, sy, ss = sp.symbols('sx sy ss', real=True)
zc = sx + sp.I * sy
qs = [2 * zc, -zc**2 * (zc**2 + 3 * ss), zc**2 - ss, -zc * (zc**2 + ss), zc**2 + ss]
hs = [-12 * ss**3, 1, 11 * ss**2, -6 * ss, 14 * ss**2]
fs = sp.expand(sum(h * sp.expand(q * sp.conjugate(q)) for h, q in zip(hs, qs)))
rs = sx**2 + sy**2
ts = (rs - ss)**2 * ((rs - 5 * ss)**2 + 12 * ss * sx**2)
report('E4', sp.expand(fs - ts) == 0 and sp.expand(fs.subs(ss, 1) - ts.subs(ss, 1)) == 0,
       'q_s^dag H_s q_s = (r-s)^2 [(r-5s)^2 + 12 s x^2] identically in (x, y, s); H_1 = H_5')

# E5
Xd = [t, 3 * t, 3, 3, 1]
trH = sum(H5[i, i] * Xd[i] for i in range(5)) / (7 * (t + 3))
trace = sum(NORMS[i] * Xd[i] for i in range(5)) / (7 * (t + 3))
tm = (15 + 9 * sp.sqrt(17)) / 16
tp = 2 * (sp.sqrt(37) - 1) / 3
report('E5', sp.simplify(trH - (29 - 9 * t) / (7 * (t + 3))) == 0 and sp.simplify(trace) == 1 and tm > sp.Rational(29, 9),
       'tr(H_5 X) = (29-9t)/(7(t+3)), trace rho_t = 1, t_- = %s > 29/9 so tr(H_5 X) < 0 on the PPT interval' % sp.N(tm, 8))


# E6
def pt(M, c):
    out = sp.zeros(8, 8)
    for i in range(8):
        for j in range(8):
            bi = [(i >> (2 - q)) & 1 for q in range(3)]
            bj = [(j >> (2 - q)) & 1 for q in range(3)]
            bi2, bj2 = list(bi), list(bj)
            bi2[c], bj2[c] = bj[c], bi[c]
            out[bi2[0] * 4 + bi2[1] * 2 + bi2[2], bj2[0] * 4 + bj2[1] * 2 + bj2[2]] = M[i, j]
    return out


def proj(x):
    xm = sp.Matrix(x)
    return xm * xm.applyfunc(conj).T


num = (t * proj(u) + 3 * t * proj(v) + 3 * proj(ww) + 3 * proj(e) + proj(o)).applyfunc(red)
for c in range(3):
    Mc = pt(num, c)
    herm = (Mc - Mc.applyfunc(conj).T).applyfunc(red) == sp.zeros(8, 8)
    cp = sp.Poly(red(sp.expand((lam * sp.eye(8) - Mc).det(method='berkowitz'))), lam)
    co = [sp.expand(red(a)) for a in cp.all_coeffs()]
    ek = [sp.Poly(sp.expand((-1)**q * co[q]), t) for q in range(1, 9)]
    # Independently intersect all eight exact polynomial positivity sets.
    # This checks the entire t>0 domain, not merely two points outside it.
    allowed = sp.Interval.open(0, sp.oo)
    for p in ek:
        intervals = sp.solvers.inequalities.solve_poly_inequality(p, '>=')
        allowed = sp.Intersection(allowed, sp.Union(*intervals))
    det = ek[7]
    ends = sp.simplify(det.eval(tm)) == 0 and sp.simplify(det.eval(tp)) == 0
    report(f'E6[cut {c}]', herm and ends and allowed == sp.Interval(tm, tp),
           'PT Hermitian; exact intersection of all e_k(t) >= 0 for t>0 is [t_-,t_+]; exact polynomial inequalities and endpoint determinants')

# E7 negative control: lowering the w-weight of H_5 from 11 to 10 breaks nonnegativity on the curve
Hc = sp.diag(-12, 1, 10, -6, 14)
zv = sp.I
fc = red(sum(conj(cv0[i]) * Hc[i, j] * cv0[j] for i in range(5) for j in range(5)))
valc = sp.simplify(red(sp.expand(fc.subs({z: zv, zb: sp.conjugate(zv)}))))
report('E7', valc < 0, 'negative control: the same form with diag(-12,1,10,-6,14) is %s < 0 at z = i, so the check can fail' % valc)

print('ALL PASS' if all(results) else 'SOME FAIL')
raise SystemExit(0 if all(results) else 1)
