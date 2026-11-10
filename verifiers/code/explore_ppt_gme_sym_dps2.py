#!/usr/bin/env python3
"""Generator (numerical, NOT a proof) of the symmetry-reduced level-two DPS witness data checked by certify_ppt_gme_dps2_sym.py.

Statement in noise_certificates.tex, cor:ppt3-noise, cor:ppt3-noise-reused, cor:ppt3-robust, cor:ppt3-extreme5.
Writes ppt_gme_sym_certificate.json.  Replaces explore_ppt_gme_dps2_witness.py and explore_ppt_gme_robust_state.py.

All states and witnesses are invariant under the cyclic qubit permutation C, the global parity Z(x)Z(x)Z and the antiunitary
T = K.S_BC (swap qubits B and C, then complex-conjugate in the computational basis).  The level-two certificate is imposed with
the same symmetries (Z1, Z2 even under the parity of A1 A2 B and under T = K.S_BC; the optimal value is unchanged), so that

  Y  = sum_k y_k E_k           10 real parameters, E_k built from the sector vectors a1=|000>, a2=|011>+|101>+|110>, e, ebar,
                               v=|111>, w=|001>+|010>+|100>, o, obar  (e, o the cube-root-of-unity charge vectors);
  Zi = F_even R_i,even F_even^+ + F_odd R_i,odd F_odd^+      two real symmetric 8x8 blocks (72 parameters),

with F the T-real basis |a1 a2> x {|00>,|11>,|01>+|10>, sqrt(-3)(|01>-|10>)} of the four qubits A1 A2 B C.  Zi >= 0 iff both
blocks are >= 0.  The objective is tr(Y rho) with tr Y = 8; the margins GAMMA (blocks R >= GAMMA I) and DELTA (P >= DELTA G0) make
the rounded data (integers / DEN) strictly positive.  Entry points: witnesses for the family rho_t at the three design values
(with an interval constraint at the bracket endpoint for t = 327/100), the surface witness Y_* (optimal for rho_* under
tr(Y sigma) <= -0.004076 and -I/4 <= Y <= 3.1 I), the seesaw for the robust state (start and mixing as in the remark of
noise_certificates.tex) and its margin witness.  Deterministic, no random numbers.
"""
import os, sys, json, time
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import itertools
from fractions import Fraction
import numpy as np, scipy.sparse as sp, cvxpy as cp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "ppt_gme_sym_certificate.json")
DEN = 10**7                 # every stored integer is in units of 1/DEN
GAMMA, DELTA = 1e-6, 1e-6   # margins of the blocks R and of P
EPS = 2.5e-4                # white-noise weight mixed into the seesaw optimum (strict PPT)
TS = [Fraction(327, 100), Fraction(10, 3), Fraction(169, 50)]
KW = dict(solver=cp.CLARABEL, tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-9, max_iter=300)
OM = np.exp(2j * np.pi / 3)


def ket(s):
    v = np.zeros(2**len(s), complex); v[int(s, 2)] = 1; return v
def outer(x, y=None): return np.outer(x, (x if y is None else y).conj())

a1, a2 = ket('000'), ket('011') + ket('101') + ket('110')
e = ket('011') + OM * ket('101') + OM**2 * ket('110'); eb = e.conj()
v, w = ket('111'), ket('001') + ket('010') + ket('100')
o = ket('100') + OM * ket('010') + OM**2 * ket('001'); ob = o.conj()
# E_k: invariant Hermitian basis (order of the JSON field "y"); traces 1,3,3,3,1,3,3,3 on diagonal parameters 0,2,3,4,5,7,8,9
YB = [outer(a1), outer(a1, a2) + outer(a2, a1), outer(a2), outer(e), outer(eb),
      outer(v), outer(v, w) + outer(w, v), outer(w), outer(o), outer(ob)]
TRW = [1, 0, 3, 3, 3, 1, 0, 3, 3, 3]
CYC = np.zeros((8, 8))
for n in range(8): CYC[(n >> 1) + 4 * (n & 1), n] = 1


def combo(y): return sum(c * M for c, M in zip(y, YB))
def sector_params(M):
    """Real coordinates of a cyclic/parity/T-invariant Hermitian matrix in the basis YB (least squares)."""
    A = np.array([np.concatenate([B.real.ravel(), B.imag.ravel()]) for B in YB]).T
    return np.linalg.lstsq(A, np.concatenate([M.real.ravel(), M.imag.ravel()]), rcond=None)[0]


def zbasis():
    cols, par = [], []
    for a, b in itertools.product(range(2), repeat=2):
        for kind in range(4):
            f = np.zeros(16, complex); base = 8 * a + 4 * b
            if kind == 0: f[base] = 1
            elif kind == 1: f[base + 3] = 1
            elif kind == 2: f[base + 1] = f[base + 2] = 1
            else: f[base + 1], f[base + 2] = 1j * np.sqrt(3), -1j * np.sqrt(3)
            cols.append(f); par.append((a + b + (kind >= 2)) % 2)
    return np.array(cols).T, np.array(par)
FZ, PZ = zbasis()
FEV, FOD = FZ[:, PZ == 0], FZ[:, PZ == 1]


def pt_matrix(dims, which):
    D = int(np.prod(dims)); rows, cols = [], []
    for ia in itertools.product(*[range(d) for d in dims]):
        for ib in itertools.product(*[range(d) for d in dims]):
            ja = tuple(ib[k] if k in which else ia[k] for k in range(len(dims)))
            jb = tuple(ia[k] if k in which else ib[k] for k in range(len(dims)))
            rows.append(np.ravel_multi_index(ja, dims) + D * np.ravel_multi_index(jb, dims))
            cols.append(np.ravel_multi_index(ia, dims) + D * np.ravel_multi_index(ib, dims))
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(D * D, D * D))
PT_B, PT_A2 = pt_matrix([2, 2, 4], {2}), pt_matrix([2, 2, 4], {1})
PT3 = pt_matrix([2, 2, 2], {0})
Pi = np.zeros((16, 16))
for a1_, a2_, b_ in itertools.product(range(2), range(2), range(4)):
    Pi[a1_ * 8 + a2_ * 4 + b_, a1_ * 8 + b_ * 2 + a2_] = 1
V0 = np.zeros((16, 12)); _c = 0
for pair in [((0, 0),), ((0, 1), (1, 0)), ((1, 1),)]:
    for b_ in range(4):
        for (x, y_) in pair: V0[x * 8 + y_ * 4 + b_, _c] = 1
        _c += 1
G0 = V0.T @ V0


def rho_t(t):
    u = ket('000') - ket('011') - ket('101') - ket('110')
    return (t * outer(u) + 3 * t * outer(v) + 3 * outer(w) + 3 * outer(e) + outer(o)) / (7 * (t + 3))
def surface(m, n, a):                    # rho(m,n,a) of thm:ppt3-orbit, trace one (weights as in verify_ppt_gme_rank5_flat_profile)
    al, W, h = m * m - n * n, 2 * m * n, m * m + n * n
    A = a * al; K = al**2 - W**2; O = K / W; E = K * (h + W)**2 / (al**2 * A); be = al**2 / A
    G = al * (E - be) / W; V = (al**2 * E**2 - 2 * E * be * K + be**2 * K) / (W * K)
    u = ket('000') - al / A * (ket('011') + ket('101') + ket('110'))
    N = A * outer(u) + W * outer(w) + V * outer(v) + G * (outer(w, v) + outer(v, w)) + E * outer(e) + O * outer(o)
    return N / np.trace(N).real


def pt_apply(M, K, D): return (K @ M.reshape(-1, order='F')).reshape(D, D, order='F')


def solve_witness(rho, delta=DELTA, gamma=GAMMA, lo=None, hi=None, extra=()):
    """min tr(Y rho): Y in span(YB), tr Y = 8, blocks R >= gamma I, P >= delta G0, optional lo I <= Y <= hi I."""
    y = cp.Variable(10)
    Y = sum(y[k] * YB[k] for k in range(10))
    R1 = [cp.Variable((8, 8), symmetric=True) for _ in range(2)]
    R2 = [cp.Variable((8, 8), symmetric=True) for _ in range(2)]
    Z1 = FEV @ R1[0] @ FEV.conj().T + FOD @ R1[1] @ FOD.conj().T
    Z2 = FEV @ R2[0] @ FEV.conj().T + FOD @ R2[1] @ FOD.conj().T
    Yf = Pi @ cp.kron(Y, np.eye(2)) @ Pi.T
    pt = lambda M, Z: cp.reshape(M @ cp.vec(Z, order='F'), (16, 16), order='F')
    P = V0.T @ (Yf - pt(PT_B, Z1) - pt(PT_A2, Z2)) @ V0
    cons = [R - gamma * np.eye(8) >> 0 for R in R1 + R2] + [P - delta * G0 >> 0, y @ np.array(TRW) == 8]
    if lo is not None: cons.append(Y - lo * np.eye(8) >> 0)
    if hi is not None: cons.append(hi * np.eye(8) - Y >> 0)
    cons += [f(y, Y) for f in extra]
    prob = cp.Problem(cp.Minimize(cp.real(cp.trace(Y @ rho))), cons)
    prob.solve(**KW)
    return prob, y.value, [R.value for R in R1], [R.value for R in R2]


def tri(R): return [int(R[i, j]) for i in range(8) for j in range(i, 8)]
def rounded(prob_out, den=DEN):
    """Integer data of one witness; y[0] is adjusted so that tr Y = 8 exactly."""
    _, y, R1, R2 = prob_out
    yi = [int(round(c * den)) for c in y]
    yi[0] += 8 * den - sum(c * t for c, t in zip(yi, TRW))
    sym = lambda R: np.rint((R + R.T) / 2 * den)
    return {"y": yi, "Z1": [tri(sym(R)) for R in R1], "Z2": [tri(sym(R)) for R in R2]}


def numeric_g(rec, rho, den=DEN):
    return -np.trace(combo(np.array(rec["y"]) / den) @ rho).real


def seesaw(iters=300, tol=1e-12):
    """Alternate the witness step and the fully-PPT state step (cyclic, parity, T symmetric), as in Proposition prop:ppt3-noise-extremal."""
    u = ket('000') - 0.974 * (ket('011') + ket('101') + ket('110'))
    x = ket('001') + ket('010') + ket('100') - 0.750 * ket('111')
    wt = [0.368, 0.247, 0.156, 0.062, 0.167]
    rho = sum(c * outer(q) / np.vdot(q, q).real for c, q in zip(wt, [u, x, e, eb, o]))
    r = sector_params(rho); hist = []
    for it in range(iters):
        _, yv, _, _ = solve_witness(combo(r), delta=0.0, gamma=0.0)
        Y = combo(yv)
        rv = cp.Variable(10); Rm = sum(rv[k] * YB[k] for k in range(10))
        cons = [Rm >> 0, rv @ np.array(TRW) == 1, cp.reshape(PT3 @ cp.vec(Rm, order='F'), (8, 8), order='F') >> 0]
        pr = cp.Problem(cp.Minimize(cp.real(cp.trace(Y @ Rm))), cons); pr.solve(**KW)
        r = rv.value; hist.append(pr.value)
        if it % 10 == 0: print(f"  seesaw {it}: tr(Y rho) = {pr.value:.10e}, tolerance {-pr.value / (1 - pr.value):.7e}", flush=True)
        if it >= 25 and abs(hist[-26] - hist[-1]) < tol: break
    print(f"seesaw end after {it + 1} iterations: tolerance {-hist[-1] / (1 - hist[-1]):.7e}", flush=True)
    return r, hist[-1]


def main():
    print("explore_ppt_gme_sym_dps2.py -- numerical generator; the exact verifier is certify_ppt_gme_dps2_sym.py")
    t0 = time.time()
    data = {"denominator": DEN, "family": [], "surface": None, "robust": None}
    for t in TS:
        rho = rho_t(float(t)); extra = ()
        if t == Fraction(327, 100):       # uniform over the PPT interval: also detect rho_s at the lower bracket end
            lowend = rho_t(3.2567)
            extra = (lambda y, Y: cp.real(cp.trace(Y @ lowend)) <= -2.5e-4,)
        out = solve_witness(rho, extra=extra)
        rec = rounded(out); rec["t"] = [t.numerator, t.denominator]; data["family"].append(rec)
        g = numeric_g(rec, rho)
        print(f"t={float(t):.5f}: tr(Y rho) = {-g:.6e}, white-noise p* = {g / (g + 1):.6e}  [{time.time() - t0:.0f}s]", flush=True)
    rs, sg = surface(19, 7, 201 / 104), surface(8, 3, 2)
    sig = lambda y, Y: cp.real(cp.trace(Y @ sg)) <= -0.004076
    out = solve_witness(rs, lo=-0.2499, hi=3.0999, extra=(sig,))
    rec = rounded(out); data["surface"] = rec
    Yv = combo(np.array(rec["y"]) / DEN); ev = np.linalg.eigvalsh(Yv)
    print(f"Y_*: g(rho_*) = {numeric_g(rec, rs):.6e}, g(sigma) = {numeric_g(rec, sg):.6e}, spec Y in [{ev[0]:.4f}, {ev[-1]:.4f}]  [{time.time() - t0:.0f}s]", flush=True)
    r, val = seesaw()
    ro = combo(r)
    ev, U = np.linalg.eigh(ro)
    print("seesaw optimum: eigenvalues of rho", np.round(ev, 6))
    for j in np.argsort(ev)[::-1][:2]:
        vec = U[:, j] / U[:, j][np.argmax(abs(U[:, j]))]
        print(f"   eigenvalue {ev[j]:.6f}: eigenvector (largest entry 1)", np.round(vec, 5))
    print("   weights of e, ebar, o, obar:", [round(float(r[k] * 3), 6) for k in (3, 4, 8, 9)])
    rt = (1 - EPS) * r + EPS / 8 * np.array([1, 0, 1 / 3, 1 / 3, 1 / 3, 1, 0, 1 / 3, 1 / 3, 1 / 3])
    ri = [int(round(c * DEN)) for c in rt]
    ri[0] += DEN - sum(c * t for c, t in zip(ri, TRW))
    rho_q = combo(np.array(ri) / DEN)
    pt8 = lambda M: pt_apply(M, PT3, 8)
    print(f"rounded state: smallest eigenvalue {np.linalg.eigvalsh(rho_q)[0]:.4e}, smallest PT eigenvalue {np.linalg.eigvalsh(pt8(rho_q))[0]:.4e}")
    out = solve_witness(rho_q, lo=-0.4149, hi=2.8575)
    rec = rounded(out); rec["rho"] = ri; data["robust"] = rec
    g = numeric_g(rec, rho_q); Yv = combo(np.array(rec["y"]) / DEN); ev = np.linalg.eigvalsh(Yv)
    print(f"robust witness: tr(Y rho) = {-g:.6e}, white-noise tolerance {g / (g + 1):.6e}, spec Y in [{ev[0]:.4f}, {ev[-1]:.4f}]")
    with open(OUT, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    print("wrote", os.path.basename(OUT), f"[{time.time() - t0:.0f}s]")


if __name__ == '__main__':
    main()
