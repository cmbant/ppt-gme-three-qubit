#!/usr/bin/env python3
"""Numerical evidence: noise tolerance of the three-qubit fully-PPT GME family rho_t
(family.tex, eq:ppt3-rho).

RELAXATION (outer approximation of the biseparable cone, unnormalised):
    rho = S_A + S_B + S_C,   S_c >= 0 (8x8), S_c in DPS2_c, where for the cut c (qubit c | other two qubits, 2 x 4):
    S_c = tr_{A2} Omega,  Omega on A1 (x) A2 (x) B  (2x2x4 = 16 dim; the qubit side is extended ONCE, i.e. two
    copies of the qubit = level-2 of DPS), with
      * Omega = P_sym Omega P_sym  (supported on Sym^2(C^2) (x) C^4, 12-dim; Omega = V W V^dagger, W>=0 12x12),
      * Omega^{T_B} >= 0 (full transpose of the 4-dim side) and Omega^{T_{A2}} >= 0.
    (Omega^{T_{A1}} >= 0 follows from the symmetry.)  Level-1 relaxation: S_c >= 0, S_c^{T_B} >= 0.
    Cyclic reduction: rho_t and white / local-depolarising noise are invariant under the cyclic qubit permutation C,
    and the relaxed set is C-invariant and convex, so WLOG S_B, S_C are C-images of S_A. (Checked against the
    unreduced 3-cut problem at t=10/3.)  Sanity states that are not cyclic use the full 3-cut problem.
DISTANCES:  robustness  r(rho) = min{ r>=0 : rho + r I/8 in cone }.  Then (1-p) rho + p I/8 is in the cone iff
    p <= r/(1+r), so p* = r*/(1+r*) with no bisection.  GME is certified for the relaxed set when r > 0.
    The local-depolarising threshold q* is found by bisection on q with the test r(rho_q) > RTOL.
    Trace distance d_T = min (1/2)||rho - sigma||_1 over sigma in cone, tr sigma = 1.
SOLVER: CLARABEL, single threaded (status "optimal_inaccurate" = Clarabel AlmostSolved, final gap ~1e-8 abs; the
witness_lower_bound() weak-duality check below gives an independent bound on r).  No random numbers are used.
"""
import os, sys, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import itertools
import numpy as np, scipy, scipy.sparse as sp, cvxpy as cp

T0 = time.time()
TOL = 1e-8           # Clarabel tolerances (Clarabel stalls at ~1e-8 gap here: status optimal_inaccurate)
RTOL = 1e-6          # r above this counts as 'infeasible' (certified GME); solver error is ~1e-9 or smaller
SOLVER_KW = dict(solver=cp.CLARABEL, tol_gap_abs=TOL, tol_gap_rel=TOL, tol_feas=TOL, max_iter=200)

# ---------------- linear algebra helpers ----------------
def ket(bits):
    v = np.zeros(8, complex); v[int(bits, 2)] = 1; return v
om = np.exp(2j * np.pi / 3)
u = ket('000') - ket('011') - ket('101') - ket('110')
v = ket('111'); w = ket('001') + ket('010') + ket('100')
e = ket('011') + om * ket('101') + om**2 * ket('110')
o = ket('100') + om * ket('010') + om**2 * ket('001')
def rho_t(t):
    out = lambda x: np.outer(x, x.conj())
    return (t*out(u) + 3*t*out(v) + 3*out(w) + 3*out(e) + out(o)) / (7*(t+3))
tm = (15 + 9*np.sqrt(17))/16; tp = 2*(np.sqrt(37)-1)/3

def perm_op(dims, perm):
    """unitary U with U (x1 x ... xn) = (x_{perm[0]} x ...) i.e. reorders tensor factors."""
    n = len(dims); D = int(np.prod(dims)); U = np.zeros((D, D))
    for idx in itertools.product(*[range(d) for d in dims]):
        a = np.ravel_multi_index(idx, dims)
        new = [idx[p] for p in perm]
        b = np.ravel_multi_index(new, [dims[p] for p in perm])
        U[b, a] = 1
    return U
def pt_matrix(dims, which):
    """matrix Pm with vec_F(PT_which(M)) = Pm vec_F(M), partial transpose on the factors in `which`."""
    D = int(np.prod(dims)); rows = []; cols = []
    for ia in itertools.product(*[range(d) for d in dims]):
        for ib in itertools.product(*[range(d) for d in dims]):
            ja = tuple(ib[k] if k in which else ia[k] for k in range(len(dims)))
            jb = tuple(ia[k] if k in which else ib[k] for k in range(len(dims)))
            r = np.ravel_multi_index(ia, dims); c = np.ravel_multi_index(ib, dims)
            r2 = np.ravel_multi_index(ja, dims); c2 = np.ravel_multi_index(jb, dims)
            rows.append(r2 + D*c2); cols.append(r + D*c)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(D*D, D*D))
def linmap_matrix(f, n_in, n_out):
    """complex-linear map f: C^{n_in x n_in} -> C^{n_out x n_out}, as a matrix on vec_F."""
    M = np.zeros((n_out*n_out, n_in*n_in), complex)
    for k in range(n_in*n_in):
        E = np.zeros(n_in*n_in, complex); E[k] = 1
        M[:, k] = f(E.reshape(n_in, n_in, order='F')).reshape(-1, order='F')
    return M

# Level-2 DPS ingredients for a cut qubit | 4-dim (ordering A1 A2 B)
dimsE = [2, 2, 4]
SW = perm_op(dimsE, [1, 0, 2]); Psym = (np.eye(16) + SW) / 2
ev, evec = np.linalg.eigh(Psym); V = evec[:, ev > 0.5]          # 16 x 12 isometry onto Sym^2 (x) C^4
assert V.shape == (16, 12)
def embed(W): return V @ W @ V.conj().T
def trA2(Om):
    return np.einsum('ajkdjl->akdl', Om.reshape(2, 2, 4, 2, 2, 4)).reshape(8, 8)   # trace out A2
M_S = linmap_matrix(lambda W: trA2(embed(W)), 12, 8)             # W -> S (cut qubit first)
PT_B = pt_matrix(dimsE, {2}); PT_A2 = pt_matrix(dimsE, {1})
M_EMB = linmap_matrix(embed, 12, 16)
M_PTB = PT_B @ M_EMB; M_PTA2 = PT_A2 @ M_EMB
PT8 = pt_matrix([2, 4], {1})                                      # level-1: PT on the 4-dim side

# qubit permutations of the 8-dim space (ABC order)
CYC = perm_op([2, 2, 2], [2, 0, 1])      # cyclic shift
CYC_lin = lambda X: CYC @ X @ CYC.conj().T
M_CYC = linmap_matrix(CYC_lin, 8, 8)
def cut_to_abc(c):
    """matrix taking S (qubit c first, rest in order) to ABC order."""
    perm = [c] + [k for k in range(3) if k != c]
    U = perm_op([2, 2, 2], perm)          # U x = reordered (cut first)
    return linmap_matrix(lambda X: U.conj().T @ X @ U, 8, 8)

def vecF(X): return np.asarray(X).reshape(-1, order='F')
def unvec(x, n): return cp.reshape(x, (n, n), order='F')

def build(rho, level=2, cyclic=True, mode='robust', Wcert=False):
    """returns (problem, objective param symbols). rho: 8x8 Hermitian numpy."""
    cuts = [0] if cyclic else [0, 1, 2]
    cons = []; total = 0
    Ws = []
    for c in cuts:
        if level == 2:
            W = cp.Variable((12, 12), hermitian=True); Ws.append(W); cons += [W >> 0]
            vW = cp.vec(W, order='F')
            S = unvec(M_S @ vW, 8)
            cons += [unvec(M_PTB @ vW, 16) >> 0, unvec(M_PTA2 @ vW, 16) >> 0]
        else:
            S = cp.Variable((8, 8), hermitian=True); Ws.append(S); cons += [S >> 0]
            cons += [unvec(PT8 @ cp.vec(S, order='F'), 8) >> 0]
        Sabc = unvec(cut_to_abc(c) @ cp.vec(S, order='F'), 8)
        if cyclic:
            Sabc1 = unvec(M_CYC @ cp.vec(Sabc, order='F'), 8)
            Sabc2 = unvec(M_CYC @ cp.vec(Sabc1, order='F'), 8)
            total = total + Sabc + Sabc1 + Sabc2
        else:
            total = total + Sabc
    I8 = np.eye(8) / 8
    if mode == 'robust':
        r = cp.Variable(nonneg=True)
        con = (total == rho + r * I8); cons.append(con)
        prob = cp.Problem(cp.Minimize(r), cons)
        return prob, r, con
    if mode == 'feas':
        con = (total == rho); cons.append(con)
        return cp.Problem(cp.Minimize(0), cons), None, con
    if mode == 'dist':
        Pm = cp.Variable((8, 8), hermitian=True); Nm = cp.Variable((8, 8), hermitian=True)
        cons += [Pm >> 0, Nm >> 0, cp.real(cp.trace(total)) == 1, rho - total == Pm - Nm]
        return cp.Problem(cp.Minimize(cp.real(cp.trace(Pm) + cp.trace(Nm)) / 2), cons), None, None

def robustness(rho, level=2, cyclic=True):
    prob, r, con = build(rho, level, cyclic, 'robust')
    prob.solve(**SOLVER_KW)
    st = prob.status
    dual = None
    try:
        es = prob.solver_stats.extra_stats
        dual = es.get('obj_val_dual', None) if isinstance(es, dict) else getattr(es, 'obj_val_dual', None)
    except Exception:
        pass
    return (r.value if r.value is not None else np.nan), st, dual

def inner_min(Y, level=2, cyclic=True):
    """m(Y) = min <Y, sigma> over sigma in the relaxed cone with tr sigma = 1."""
    cuts = [0] if cyclic else [0, 1, 2]
    total = 0; cons = []
    for c in cuts:
        if level == 2:
            W = cp.Variable((12, 12), hermitian=True); cons += [W >> 0]; vW = cp.vec(W, order='F')
            S = unvec(M_S @ vW, 8)
            cons += [unvec(M_PTB @ vW, 16) >> 0, unvec(M_PTA2 @ vW, 16) >> 0]
        else:
            S = cp.Variable((8, 8), hermitian=True); cons += [S >> 0, unvec(PT8 @ cp.vec(S, order='F'), 8) >> 0]
        Sabc = unvec(cut_to_abc(c) @ cp.vec(S, order='F'), 8)
        if cyclic:
            S1 = unvec(M_CYC @ cp.vec(Sabc, order='F'), 8); S2 = unvec(M_CYC @ cp.vec(S1, order='F'), 8)
            total = total + Sabc + S1 + S2
        else:
            total = total + Sabc
    cons.append(cp.real(cp.trace(total)) == 1)
    prob = cp.Problem(cp.Minimize(cp.real(cp.trace(Y @ total))), cons)
    prob.solve(**SOLVER_KW)
    return prob.value

def witness_lower_bound(rho, level=2, cyclic=True):
    """Independent check of r*: take the equality-constraint multiplier Y of the robustness SDP, shift it to
    Y' = Y - m(Y) I so that <Y',sigma> >= 0 on the relaxed cone (m from a second SDP), and use weak duality
    <Y', rho + r I/8> >= 0  =>  r >= -<Y',rho> / (tr Y'/8).  Best of +-Y.  Same solver accuracy as everything else."""
    prob, r, con = build(rho, level, cyclic, 'robust'); prob.solve(**SOLVER_KW)
    Y0 = np.array(con.dual_value); Y0 = (Y0 + Y0.conj().T) / 2
    best = -np.inf
    for sgn in (+1, -1):
        Y = sgn * Y0; m = inner_min(Y, level, cyclic); Yp = Y - m * np.eye(8)
        den = np.trace(Yp).real / 8
        if den > 1e-12: best = max(best, -np.trace(Yp @ rho).real / den)
    return best

def feasible(rho, level=2, cyclic=True):
    """True/False by the robustness SDP (feasible iff r<=RTOL)."""
    r, st, _ = robustness(rho, level, cyclic)
    return r <= RTOL, r, st

def tracedist(rho, level=2, cyclic=True):
    prob, _, _ = build(rho, level, cyclic, 'dist')
    prob.solve(**SOLVER_KW)
    return prob.value, prob.status

def local_depol(rho, q):
    """(D_q)^{(x)3}: D_q(X) = (1-q) X + q tr(X) I/2, applied to each qubit via tr(X) I/2 = (1/4) sum_P P X P."""
    paulis = [np.eye(2), np.array([[0,1],[1,0]]), np.array([[0,-1j],[1j,0]]), np.diag([1.,-1.])]
    out = rho
    for ax in range(3):
        ops = []
        for P in paulis:
            m = [np.eye(2)]*3; m[ax] = P
            ops.append(np.kron(np.kron(m[0], m[1]), m[2]))
        twirl = sum(P @ out @ P.conj().T for P in ops) / 4     # = tr_ax(out) (x) I/2
        out = (1-q)*out + q*twirl
    return out

def fmt(x, n=3): return f"{x:.{n}e}"
fails = []
def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)
    if not ok: fails.append(name)

# ------------------------------------------------------------------
print("explore_ppt_gme_noise_tolerance.py -- numerical evidence, NOT a proof")
print(__doc__)
print(f"python {sys.version.split()[0]}  numpy {np.__version__}  scipy {scipy.__version__}  cvxpy {cp.__version__}  "
      "")
try:
    import clarabel; print("clarabel", clarabel.__version__)
except Exception: pass
print(f"solver CLARABEL tol_gap_abs=tol_gap_rel=tol_feas={TOL}, max_iter=200; GME threshold RTOL={RTOL}; no seed (deterministic)")
print(f"t_- = {tm:.6f}  t_+ = {tp:.6f}\n")

# cyclic invariance of rho_t
rt = rho_t(10/3)
_Z = local_depol(np.eye(8)/8, 0.3); _r = rho_t(3.3); _ld = local_depol(_r, 1.0)
check("local_depol: trace preserving, q=1 gives I/8, q=0 identity", abs(np.trace(local_depol(_r, 0.2))-1) < 1e-12 and np.allclose(_ld, np.eye(8)/8) and np.allclose(local_depol(_r, 0), _r))
check("rho_t cyclic invariant & trace 1", np.allclose(CYC @ rt @ CYC.conj().T, rt) and abs(np.trace(rt)-1) < 1e-12)
# PPT check numerically
def ptq(rho, k):
    R = rho.reshape(2,2,2,2,2,2); ax = [k, k+3]; R = R.swapaxes(*ax); return R.reshape(8,8)
check("rho_{10/3} is PPT at all 3 cuts (min eig PT)", min(np.linalg.eigvalsh(ptq(rt,k)).min() for k in range(3)) > -1e-12,
      "min eig " + fmt(min(np.linalg.eigvalsh(ptq(rt,k)).min() for k in range(3))))

# sanity states
ghz = np.zeros(8, complex); ghz[0] = ghz[7] = 1/np.sqrt(2); rghz = np.outer(ghz, ghz.conj())
prod = np.outer(ket('000'), ket('000'))
I8 = np.eye(8)/8
print("\n-- (i) level-1 (fully PPT) relaxation on rho_{10/3}")
t1 = time.time()
r1, st1, _ = robustness(rt, level=1, cyclic=True)
check("(i) rho_{10/3} feasible in level-1 relaxation (r = 0)", r1 <= 1e-8, f"r={fmt(r1)} status={st1}")
print("-- (ii),(iii) sanity states, full (non-cyclic-reduced) 3-cut level-2 problem")
ok, r, st = feasible(prod, 2, False); check("(ii) product |000> feasible", ok, f"r={fmt(r)} {st}")
# GHZ + white noise is biseparable iff p >= 4/7 = 0.571 (fidelity witness F>1/2 and Duer-Cirac); p = 0.6 is above
ok, r, st = feasible(0.4*rghz + 0.6*I8, 2, False); check("(ii) 0.4 GHZ + 0.6 I/8 feasible (known biseparable)", ok, f"r={fmt(r)} {st}")
ok, r, st = feasible(rghz, 2, False);  check("(iii) GHZ infeasible (GME detected)", not ok, f"r={fmt(r)} {st}")
r2, st2, _ = robustness(rghz, 2, False)
pg = r2/(1+r2)
print(f"      GHZ white-noise threshold in the relaxation: p* = {pg:.5f}  (exact GME threshold of GHZ+white noise is 4/7 = {4/7:.5f} by the fidelity witness F>1/2"
      f"; an outer relaxation can only give p* <= 4/7)")
check("GHZ p* = 4/7 (relaxation is tight for GHZ + white noise)", abs(pg - 4/7) < 1e-5, f"p*={pg:.6f}")
ok_c, rc, _ = feasible(rghz, 2, True)
check("cyclic-reduced problem agrees with full 3-cut problem on GHZ", abs(rc - r2) < 1e-6, f"{fmt(rc)} vs {fmt(r2)}")
tf = time.time(); rfull, stf, _ = robustness(rt, 2, False); tfull = time.time()-tf
tc = time.time(); rcyc, stc, dualc = robustness(rt, 2, True); tcyc = time.time()-tc
check("cyclic-reduced = full 3-cut robustness at t=10/3", abs(rfull - rcyc) < 1e-7,
      f"full {fmt(rfull,6)} ({tfull:.1f}s), cyclic {fmt(rcyc,6)} ({tcyc:.1f}s)")
print(f"      white-noise robustness r(rho_10/3) = {rcyc:.6e}  p* = {rcyc/(1+rcyc):.6e}")

print("\n-- (iv) trace distance of rho_{10/3} to level-2 relaxed set")
d, stt = tracedist(rt, 2, True)
print(f"      min (1/2)||rho - sigma||_1 = {d:.6e}   (||.||_1 itself = {2*d:.6e}); status {stt}; EF127 quotes 3.19e-4")
check("(iv) trace distance within 25% of 3.19e-4 (either half-norm or 1-norm convention)",
      abs(d-3.19e-4)/3.19e-4 < 0.25 or abs(2*d-3.19e-4)/3.19e-4 < 0.25, f"d_T={d:.4e}")
d1, _ = tracedist(rt, 1, True)
print(f"      level-1 trace distance (should be ~0): {d1:.3e}")

# ------------------------------------------------------------------
print("\n== main table ==  (r = white-noise robustness, margin = r itself vs RTOL; q* by bisection on r(rho_q) > RTOL)")
ts = [tm + 0.005, 3.30, 10/3, 3.37, tp - 0.005]
rows = []
for t in ts:
    rho = rho_t(t)
    r, st, dual = robustness(rho, 2, True)
    pstar = r/(1+r)
    rlb = witness_lower_bound(rho, 2, True)
    # direct check: noisy states just below/above p*
    def noisy(p): return (1-p)*rho + p*I8
    rin, _, _ = robustness(noisy(0.98*pstar), 2, True); rout, _, _ = robustness(noisy(1.02*pstar), 2, True)
    # local depolarising: bisection on q
    lo, hi = 0.0, 0.001
    rhi = robustness(local_depol(rho, hi), 2, True)[0]
    while rhi > RTOL and hi < 0.5:
        lo = hi; hi *= 2; rhi = robustness(local_depol(rho, hi), 2, True)[0]
    for _ in range(16):
        mid = (lo+hi)/2
        rm = robustness(local_depol(rho, mid), 2, True)[0]
        if rm > RTOL: lo = mid
        else: hi = mid
    qstar = lo
    rq = robustness(local_depol(rho, qstar), 2, True)[0]
    rows.append((t, r, pstar, qstar, rin, rout, rq, st, rlb))
    print(f"t={t:.5f}  r={r:.4e} p*={pstar:.4e} q*={qstar:.4e} | r(0.98p*)={rin:.2e} r(1.02p*)={rout:.2e} r(q*)={rq:.2e} "
          f"status={st} witness-bound r>={rlb:.4e}  [{time.time()-T0:.0f}s]", flush=True)

print("\n  t        p* (white)   q* (local depol)   r at t (robustness)")
for t, r, p, q, *_ in rows: print(f"  {t:.5f}  {p:.4e}   {q:.4e}        {r:.4e}")
best_p = max(rows, key=lambda x: x[2]); best_q = max(rows, key=lambda x: x[3])
print(f"most robust (white noise): t={best_p[0]:.5f}; (local depolarising): t={best_q[0]:.5f}")
check("all five t values certified GME in the DPS-2 relaxation (r > RTOL)", all(r[1] > RTOL for r in rows))
check("direct feasibility brackets p*: r(0.98p*) = 0.02 r(1+..) > RTOL (infeasible) and r(1.02p*) <= RTOL (feasible)",
      all(x[4] > RTOL and x[5] <= RTOL for x in rows), "r(0.98p*), r(1.02p*) printed per row")
check("witness lower bound agrees with SDP value to 1% (independent weak-duality check)",
      all(abs(x[8]-x[1]) < 0.01*x[1] for x in rows), "see 'witness-bound' per row")
print(f"\nSANITY FAILURES: {fails if fails else 'none'}")
print(f"total runtime {time.time()-T0:.0f} s")
