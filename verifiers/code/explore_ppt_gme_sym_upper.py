#!/usr/bin/env python3
"""Generator (numerical, NOT a proof) of the orbit-representative biseparable decompositions checked by certify_ppt_gme_noise_upper.py.

Statement in noise_certificates.tex, prop:ppt3-noise-sharp.  Writes ppt_gme_noise_upper_certificate.json.
Replaces explore_ppt_gme_noise_upper.py.

For a three-qubit state rho invariant under the cyclic permutation C, the parity Z(x)Z(x)Z and T = K.S_BC, the noisy state
(1-p) rho + p I/8 is biseparable if it equals the cyclic average of one A|BC-separable operator
    tau = sum_r sum_{g in {1, Z, T, ZT}} g( |a_r><a_r| (x) B_r ),   B_r >= 0 on BC,  a_r = (1, z_r) (unnormalised),
with Z: (a, B) -> (Z a, ZZ B ZZ) and T: (a, B) -> (conj a, S conj(B) S) (S swaps the two qubits of BC).  Averaging any feasible
tau over the group keeps feasibility, so restricting to group-closed tau loses nothing, and the data are the orbit
representatives (z_r, B_r) only.  Column generation over the representatives: for fixed directions minimise p (SDP, CVXOPT); the
symmetrised dual witness W of the equality constraint gives the most violated direction argmin_a lambda_min(<a|W|a>) (Nelder-Mead
from a fixed grid in the fundamental domain 0 <= phi <= pi/2); add up to three new orbits, drop inactive ones, repeat until no
direction is violated by more than 5e-7.  Then z_r are rounded to Gaussian rationals with denominator DEN_Z, the target is
p_up = ceil(10^6 (p_opt + 1.5e-6)) / 10^6 (just above the optimum, leaving room for the solver residual), the SDP maximises the
reserved separable noise delta at fixed p_up, and B_r (spectrum clipped at 3/DEN_B) are rounded to Gaussian rationals with denominator
DEN_B; the first (DEN_Z, DEN_B) for which the Gurvits-Barnum residual test of the exact verifier passes with a factor-two margin is
kept.  Cases: rho_* = rho(312,266,603)/tr, rho_(10/3) and the rational full-rank state of ppt_gme_sym_certificate.json.  For comparison
it prints the optimal level-two white-noise bound of each state (explore_ppt_gme_sym_dps2.py), a lower bound on the threshold.
Deterministic.
"""
import json, math, os, sys, time
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"
import numpy as np, cvxpy as cp
from scipy.optimize import minimize
import explore_ppt_gme_sym_dps2 as G

HERE = G.HERE
OUT = os.path.join(HERE, "ppt_gme_noise_upper_certificate.json")
KW = dict(solver=cp.CVXOPT)   # CLARABEL returned slightly infeasible B here
ZZ = np.kron(np.diag([1., -1]), np.diag([1., -1])); SW = np.eye(4)[[0, 2, 1, 3]]
PC = G.CYC
outer = G.outer


def cyc_avg(X): return (X + PC @ X @ PC.T + PC.T @ X @ PC) / 3
def cj(M): return M.conj() if isinstance(M, np.ndarray) else cp.conj(M)
def group_images(a, B=None):
    """The four images g(|a><a| (x) B) as pairs (a', f) with B' = f(B) (B a numpy array or a cvxpy expression)."""
    za = np.array([a[0], -a[1]])
    return [(a, lambda M: M), (za, lambda M: ZZ @ M @ ZZ), (a.conj(), lambda M: SW @ cj(M) @ SW), (za.conj(), lambda M: SW @ cj(ZZ @ M @ ZZ) @ SW)]
def sym_group(M):
    """Average of an 8x8 matrix over {1, Z, T, ZT} (acting as X -> P X P and X -> S X-bar S on three qubits)."""
    S3, P3 = np.kron(np.eye(2), SW), np.kron(np.diag([1., -1]), ZZ)
    X = (M + P3 @ M @ P3) / 2
    return (X + S3 @ X.conj() @ S3) / 2


def merged_images(z):
    """Images of the representative (1, z) grouped by direction: list of (a', [f, ...]); equal directions share one term."""
    groups = []
    for ai, f in group_images(np.array([1, z])):
        for g in groups:
            if abs(abs(np.vdot(g[0], ai)) - np.linalg.norm(g[0]) * np.linalg.norm(ai)) < 1e-9 * np.linalg.norm(ai)**2:
                g[1].append(f); break
        else:
            groups.append((ai, [f]))
    return groups


def tau_expr(zs, Bs):
    out = 0
    for z, B in zip(zs, Bs):
        for ai, fs in merged_images(z):
            out = out + cp.kron(outer(ai), sum(f(B) for f in fs))
    return out


def solve(rho, zs, p_fix=None, kw=KW, fallback=True):
    Bs = [cp.Variable((4, 4), hermitian=True) for _ in zs]
    tau = tau_expr(zs, Bs)
    if p_fix is None:
        p = cp.Variable(); eq = cyc_avg(tau) == (1 - p) * rho + p * np.eye(8) / 8
        prob = cp.Problem(cp.Minimize(p), [B >> 0 for B in Bs] + [eq])
    else:
        p, d = p_fix, cp.Variable()
        eq = cyc_avg(tau) + d * np.eye(8) / 8 == (1 - p) * rho + p * np.eye(8) / 8
        prob = cp.Problem(cp.Maximize(d), [B >> 0 for B in Bs] + [eq, d >= 0])
    try:
        prob.solve(**kw)
    except cp.error.SolverError:                       # CVXOPT stalls at some directions; SCS fallback unless the caller skips
        if not fallback: raise
        prob.solve(solver=cp.SCS, eps=1e-9, max_iters=400000)
    if eq.dual_value is None or prob.value is None or not np.isfinite(prob.value):   # infeasible at the requested p
        return -1.0, np.zeros((8, 8)), [np.zeros((4, 4), complex) for _ in zs]
    W = eq.dual_value; W = cyc_avg(sym_group((W + W.conj().T) / 2))
    return prob.value, W, [B.value for B in Bs]


def qb(th, ph): return np.array([np.cos(th / 2), np.exp(1j * ph) * np.sin(th / 2)])
def lam(W, th, ph):
    a = qb(th, ph); M = np.kron(a.conj()[None, :], np.eye(4)) @ W @ np.kron(a[:, None], np.eye(4))
    return np.linalg.eigvalsh((M + M.conj().T) / 2)[0]
def to_z(th, ph): a = qb(th, ph); return a[1] / a[0]


def same_orbit(z1, z2, tol=0.99999):
    a, b = np.array([1, z1]), np.array([1, z2])
    return any(abs(np.vdot(a, g)) / (np.linalg.norm(a) * np.linalg.norm(g)) > tol for g, _ in group_images(b))


def th_ph(z): a = np.array([1, z]) / np.sqrt(1 + abs(z)**2); return 2 * np.arccos(a[0].real), np.angle(z)


def refine(rho, zs, maxfev=600):
    """Nelder-Mead on the directions (theta_r, phi_r) of the representatives, minimising the SDP value p."""
    def f(x):
        try: return solve(rho, [to_z(x[2 * i], x[2 * i + 1]) for i in range(len(zs))], fallback=False)[0]
        except Exception: return 1.0
    x0 = np.array([c for z in zs for c in th_ph(z)])
    res = minimize(f, x0, method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-13, maxfev=maxfev, adaptive=True))
    print(f'  refine: {len(zs)} reps, p = {res.fun:.9f} ({res.nfev} SDPs)', flush=True)
    return res.fun, [to_z(res.x[2 * i], res.x[2 * i + 1]) for i in range(len(zs))]


# Orbit representatives (theta, phi) of the Bloch directions, read off a first unsymmetrised column-generation run (24, 23 and 10
# near-duplicate directions that cluster at these few orbit points; phi = pi/2 and 0 are the fixed lines of the group action).
SEEDS = {
    'rho_star': [(2.2185, np.pi / 2), (1.1906, np.pi / 2), (1.1070, 0.0)],
    'rho_10/3': [(2.3910, np.pi / 2), (1.5452, np.pi / 2), (1.4050, 0.0)],
    'robust_state': [(2.1989, np.pi / 2), (0.5182, np.pi / 2), (1.6625, 0.7640)],
}


def colgen(rho, name):
    """Move the seed representatives to the optimum (Nelder-Mead on the SDP value), then check by one global search of the dual
    landscape lambda_min(<a|W|a>) that no further direction is violated (optimality of the direction set)."""
    p, zs = refine(rho, [to_z(*t) for t in SEEDS[name]], maxfev=800)
    p, W, Bs = solve(rho, zs)
    W = W if np.trace(W @ rho).real < 0 else -W
    starts = [(t0, f0) for t0 in np.linspace(0.1, 3.0, 8) for f0 in np.linspace(0, np.pi / 2, 5)]
    worst = min(minimize(lambda x: lam(W, *x), s, method='Nelder-Mead', options=dict(xatol=1e-11, fatol=1e-15)).fun for s in starts)
    print(f'  {len(zs)} representatives, p = {p:.9f}, most violated direction {worst:.2e}', flush=True)
    return p, zs


def rnd_c(z, den): return [int(round(z.real * den)), int(round(z.imag * den))]


def gb_ratio(rho, p, zs, Bs):
    """||X - cI||_2 / c for X = (1-p) rho + p I/8 - Cyc(tau), the test of the exact verifier (floating point here)."""
    tau = sum(np.kron(outer(ai), f(B)) for z, B in zip(zs, Bs) for ai, f in group_images(np.array([1, z]), B))
    X = (1 - p) * rho + p * np.eye(8) / 8 - cyc_avg(tau); c = np.trace(X).real / 8
    return np.linalg.norm(X - c * np.eye(8)) / c if c > 0 else np.inf


def build_case(name, rho, level2):
    print(f'{name}: directions', flush=True)
    p_opt, zs = colgen(rho, name)
    pn, pd = math.ceil((p_opt + 1.5e-6) * 10**6), 10**6
    g = math.gcd(pn, pd); pn, pd = pn // g, pd // g
    print(f'{name}: optimal level-two (DPS) white-noise bound p = {level2:.8f} (lower bound); last column-generation p = {p_opt:.9f}; p_up = {pn}/{pd}', flush=True)
    for den_z in (10**4, 10**5, 10**6):
        zr = [complex(*rnd_c(z, den_z)) / den_z for z in zs]
        d, _, Bs = solve(rho, zr, p_fix=pn / pd)
        print(f'  den_z = 1e{round(math.log10(den_z))}: reserved delta = {d:.3e}', flush=True)
        if d <= 0: continue
        keep = [(z, B) for z, B in zip(zr, Bs) if np.trace(B).real > 1e-9]
        for den_b in (10**7, 10**8, 10**9, 10**10):
            zq, Bq = [], []
            for z, B in keep:
                ev, U = np.linalg.eigh((B + B.conj().T) / 2)
                B = (U * np.maximum(ev, 3.0 / den_b)) @ U.conj().T
                Bq.append(np.array([[complex(*rnd_c(B[i, j], den_b)) / den_b for j in range(4)] for i in range(4)])); zq.append(z)
            ratio = gb_ratio(rho, pn / pd, zq, Bq)
            print(f'  den_z = 1e{round(math.log10(den_z))}, den_B = 1e{round(math.log10(den_b))}: {len(zq)} reps, reserved delta = {d:.3e}, GB ratio {ratio:.3f}', flush=True)
            if ratio < 0.5:
                rows = [dict(z=rnd_c(z, den_z), B=[[rnd_c(B[i, j], den_b) for j in range(4)] for i in range(4)]) for z, B in zip(zq, Bq)]
                return dict(name=name, p_up=[pn, pd], den_z=den_z, den_b=den_b, reps=rows)
    raise RuntimeError('no rounding passed')


def main():
    t0 = time.time()
    data = json.load(open(os.path.join(HERE, 'ppt_gme_sym_certificate.json')))
    robust = G.combo(np.array(data['robust']['rho']) / data['denominator'])
    cases = [('rho_star', G.surface(19, 7, 201 / 104)), ('rho_10/3', G.rho_t(10 / 3)), ('robust_state', robust)]
    out = dict(cases=[])
    for name, rho in cases:
        prob, _, _, _ = G.solve_witness(rho, delta=0.0, gamma=0.0)
        out['cases'].append(build_case(name, rho, -prob.value / (1 - prob.value)))
    with open(OUT, 'w') as fh:
        json.dump(out, fh, separators=(',', ':'))
    print(f'wrote {os.path.basename(OUT)}; runtime {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
