# Fully-PPT three-qubit genuine entanglement: companion material

**Compiled PDF:** [paper/fully_ppt_three_qubit_gme.pdf](paper/fully_ppt_three_qubit_gme.pdf)

The public repository is [cmbant/ppt-gme-three-qubit](https://github.com/cmbant/ppt-gme-three-qubit).

[![DOI](https://zenodo.org/badge/1452542690.svg)](https://doi.org/10.5281/zenodo.24113583)

This is the companion to the paper *Genuinely multipartite bound entanglement of three qubits*. The paper gives the constructions and proofs; this tree also holds the exact certificates and the Lean
formalisation that back them.

## What the paper shows

- A mixture of five orthogonal GHZ/W-type vectors gives a one-parameter family of rank-five three-qubit states that is
  positive under partial transposition (PPT) across all three cuts on an exact interval, and genuinely multipartite
  entangled (GME) throughout it. To our knowledge no fully-PPT GME three-qubit state was known; rank five is the minimum possible.
- One Hermitian form on the five-dimensional support, nonnegative on its product-vector curves by an explicit square
  factorisation, excludes biseparable mixtures across different cuts.
- Extreme fully-PPT GME states with partial-transpose rank profiles `(5,6,6,6)` and `(5,5,5,5)` exist, with an analytic
  surface of states covering the whole PPT window for the latter.
- Every explicit state, including the two extreme states and the noise-robust full-rank state, is a mixture of pure states
  with zero three-tangle, so it lies in the mixed W class, as Acín et al. conjectured for PPT entangled three-qubit states.
- One rank-five support carries both biseparable and GME fully-PPT states, separated by an exact transition inside the
  strictly PPT region.
- Exact-rational symmetric-extension (level-two) witnesses and explicit biseparable decompositions bound the white-noise
  thresholds from both sides: an explicit full-rank fully-PPT state stays GME under 2.33% white noise, while
  2.336% destroys it (known to within half a per cent in relative terms).

## Layout

| folder | contents |
|---|---|
| `paper/` | `fully_ppt_three_qubit_gme.pdf` and its TeX sources (the master and its part files); the PDF is built from this copy |
| `verifiers/` | `code/`: the exact verifiers named in the paper's verification appendix, the helper modules they import and the certificate data (JSON); `verification/`: their captured output; `MAP.md`: check tag, script, log; `requirements.txt`: pinned versions |
| `lean/` | a standalone Lean 4 / Mathlib formalisation of the main exact statements (below) |

## Where each result is checked

The table lists each theorem-like statement of the paper with its number and page, its claim class and the tags of the
checks it rests on. Each tag is a row of [`verifiers/MAP.md`](verifiers/MAP.md), which gives the script that backs it and
says what it checks. A dash in the checks column means the proof is analytic and has no finite check. "Proved" means proved
on paper; the exact-certificate classes mean a finite exact computation (rational or algebraic arithmetic, no floating
point) certifies the statement; "numerical evidence" is not a proof.

| result | TeX label | statement | claim class | checks (tags of `verifiers/MAP.md`) |
|---|---|---|---|---|
| Theorem 1, p. 3 | `thm:ppt3-family` | the rank-five fully-PPT GME family | proved | CF-1, CF-2, CF-3, IN-1, IN-2, IN-3, IN-4, IN-5, IN-6, IN-7 |
| Proposition 2, p. 4 | `prop:ppt3-witness` | supported witnesses | proved analytically | CF-2, CF-3, IN-2, IN-3, IN-4, SW-1, SW-2 |
| Corollary 3, p. 5 | `cor:ppt3-symmetrisation` | reflection symmetry destroys the construction | proved | – |
| Corollary 4, p. 5 | `cor:ppt3-chirality` | intrinsic local chirality | proved | SC-1, SC-2, SC-3, SC-4, SC-5 |
| Corollary 5, p. 6 | `cor:ppt3-extreme6` | an extreme fully-PPT GME state of profile `(5,6,6,6)` | exact-symbolic certificate, computer-assisted | CF-6 |
| Lemma 6, p. 7 | `lem:ppt3-h0` | a `\tau`-independent supported witness | proved analytically | SW-7 |
| Theorem 7, p. 7 | `thm:ppt3-orbit` | fully-PPT GME states of profile `(5,5,5,5)` | proved analytically | FP-1, SW-4, SW-5, SW-6, SW-7 |
| Corollary 8, p. 7 | `cor:ppt3-extreme5` | an explicit extreme state with rational weights and profile `(5,5,5,5)` | exact-rational certificate, computer-assisted | FP-2, FP-3, FP-4, FP-5, SW-3, DS-6, SB-1, SB-2 |
| Proposition 9, p. 8 | `prop:ppt3-moduli` | the profile-`(5,5,5,5)` states are not isolated | exact-symbolic certificate, computer-assisted | MO-1, MO-2, MO-3, MO-4, SB-1, SB-3, SB-4 |
| Theorem 12, p. 9 | `thm:ppt3-w-class` | the constructed fully-PPT states belong to the mixed `W` class | proved analytically | WC-1, WC-2, WC-3, WC-4 |
| Theorem 14, p. 12 | `thm:ppt3-phase-coexistence` | biseparable and GME states on one PPT support | proved | PG-2 |
| Proposition 15, p. 12 | `prop:ppt3-exposed-transition` | an exposed phase boundary inside the strict-PPT region | proved | PG-4 |
| Corollary 16, p. 12 | `cor:ppt3-volume` | positive volume in any finite dimension | proved | – |
| Lemma 17, p. 13 | `lem:ppt3-dps` | certificate and stability principle | proved | – |
| Corollary 18, p. 13 | `cor:ppt3-noise` | noise and preparation radius of the family | exact-rational certificate, with proved PPT stability | DS-1, DS-2, DS-3, DS-6 |
| Corollary 19, p. 14 | `cor:ppt3-noise-reused` | noise bounds for the surface state from one exact witness | exact-symbolic certificate | DS-6 |
| Corollary 20, p. 14 | `cor:ppt3-robust` | a full-rank fully-PPT GME state tolerating `2.3\%` white noise | exact-rational certificate, with proved PPT stability | DS-4, DS-6 |
| Proposition 21, p. 15 | `prop:ppt3-noise-sharp` | the white-noise thresholds are essentially sharp | exact-rational/algebraic certificate, computer-assisted | DS-7 |
| Proposition 22, p. 15 | `prop:ppt3-noise-extremal` | an extremal reduction for white-noise robustness | proved | – |
| Corollary 24, p. 16 | `cor:ppt3-ppt-mixture` | PPT mixtures and outer relaxations | proved | – |
| Theorem 26, p. 16 | `thm:ppt3-universal-moments` | three images of one quartic moment cone | proved | – |
| Conjecture 27, p. 17 | `conj:ppt3-flat` | flat rank-five states | open | – |
| Lemma 29, p. 28 | `lem:ppt3-product-kernel` | product kernel | proved | – |
| Lemma 30, p. 29 | `lem:ppt3-kernel-penalty` | kernel-penalty witness lift | proved | – |
| Proposition 31, p. 29 | `prop:ppt3-kernel-quant` | quantitative form | proved | – |


To check one result, run the scripts of every tag in its row. Appendix G of the paper lists the
proof dependencies by group.

## Running the verifiers

```bash
python3 -m pip install -r verifiers/requirements.txt
cd verifiers/code
python3 verify_ppt_gme_moduli.py     # exits 0; output equals ../verification/verify_ppt_gme_moduli.txt
```

Each `verify_*` and `certify_*` script prints its checks and exits with status 0 if all pass; its output is committed as
`verification/<script>.txt` (timing lines may differ). The scripts use exact arithmetic (`sympy`, rational and algebraic
numbers) and read the certificate data in `code/` by relative path, so run them from `verifiers/code/`. The `explore_*`
scripts are numerical evidence only (they need `cvxpy` and a conic solver), search for the certificates and rewrite the
data files they produce: run them on a copy of `code/`. The versions in `requirements.txt` are those the logs were
captured with.

## Claim classes

Every statement in the paper is proved unless tagged. The classes are: **proved**; **proved under explicit hypothesis**;
**exact-symbolic / exact-rational certificate** (a finite exact computation certifies a finite statement; it does not by
itself prove a statement for all parameters); **numerical evidence**; **open / conjectural**. A script with no committed
log is unrun, and a claim resting on it would be at most numerical evidence.

## Lean

```bash
cd lean
lake exe cache get   # downloads the Mathlib build cache (network)
lake build
```

The project has two libraries: the paper's development `lean/PPTThreeQubit/` and the general Gurvits–Barnum separable ball
`lean/GurvitsBarnum/`, each with an aggregate (`lean/PPTThreeQubit.lean`, `lean/GurvitsBarnum.lean`); the modules are
`PPTThreeQubit.Basic`, `PPTThreeQubit.Chirality`, `PPTThreeQubit.ChiralityBasics`, `PPTThreeQubit.ChiralityRays`, `PPTThreeQubit.ChiralityStabiliser`, `PPTThreeQubit.ChiralityUnitary`, `PPTThreeQubit.Coexistence`, `PPTThreeQubit.CoexistencePPT`, `PPTThreeQubit.CoexistencePPTData`, `PPTThreeQubit.Criterion`, `PPTThreeQubit.CurveParam`, `PPTThreeQubit.CutEquiv`, `PPTThreeQubit.ExtremeCriterion`, `PPTThreeQubit.FamilyGME`, `PPTThreeQubit.FamilyPPT`, `PPTThreeQubit.FamilyPPTFactor`, `PPTThreeQubit.Flat`, `PPTThreeQubit.FlatPauli`, `PPTThreeQubit.KernelPenalty`, `PPTThreeQubit.KernelQuantCuts`, `PPTThreeQubit.MomentCone`, `PPTThreeQubit.MomentConeCaratheodory`, `PPTThreeQubit.MomentConeGME`, `PPTThreeQubit.PPTSeparable23`, `PPTThreeQubit.PartialTranspose`, `PPTThreeQubit.PencilKernel`, `PPTThreeQubit.Pinching`, `PPTThreeQubit.ProductKernel`, `PPTThreeQubit.ProductKernelThree`, `PPTThreeQubit.RankProfile`, `PPTThreeQubit.SigmaExtreme`, `PPTThreeQubit.SigmaExtremeAlg`, `PPTThreeQubit.SigmaExtremeData`, `PPTThreeQubit.Surface`, `PPTThreeQubit.SurfaceExtreme`, `PPTThreeQubit.SurfacePPT`, `PPTThreeQubit.SurfacePPTFactor`, `PPTThreeQubit.SurfaceWitness`, `PPTThreeQubit.SwapSymmetric`, `PPTThreeQubit.TraceDistance`, `PPTThreeQubit.WClass`, `PPTThreeQubit.WClassStates`, `PPTThreeQubit.WClassTilde`, `PPTThreeQubit.Witness`, `GurvitsBarnum.Ball`, `GurvitsBarnum.Basic`, `GurvitsBarnum.BlockBound`, `GurvitsBarnum.Choi`, `GurvitsBarnum.Closed`, `GurvitsBarnum.Contraction`, `GurvitsBarnum.ContractionBase`, `GurvitsBarnum.ContractionPolar`, `GurvitsBarnum.ContractionUnitary`, `GurvitsBarnum.Separation`, `GurvitsBarnum.Witness`, `FormalResearch.ConvexGeometry.ClosedHull`, `FormalResearch.ConvexGeometry.ConeCaratheodory`, `FormalResearch.MatrixAnalysis.ExtremeCriterion`, `FormalResearch.MatrixAnalysis.HelstromPure`, `FormalResearch.MatrixAnalysis.KernelQuant`, `FormalResearch.MatrixAnalysis.TraceDistance`, `FormalResearch.MatrixAnalysis.TraceNorm`, `FormalResearch.MatrixAnalysis.TraceNormGeneral`, `FormalResearch.MatrixAnalysis.VectorNorm`. The
toolchain is pinned in `lean-toolchain` and Mathlib in `lakefile.toml` and `lake-manifest.json`. Each module's docstring
names the paper label it formalises (`grep` it in `paper/*.tex`). Biseparability is the finite-decomposition definition (a
sum of rank-one projectors onto vectors that are products across some cut) and GME is its negation. The headline theorems
depend only on `propext`, `Classical.choice` and `Quot.sound` (pinned `#print axioms` guards); there is no `sorry` and no
`native_decide`. The data module `CoexistencePPTData.lean` is generated by `PPTThreeQubit/gen_coex_ppt.py` (exact
arithmetic over `Q(omega)`, `omega^2 = -1 - omega`) and checked by the kernel.

| paper label | Lean declarations | scope |
|---|---|---|
| `prop:ppt3-witness` | `prodA_classification`, `prodA_of_qs`, `prodA_of_v`, `formH_qs`, `formH_nonneg_of_cutProduct` | every `s > 0`, all three cuts |
| support argument (proof of GME) | `trace_Hs_nonneg_of_biseparable`, `trace_nonneg_of_biseparable` | any supported form nonnegative on cut products |
| `thm:ppt3-family` | `ppt3_family_ppt_iff`, `ppt3_family_ppt_gme`, `ppt3_family_gme`, `tInterval_order`, `ppt3_family_explicit_member_ppt` | exact PPT window `[t_-, t_+]` (both directions) and GME; ranks: `rank_rhoT` (`rank ρ_t = 5`), `rankProfile_rhoT_interior` (`(5,8,8,8)`), `rankProfile_rhoT_tMinus`, `rankProfile_rhoT_tPlus` (`(5,7,7,7)`); minimality of rank five not formalised |
| `thm:ppt3-phase-coexistence` | `ppt3_phase_coexistence`, `ppt3_phase_coexistence_full`, `sigmaBS_strict_ppt` | explicit biseparable decomposition; strict PPT |
| `prop:ppt3-exposed-transition` | `ppt3_exposed_state`, `ppt3_exposed_transition`, `ppt3_exposed_transition_full`, `exposed_line_fully_ppt` | strict PPT along the whole line |
| `lem:ppt3-h0` | `formH0_qs_identity`, `formH0_nonneg_of_cutProduct`, `trace_H0_nonneg_of_biseparable` | |
| `thm:ppt3-orbit` (b), PPT and GME | `rhoSurface_filter`, `ppt3_orbit_ppt`, `ppt3_orbit_gme`, `ppt3_orbit_ppt_gme`, `rankProfile_rhoSurface` | every `1/sqrt 3 < tau < 1`, `a > 0`, with the rank profile `(5,5,5,5)`; part (c) not formalised |
| `cor:ppt3-extreme5` (PPT, GME) | `sigma5_eq_member`, `sigma5_gme`, `sigma5_ppt_gme`, `rankProfile_sigma5` | the level-two noise bound not formalised |
| `cor:ppt3-extreme5`, extremality of the rational member | `sigma5_extreme` (via `ppt3_extreme_of_sandwich`), `sigma5_null`, `face_solve`, `SigmaExtremeData.cert`, `SigmaExtremeData.kernel_check` | the single state `tau = 48/55`, `Q = 11/5`: the 24 face equations are solved by an exact certificate over `Q(omega)` checked by the kernel. Extremality for all `Q`, and of `sigma_6` (ranks `rankProfile_sigma6`), is not formalised |
| supported-face criterion (proofs of `cor:ppt3-extreme6`, `cor:ppt3-extreme5`, `thm:ppt3-orbit` extremality) | `mem_extremePoints_iff_faceDirection`, `mem_extremePoints_of_sandwich`, `ppt3_extreme_iff`, `ppt3_extreme_of_sandwich` | a state of the unit-trace set cut out by finitely many real-linear cut maps is extreme iff every face direction is a multiple of it; three-qubit form with `ptA/ptB/ptC` |
| `cor:ppt3-chirality` (a), (b) | `ppt3_chirality_conj`, `ppt3_chirality_timeReverse`, `ppt3_chirality_window`; `fullProd_mem_Rs_iff`, `no_local_filter_conj`, `no_local_filter_state`, `stabiliser_iff` | (a): no local unitary maps `rho_t` to its conjugate, for every `t > 0`; (b): `R_s` has exactly three fully product rays and no invertible local filter relates it to its conjugate; the abstract isomorphism with `S_3` is not formalised |
| `cor:ppt3-symmetrisation` | `isSeparable_of_swap`, `no_gme_of_swap_invariant`, `isBiseparable_swap_average`, `isBiseparable_permAvg`; `bsDist_le_permAvg`, `bsDist_le_cyclic`, `bsDist_le_swap_average`, `bsDist_le_tnormG_P3P1` | PPT and swap-invariant implies separable; the distance bounds `eq:ppt3-asymmetry-upper` and `eq:ppt3-pinch` with the trace norm `tnormG X = sup {Re tr(UX) : U†U ≤ 1}` |
| `lem:ppt3-kernel-penalty`, `prop:ppt3-kernel-quant` | `exists_kernel_penalty`, `kernel_penalty_witness`; `kernel_quant_form`, `kernel_quant_matrix`, `kernel_quant_cutProduct`, `kernel_quant_distance`, `kernel_lift_distance` | the lifted witness is non-negative on every cut-product vector for `lambda = (L + L^2/alpha)/mu`, and `D(rho, BS) >= g/(2L + lambda beta)`; the explicit constants of `rem:ppt3-hand-radius` are not formalised |
| `thm:ppt3-w-class` (mixed `W` class) | `hyperdet`, `IsWClass`, `isWClass_rhoT`, `isWClass_sigma6`, `isWClass_rhoSurface`, `isWClass_sigma5`, `isWClass_sigmaLine`, `isWClass_rhoTilde`, `ppt3_w_class` | the cone of non-negative combinations of `|phi><phi|` with vanishing hyperdeterminant; membership of the unnormalised operators, diagonal local filters only; `rem:ppt3-w-class-scope` is not formalised |
| moment cone: Carathéodory and closedness | `ConvexGeometry.ConeCar.exists_coneRep_card_le`, `biseparable_card_le_25`, `isClosed_inMomentCone`, `isClosed_isBiseparable`, `isBiseparable_card_le_64` | conic Carathéodory bound; the 25-term bound given `CutCurve` with invertible `T_c` |
| flat states fully PPT (Section 7.3) | `Flat.flat_form`, `Flat.flat_pt_conj`, `Flat.flat_fully_ppt` | flatness defined by maximally mixed marginals |
| local filters preserve biseparability | `isBiseparable_conj_localDiag_iff`, `isBiseparable_smul_pos` | diagonal product filters |
| `lem:ppt3-product-kernel` | `ppt_product_kernel`, `kernel_free_of_gme` | unconditional; the 2⊗3 PPT separability theorem (Horodecki) it uses is `pptSeparable23`, transported from the vendored math-lab proof (`lean/ThirdParty/math-lab`) |
| `thm:ppt3-universal-moments` | `ppt3_universal_moments`, `biseparable_iff_mem_cone`, `dual_iff`, `ppt3_universal_moments_of_gme` | unconditional (the 25-term bound is in the Carathéodory row) |
| Gurvits–Barnum ball (PRA 66, 062311 (2002); used in `prop:ppt3-noise-sharp`) | `gurvits_barnum`, `isSeparable_of_sum_norm_sq_sub_smul_one_le`, `sum_norm_sq_le_trace_sq_of_blockPositive` | any finite bipartite dimensions; Hermitian `Δ` with `‖Δ‖₂ ≤ 1` gives `I + Δ` separable (an explicit finite sum of product projectors); the dual inequality `‖W‖₂ ≤ tr W` for block-positive `W` |

Not formalised: extremality beyond the single rational state `sigma` of `cor:ppt3-extreme5` (all `Q`; `sigma_6`), which in the paper rests on the `2x4` low-rank theorem
or the exact verifiers; minimality of rank five; the level-two noise certificates; the explicit biseparable decompositions behind the
white-noise upper bounds (the Gurvits–Barnum step they end with is formalised); `rem:ppt3-hand-radius`; the endpoint remark. These rest on the exact
verifiers above.
A `.lean` file is not a verification until it builds: build the project before relying on it.

## Rebuilding the PDF

```bash
cd paper
pdflatex fully_ppt_three_qubit_gme.tex
pdflatex fully_ppt_three_qubit_gme.tex
```

## Licence

The code in `verifiers/` and `lean/` is MIT-licensed, except `lean/ThirdParty/math-lab/`, a vendored Apache-2.0 proof of
the 2⊗3 PPT criterion that keeps its own licence and notice files (`lean/ThirdParty/math-lab/NOTICE.md` lists its source
and our two compatibility patches). The paper in `paper/` is CC BY 4.0. See `LICENSE`.
