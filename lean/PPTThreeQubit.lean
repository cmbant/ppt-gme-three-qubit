module

public import PPTThreeQubit.Basic
public import PPTThreeQubit.Witness
public import PPTThreeQubit.Criterion
public import PPTThreeQubit.FamilyGME
public import PPTThreeQubit.PartialTranspose
public import PPTThreeQubit.FamilyPPT
public import PPTThreeQubit.Coexistence
public import PPTThreeQubit.FlatPauli
public import PPTThreeQubit.Flat
public import PPTThreeQubit.SurfaceWitness
public import PPTThreeQubit.Surface
public import PPTThreeQubit.SurfaceExtreme
public import PPTThreeQubit.SurfacePPT
public import PPTThreeQubit.CoexistencePPT
public import PPTThreeQubit.ProductKernel
public import PPTThreeQubit.CutEquiv
public import PPTThreeQubit.ProductKernelThree
public import PPTThreeQubit.PencilKernel
public import PPTThreeQubit.MomentCone
public import PPTThreeQubit.CurveParam
public import PPTThreeQubit.MomentConeGME
public import FormalResearch.ConvexGeometry.ConeCaratheodory
public import PPTThreeQubit.MomentConeCaratheodory
public import PPTThreeQubit.KernelPenalty
public import PPTThreeQubit.SwapSymmetric
public import PPTThreeQubit.PPTSeparable23
public import PPTThreeQubit.ExtremeCriterion
public import PPTThreeQubit.SigmaExtreme
public import PPTThreeQubit.TraceDistance
public import PPTThreeQubit.Pinching
public import PPTThreeQubit.KernelQuantCuts
public import PPTThreeQubit.ChiralityBasics
public import PPTThreeQubit.ChiralityRays
public import PPTThreeQubit.Chirality
public import PPTThreeQubit.ChiralityStabiliser
public import PPTThreeQubit.ChiralityUnitary
public import PPTThreeQubit.RankProfile
public import PPTThreeQubit.WClass
public import PPTThreeQubit.WClassStates
public import PPTThreeQubit.WClassTilde

/-!
# Fully-PPT three-qubit GME: the rank-five family

Lean companion of `paper/family.tex`
(`prop:ppt3-witness`, `thm:ppt3-family` GME part).  Headlines:
`prodA_classification` / `prodA_of_qs` (`prop:ppt3-witness`(a)), `formH_qs`
(`eq:ppt3-sos`), `formH_nonneg_of_cutProduct` (`prop:ppt3-witness`(b), all cuts),
`trace_Hs_nonneg_of_biseparable` (support criterion), `trace_H5_Xt` (`eq:ppt3-h5-value`),
`ppt3_family_gme` (`thm:ppt3-family`, GME part).

Moment-cone part (`relaxations_and_scope.tex`, `thm:ppt3-universal-moments`; `appendix_moment_cone.tex`,
`lem:ppt3-product-kernel`): `pptSeparable23` (the Horodecki 2⊗3 theorem, transported in
`PPTSeparable23.lean` from the vendored math-lab proof), `ppt_product_kernel`, `kernel_free_of_gme`, `biseparable_iff_mem_cone`, `dual_iff`, `cofCoeff_indep`, `curve_cut`,
`ppt3_universal_moments` (unconditional), `ppt3_universal_moments_of_gme`; all unconditional.

Extremality (`ExtremeCriterion.lean`; `appendix_family.tex`, `appendix_surface.tex`, the supported-face
criterion of Garberg--Irgens--Myrheim used for `cor:ppt3-extreme6` and `cor:ppt3-extreme5`):
`ExtremeCriterion.mem_extremePoints_iff_faceDirection` (generic: a state of the unit-trace set cut out by
finitely many real-linear cut maps is extreme iff every face direction is a multiple of it),
`ExtremeCriterion.mem_extremePoints_of_sandwich` (`σ = B X B†` form), and the three-qubit forms
`ExtremeCriterion.ppt3_extreme_iff`, `ExtremeCriterion.ppt3_extreme_of_sandwich`; `SigmaExtreme.lean` applies it to
the rational surface state (`sigma5_extreme`) with an exact `ℚ(ω)` rank certificate (`SigmaExtremeData.lean`).

Moment cone, Carathéodory and closedness (`ConvexGeometry/ConeCaratheodory.lean`, `MomentConeCaratheodory.lean`;
`thm:ppt3-universal-moments`): `ConvexGeometry.ConeCar.exists_coneRep_card_le` (conic Carathéodory), `ConvexGeometry.ConeCar.isClosed_cone_of_isCompact`,
`biseparable_card_le_25`, `isClosed_inMomentCone`, and for all three-qubit operators `isClosed_isBiseparable`,
`isBiseparable_card_le_64`.

Swap symmetry (`SwapSymmetric.lean`; `family.tex`, `cor:ppt3-symmetrisation`): `isSeparable_of_swap` (PPT and
swap-invariant on `ℂ² ⊗ ℂ⁴` implies separable), `isBiseparable_of_swapBC/AC/AB`, `no_gme_of_swap_invariant` (no
fully-PPT GME state is invariant under the exchange of two qubits), `isBiseparable_swap_average` (pinching),
`isBiseparable_permAvg` (permutation average of a fully-PPT state).

Kernel-penalty lift (`KernelPenalty.lean`; `kernel_penalty_lemma.tex`, `lem:ppt3-kernel-penalty`):
`exists_kernel_penalty` (any compact free set), `exists_kernel_penalty_cutProduct`, `kernel_penalty_witness`
(the nonstrict form with `W₀ + εI + λK†K` and the trace identity on states supported on `ker K`).

Chirality (`ChiralityBasics.lean`, `ChiralityRays.lean`, `Chirality.lean`, `ChiralityStabiliser.lean`;
`family.tex` `cor:ppt3-chirality` (b), proof in `appendix_family.tex`): `det_phiS`, `fullProd_mem_Rs_iff`
(`R_s` contains exactly three fully product rays), `eq_smul_one_of_eigen` (Möbius maps are determined by
three points), `no_local_filter_conj`, `no_local_filter_conj_onto`, `no_local_filter_state`,
`stabiliser_iff`, `exists_filter_perm`, `eq_smul_of_perm`.

Local-unitary chirality, ranks and the mixed `W` class: `ChiralityUnitary.lean` (`ppt3_chirality_conj`,
`ppt3_chirality_timeReverse`: no local unitary relates `ρ_t` to `ρ̄_t` or to its spin time reversal;
`cor:ppt3-chirality` (a)), `RankProfile.lean` (`rankProfile_rhoT_interior` `(5,8,8,8)`,
`rankProfile_rhoT_tMinus`/`tPlus` `(5,7,7,7)`, `rankProfile_rhoSurface` `(5,5,5,5)`, `rankProfile_sigma6`
`(5,6,6,6)`), `WClass.lean`, `WClassStates.lean`, `WClassTilde.lean` (`IsWClass`, `ppt3_w_class`,
`isWClass_rhoTilde`; `thm:ppt3-w-class`).
-/

@[expose] public section

