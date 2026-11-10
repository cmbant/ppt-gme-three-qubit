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
public import PPTThreeQubit.PPTSeparable23

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
-/

@[expose] public section

