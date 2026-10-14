# ExaModelsSIF — implemented & verified

Purpose of this file (requested by Sungho, 2026-10-12): record what is already
implemented and HOW it was verified, so verified behavior is not re-checked.
Curate on every change: entries here state the CURRENT truth; re-verify only
when the implementation or the ExaModels dependency changes.

## Dependency (integration stack, 2026-10-13)

ExaModels branch sif-integration (pushed to madsuite-org/ExaModels.jl):
compile-time-prototypes + main (PR #330 bselect/_amul/mod) + add_reduce
(issue #331, Egret) + PR #332 scalar-oracle Hessian fix. Worktree
~/git/worktrees/ExaModels-sif; Sungho merges the PRs on his own schedule.
Wide groups lower through add_reduce with ONE buffered stage per element
TYPE spanning all instances -- stages live in the core tuple type, so the
stage count must never scale with the problem (a per-position layout hung
FMINSURF for 90 min; transposed it builds in 21 s). External-Fortran
problems build via the CUTEst-backed oracle tier; params kwarg implements
sifdecoder -param. Reference adjudications: VANDERM* noise is CUTEst's own
(BigFloat: ours 7e-15, CUTEst 2.8e-6) -- REF_NOISE tol 1e-4.

## Dependency

- ExaModels branch `compile-time-prototypes` (madsuite-org/ExaModels.jl),
  worktree ~/git/worktrees/ExaModels-ctp, branch ctp-plus-330 (compile-time-prototypes merged with main incl. PR #330: bselect/_amul/mod), commit `a683815f6b74cbb2017f862ffd2a4e144484b5aa`.
  Provides `add_expr(...; buffered = true)` (cached subexpressions) used for
  wide-group lowering. Re-pin + re-run pilot when this moves.

## Verified (do not re-check unless the code or dependency changes)

Verification harness: test/pilot.jl — compares nvar/ncon, x0, lvar/uvar,
lcon/ucon, obj, grad, cons, Jacobian, Lagrangian Hessian against CUTEst.jl
(Fortran SIFDecode reference) at x0 and a perturbed point, tol 1e-6.
Population: the 30 PILOT problems in test/pilot.jl (22 original + 8 second
wave), all PASS (2026-10-12, commit of record in git log).

| Feature | Exercised by (within the verified population) |
|---|---|
| Fixed-field card parsing, section dispatch | all 23 |
| DO/OD/ND loops (incl. ND-closes-all, OD-closes-one) | DEGDIAG, GENHS28, HS45, GPP* |
| Integer/real parameter arithmetic (IE/IA/RI/RD/RE/RA/R=/R- etc.) | DEGDIAG, HS45, GENHS28 |
| Array names + inline index substitution | DEGDIAG, GENHS28, GPP* |
| Default bounds [0, inf), FR/MI/PL/UP/LO/FX + X/Z variants | HS5, HS21, ALSOTAME, DEGDIAG (ZL) |
| START POINT incl. 'DEFAULT' | HS35, HS45, GENHS28 |
| Linear groups, CONSTANTS (folded into constraint value), RANGES conv. | TAME, BOOTH, ZANGWIL3, HS35 |
| Group SCALE divides the group value (GF(alpha)/s) | ROSENBR |
| Nontrivial group types on objective and constraints | ROSENBR, TAME, ALSOTAME, HS1/HS2 |
| Elements: types, multiple per group, weights, ZE weight-from-param | HS35, HS14, HS45 |
| Element parameters (EP/P cards) | HS11 |
| Internal variables (IV + R range rows) | HS26, HS28, GENHS28 |
| QUADRATIC section, diagonal halved (:half convention) | DEGDIAG |
| Fortran F-line translation (**, D-exponents, intrinsics, temporaries, GLOBALS) | all nonlinear pilots, ALSOTAME (A-lines) |
| Wide nontrivial groups via buffered subexpression chunks (CHUNK=16) | BROWNAL (199-200 lin terms/group) |
| Element-sum chunking for element-wide nontrivial groups (ELEM_BUDGET=40) | exercised in sweep (MSS3-class) |
| Variable 'SCALE' factors: parsed and IGNORED (LANCELOT hints; CUTEst evaluates unscaled) | MEYER3 |
| Parameter-only conditional temporaries (L temps, A logical, I/E cards) | DECONVC |
| Z variable cards; INTEGER/CONTINUOUS markers ignored (continuous relax.) | CATENA, AVGASB |
| First-named-set rule (CONSTANTS/RANGES/BOUNDS/START POINT) | HS91-class files, DECONVC |
| IR code; RA/RS/RM/RD and IA/IS/IM/ID with parameter operand in field 5 | HS111LNP, DRUGDIS, LUKVLE14 |
| Inline $-comments; params named like (N-2)/3 | GRIDNETG, LUKVLE14 |
| Constraint & variable order: CUTEst reorders BOTH (linear/equation-first rows, nonlinear-first variables). Harness permutes by exact SIFDecode names (varnames/cnames) — never guess the order | HS14, DECONVC, GPP |

*GPP: decode verified (loop semantics); not yet CUTEst-validated.

Calibrated conventions (validated by the population above, do not re-derive):
- Group constant b: constraint value includes -b; bounds built around 0.
- DO-loop variable keeps its LAST value after the loop (loopfinal=:last).
- CUTEst.jl non-default efirst/lfirst flags crash (xlrQh); never use them.

Measured SIF semantics, 2026-10-13 (each validated by the named problem
flipping to PASS against the Fortran reference):
- MI/XM/ZM bound: lower -Inf AND an upper still at the default +Inf drops
  to 0 (MATRIX2).
- A BOUNDS 'DEFAULT' card also governs variables created AFTER the BOUNDS
  section, e.g. by ZV in ELEMENT USES (NUFFIELD; State carries running
  lvar_default/uvar_default used at variable creation).
- A RANGES 'DEFAULT' card leaves equality rows untouched; G/L rows get the
  range (DEGENQP, DEGENQPC).
- ** with an I-declared (integer) temporary exponent is Fortran INTEGER
  power, defined for a negative base: lowered to _npow (registered, zero
  exponent-partials), never _fpow; integer temp assignments truncate
  (MANCINO value to 3e-13 of ref, MANCINONE). _npow propagates non-finite
  exponents (the ExaModels NaNSource structure probe would otherwise throw
  in round(Int, NaN)).
- LOG on numbers AND model expressions is quiet (_flog registered with
  ExaModels; a native log node THROWS DomainError at evaluation): LHAIFAM
  values. DLOG/DSQRT quiet alike.

Hessian reference bugs (hand-coded H lines contradict fd of the reference's
own gradient while ours matches; HESS_REF_BUGS gate in pilot.jl): GULF,
WATSON, HIMMELBB, HIMMELBF, MAXLIKA, ECKERLE4LS, HS70, ODC, and (y-weighted
Lagrangian fd, 2026-10-13) GULFNE, HIMMELBFNE, WATSONNE, ROBOTARM (fd equals
OUR worst entry to 12 digits, 18.236447894217633 vs ref 18.2828).

Microbenchmark (this host, 002): single traced kernel compile time vs linear
width L: 16→3.8s, 32→17.7s, 48→53s (superquadratic). Buffered chunks (C=16):
L=96→5.8s, L=200→16s, L=400→59s; obj/grad exact, Hessian to fd-noise (1e-8).

## Known unsupported (sweep buckets, to be fixed; see sweep_results.tsv)

- Conditional expressions in A-lines (.LE. etc.)
- Integer-variable markers
- MOD intrinsic on model expressions (ExaModels lacks mod node)
- Element part of a nontrivial group wider than 64 traced fields
- Groups named with embedded $ comments ("$ flow to" bucket) — parse issue

## Environment

julia --project=. always; Manifest.toml committed. Canonical SIF collection:
the CUTEst_jll artifact (1548 problems) — exactly what the Fortran reference
decodes. The old ~/git/pkg/cutest-sif mirror (lpoo/cutest-sif, 1149 files) is
stale: 390 files differ, 411 problems missing, 12 removed upstream. Progress
denominators are /1548 from round 4 on. Full-collection sweep:
bash test/sweep_driver.sh N (restartable, shardable, resumes from
sweep_results*.tsv; SWEEP_ONLY=<file> restricts to listed names).
