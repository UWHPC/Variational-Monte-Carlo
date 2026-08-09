# Variational Monte Carlo correctness and backend plan

## Purpose

This document is the implementation plan for turning this repository into a
scientifically defensible, polished VMC code with genuinely separate CPU and
CUDA backends.

The intended runtime contract is:

```text
vmc --config run.cfg           # complete simulation on host
vmc --cuda --config run.cfg    # complete simulation on the CUDA device
```

For the CUDA backend, “complete simulation on device” means that configuration
and immutable tables may be transferred to the GPU once, and final summaries or
explicit checkpoints may be transferred back, but initialization of walker
state, warmup, proposals, accept/reject decisions, wavefunction updates, local
energies, adaptation, and statistical accumulation remain on the device. There
must be no per-move host dereferences, scalar copies, or host/device state
ping-pong.

The plan is based on an audit performed on 2026-08-08 against repository commit
`cd3558e` plus the current uncommitted xpu migration, and against xpu commit
`5ab6884`. No implementation task below should be marked complete merely because
the stock test suite passes; each item has its own acceptance criteria.

## Evidence standard and terminology

Every finding below is classified implicitly by the evidence stated with it:

- **Confirmed source defect** means the result follows directly from the code,
  an algebraic comparison, or a compiler diagnostic.
- **Measured defect or risk** means a dedicated diagnostic reproduced it and
  the observed values are recorded here.
- **Design requirement** means it follows from the requested CPU/GPU contract or
  xpu's documented memory model.
- **Proposed fix** is a recommended implementation. It is not treated as proven
  until its listed acceptance criteria pass.
- **Not yet quantified** identifies questions for which the audit found a real
  mathematical or engineering concern but did not establish its effect on the
  final energy.

This distinction is intentional. In particular, no current VMC CUDA result has
been declared correct: the current CUDA target does not compile, so it could not
be run end to end.

## Verified baseline

### Repository and build state

- The current worktree contains an unfinished migration from the repository's
  old allocation and view types to xpu. Existing user changes are intentionally
  not enumerated as completed work in this plan.
- The clean reference's `src/utilities/aligned_soa.cuh` allocates CUDA storage
  with `cudaMallocManaged`. The current migration deletes that utility; xpu
  commit `5ab6884` uses ordinary `cudaMalloc`, which is aligned with the target
  device-residency contract.
- A fresh CPU build of the current tree with GCC 15 and FP64 fails to compile.
  Confirmed failures include use of `Particles::size()` after the API changed to
  `count()`, use of obsolete `.x_/.y_/.z_` members on `xpu::soa_view`, and an
  include of deleted `src/utilities/ptr3d.hpp`.
- A fresh CUDA 13.3 build for compute capability 8.6 also fails. In addition to
  the preceding migration errors, VMC calls `xpu::cuda_check`, while xpu commit
  `5ab6884` defines `xpu::cu_check` in `xpu/config.hpp`.
- The backend is currently selected at configure time by `VMC_ENABLE_CUDA` in
  `CMakeLists.txt`; `src/main.cu` ignores `argc` and `argv` and always loads
  `config.cfg`. There is no runtime `--cuda` interface.

### Tests and diagnostics

The following results were reproduced in clean temporary build directories:

| Artifact | Configuration | Result |
|---|---|---|
| Clean VMC commit `cd3558e` | CPU, FP64, Debug | 94/94 stock tests pass |
| Clean VMC commit `cd3558e` | CPU, FP32, Debug | 94/94 stock tests pass |
| Clean VMC commit `cd3558e` | CPU, FP64, ASan+UBSan | 94/94 stock tests pass |
| Clean VMC commit `cd3558e` | CPU, FP64, strict GCC 15 warning set | builds with zero project warnings |
| Clean VMC commit `cd3558e` | CUDA 13.3, FP64, SM 8.6 | builds, with 16 host/device-call diagnostics from `metropolis_step` |
| Current worktree | CPU, FP64, Debug | compilation fails |
| Current worktree | CUDA 13.3, FP64, SM 8.6 | compilation fails |
| xpu `5ab6884` | CPU smoke test | passes |
| xpu `5ab6884` | CUDA smoke test on RTX 3070 | passes |
| xpu `5ab6884` | CUDA compute-sanitizer memcheck | 0 reported errors |

LeakSanitizer was disabled for the sanitizer run because it could not operate
under the traced Catch2 discovery environment. Leak freedom is therefore not a
verified property.

### Physical model audited

The audited model is a fully spin-polarized three-dimensional homogeneous
electron gas in a cubic periodic cell, in Hartree atomic units, using:

- a real plane-wave Slater determinant;
- a two-body Padé Jastrow factor;
- single-electron random-scan Metropolis moves;
- a neutralizing-background Ewald interaction; and
- the local-energy identity
  \(E_L=-\frac12\sum_i[\nabla_i^2\log\Psi+|\nabla_i\log\Psi|^2]+V\).

Generalization to an unpolarized gas, arbitrary spin populations, twists,
backflow, or non-cubic cells is not part of the current verified model and must
not be implied by documentation or output.

## Priority and phase definitions

- **P0 — stop-the-line:** current results or builds cannot be trusted until the
  item is resolved.
- **P1 — required for a credible release:** necessary for the requested backend
  contract, reproducibility, or scientific validation.
- **P2 — polish and performance:** important for a codebase suitable for public
  use, after correctness gates pass.

Priority is scoped to the affected product or release gate; it is not a command
to serialize every P0 item before all P1 work. For example, a CUDA-residency P0
does not block independent CPU-oracle work, and an external-validation P1 does
not block construction of the GPU ownership boundary. The phases group related
work, but the executable schedule and actual dependency edges are defined near
the end of this document. Measurement infrastructure and an unoptimized baseline
must be established early; performance-changing optimization remains behind the
scientific correctness gates.

### Decisions that must be made explicitly

The audit did not invent values for requirements that have not yet been stated.
The following decisions should be recorded before their dependent release gate:

1. **Accuracy target:** choose total- and per-electron energy tolerances that
   determine Ewald convergence, precision policy, run length, and acceptable
   CPU/CUDA differences. Until then, convergence data must be reported instead
   of declaring a cutoff “accurate enough.”
2. **Precision policy:** decide whether FP32 is a production mode, a mixed-
   precision performance mode with FP64 critical operations, or test-only. The
   measured determinant drift means FP32 cannot be accepted solely because the
   short stock suite passes.
3. **Trial ansatz:** approve the periodic short-range cutoff and whether a
   long-range reciprocal-space Jastrow is required for the intended HEG claims.
4. **Finite-size scope:** decide whether the product will implement twists and
   finite-size corrections or explicitly limit itself to finite-cell results.
5. **Trajectory policy:** decide whether production CUDA runs require particle
   trajectories. If so, specify a checkpoint cadence and accepted transfer cost;
   per-move trajectory output is incompatible with the residency goal.
6. **GPU mapping:** treat one block per walker as the initial policy to test, not
   a promised optimum. Select the final mapping from correctness-preserving
   benchmarks.

---

## Phase 0 — Freeze the reference and restore buildability

### BUILD-001 — Record and isolate the in-progress xpu migration

**Priority:** P0

**Verified problem:** The current tree mixes deleted VMC utilities, old xpu-like
accessors, and the current xpu API. Both CPU and CUDA builds stop before tests can
be compiled.

`ALIGNED_SOA_REMOVAL.md` is a useful mechanical inventory, but it is not a valid
description of the target CUDA architecture. Its section 7 proposes explicit
device-to-host copies for scalar results, including the determinant ratio and
Jastrow delta on every Metropolis move. Those copies would avoid illegal host
dereferences of `cudaMalloc` storage, but they would preserve a host-driven move
loop and directly violate the residency contract in this plan.

**Why this matters:** No subsequent physics or backend change can be validated
against a tree that does not build. Combining API migration, scientific changes,
and GPU restructuring in one unreviewable patch also makes regressions difficult
to localize.

**Proposed fix:**

1. Treat `cd3558e` as the temporary behavioral reference, not as a declaration
   that all of its science is correct.
2. Complete the xpu API migration mechanically before redesigning algorithms:
   replace `size()` with the selected particle-count API; replace named SoA
   fields with indexed or structured component access; remove deleted includes
   and `AlignedSoA`; and rename `xpu::cuda_check` to `xpu::cu_check`.
3. Reconcile or supersede `ALIGNED_SOA_REMOVAL.md`: retain its API-call-site
   inventory, but mark its proposed per-move scalar readbacks as rejected target
   behavior. Initialization and explicit result/checkpoint transfers belong in
   the CUDA backend boundary; move-level ratios and decisions must remain on the
   device. A readback used only by a temporary diagnostic must be isolated from
   production and carry an explicit deletion condition.
4. Do not add compatibility aliases for deleted interfaces unless they are
   deliberately short-lived and carry a removal issue. Compatibility aliases
   would conceal unfinished migration work.
5. Split mechanical migration commits from scientific corrections. Until the
   device engine reaches Gate C, label a build-only CUDA artifact as
   non-production; compilation alone does not validate its address-space use.

**Done when:**

- [ ] A clean current-tree CPU configure, build, and test run succeeds.
- [ ] A clean current-tree CUDA configure and build succeeds.
- [ ] `rg` finds no use of `AlignedSoA`, `ptr3d.hpp`, obsolete SoA named fields,
      or `xpu::cuda_check`.
- [ ] Migration documentation no longer recommends per-move D2H transfers as a
      production solution, and any diagnostic-only readback is visibly tracked
      for removal.
- [ ] The build does not rely on stale objects from an existing build directory.

### BUILD-002 — Establish a reproducible build matrix

**Priority:** P1

**Verified problem:** The README advertises CMake 3.20, GCC 13, Clang, native
Windows, and macOS. `CMakeLists.txt` actually requires CMake 3.25 and GNU GCC 15,
and explicitly rejects native Windows; no macOS configuration is exercised by
CI. The only build/test workflow uses one unspecified `ubuntu-latest` Release
job. Release flags include unconditional `-march=native`, and CUDA defaults
to `CMAKE_CUDA_ARCHITECTURES=native`, so two nominally identical Release builds
can target different instruction sets without that difference appearing in the
current result metadata.

**Why this matters:** A scientific result must be reproducible from documented
requirements. A build that silently selects CUDA when a compiler happens to be
installed also does not provide the requested explicit backend contract.

**Proposed fix:**

1. Add CMake presets for at least `cpu-debug-fp64`, `cpu-release-fp64`,
   `cpu-debug-fp32`, `cuda-debug-fp64`, and `cuda-release-fp64`.
2. Either support and test Clang on CPU or state that GCC 15 is the sole supported
   compiler. Do not claim untested platforms.
3. Make CUDA artifact construction explicit. Absence of nvcc must not transform
   a requested CUDA build into a CPU build that appears successful.
4. Keep the xpu commit pinned and print the VMC/xpu revisions in configure and run
   metadata.
5. Provide a portable release baseline and a separately named native-tuned
   preset. Make CPU ISA and CUDA architecture lists explicit and record the
   effective compiler/link flags in build provenance.
6. Lock FetchContent inputs to reviewed immutable revisions (prefer full commit
   IDs over movable tags/short identifiers), and document a cached/offline
   dependency path for reproducible builds.

**Done when:**

- [ ] Every documented preset is exercised in CI or explicitly documented as a
      local-only/self-hosted GPU preset.
- [ ] README requirements exactly match enforced CMake requirements.
- [ ] A requested CUDA build fails clearly if its toolchain is unavailable.
- [ ] A build manifest distinguishes portable and native-tuned artifacts and
      records their CPU/CUDA architecture targets.
- [ ] Dependency resolution records the exact fetched revisions and can be
      reproduced without depending on a tag retaining its old target.

### DEP-001 — Harden and pin the xpu dependency contract

**Priority:** P1

**Verified status:** xpu commit `5ab6884` uses aligned `operator new` on CPU and
`cudaMalloc` on CUDA. Its CPU and real-GPU smoke tests pass, and its CUDA smoke
test reports zero errors under compute-sanitizer memcheck. This verifies the
tested buffer/SoA allocation, fill, launch, and device-to-host copy path; it does
not verify VMC or all xpu APIs.

The same xpu README states that CPU and CUDA cannot be mixed in one linked
program, transfers are not wrapped, accessors are not bounds checked, and
multidimensional launch configuration is still work in progress. Its testing
section is stale: it advertises `--cpu`, while `scripts/test.sh` accepts `--cpp`,
and it says the CPU backend lacks a runtime test even though `tests/smoke.cpp` is
built and run in the CPU configuration.

**Why this matters:** VMC should depend on an explicit, tested xpu contract rather
than compensating for evolving APIs throughout physics code. Stale dependency
documentation can also cause a CPU smoke test to be skipped while appearing to
follow documented instructions.

**Proposed fix:**

1. Keep xpu pinned to a reviewed commit and update it only in a dedicated change
   with both xpu and VMC test matrices.
2. Fix xpu's testing documentation or script flags and accurately state the
   scope of its CPU/CUDA coverage.
3. Add focused xpu tests needed by VMC: zero-sized and ordinary buffers, SoA
   indexing/stride, move ownership, initialization, overflow rejection, CUDA
   error reporting, and supported launch dimensions.
4. Decide whether allocation/CUDA errors should abort, throw, or return an error.
   If xpu retains abort semantics, VMC must preflight capacities and document
   which unexpected dependency failures remain fatal.
5. Keep explicit transfer operations inside the CUDA backend's initialization,
   checkpoint, and result layers; do not scatter raw `cudaMemcpy` calls through
   physics components merely because xpu does not wrap transfers.

**Done when:**

- [ ] xpu's documented commands match its actual script and CI.
- [ ] The pinned xpu revision passes its CPU, CUDA, and compute-sanitizer tests.
- [ ] VMC has one narrow adapter/ownership layer for any xpu API expected to
      evolve.
- [ ] Updating xpu cannot silently change the selected allocation backend or
      address-space contract.

---

## Phase 1 — Correctness of the trial wavefunction and local energy

### PHYS-001 — Fix the CUDA incremental Padé gradient

**Priority:** P0

**Verified problem:** For

\[
u(r)=\frac{ar}{1+br},\qquad
u'(r)=\frac{a}{(1+br)^2},
\]

`src/jastrow_pade/jastrow_pade.cu` computes `first_deriv` with one factor of
`a_local` and then multiplies the gradient factor by `a_local` a second time.
The CUDA gradient is therefore proportional to \(a^2\), not \(a\). At the
default \(a=1/4\), its magnitude is one quarter of the correct value. The same
kernel tests `dist_old` before dividing by `dist_new`.

**Why this matters:** The gradient enters
\(|\nabla\log\Psi|^2\) in the kinetic energy. A wrong gradient changes the local
energy and therefore the reported variational energy, even when Metropolis
acceptance uses the correct Jastrow value.

**Proposed fix:**

1. Define one shared, pure pair-derivative helper returning `u_prime_over_r` and
   `u_double_prime + 2*u_prime/r`.
2. Compile that helper for CPU and device rather than maintaining two algebraic
   copies.
3. Use `first_deriv * inv_dist`; do not multiply by `a_local` again.
4. Guard the new contribution using `dist_new`, and define an explicit policy for
   coincident particles as described in PHYS-007.

**Done when:**

- [ ] CPU and CUDA pair derivatives agree with the analytic formulas over a
      parameter grid including multiple non-unit values of `a`.
- [ ] Full and incremental CUDA gradients/Laplacians agree after accepted moves.
- [ ] Independent finite differences of `log(Psi)` agree with CPU and CUDA
      derivatives within precision-specific tolerances.
- [ ] A regression test fails if the extra factor of `a` is reintroduced.

### PHYS-002 — Add periodic Slater determinant and inverse rebuilding

**Priority:** P1

**Measured robustness risk; production impact not established:** Production uses
repeated Sherman–Morrison updates without ever rebuilding the determinant and
inverse. Extending the existing deterministic drift test produced:

| Precision | First observed failure | Maintained inverse residual | Fresh residual |
|---|---:|---:|---:|
| FP32 | 21st accepted update | `9.32617e-2` | `1.46484e-3` |
| FP64 | 331st accepted update | `1.69894e-9` | `2.18279e-11` |

The FP32 predictive determinant-ratio difference was `0.0136337`; the FP64
absolute difference was `1.06214e-7`. The stock test passes because it exercises
only 16 FP32 updates and 64 FP64 updates.

That diagnostic begins from `set_stable_closed_shell_positions`, which places
all seven particles on the line
\(\mathbf r_i=\mathbf r_0+i(1.1,0.7,1.3)\). Reconstructing the actual orbital
matrix gives \(\kappa_2(D)=5626.50\). In a deterministic sample of 20,000
uniform random seven-particle configurations at the same box size, generated
with SplitMix64 seed `20260808`, the median condition number was `19.8876` and
42 of 20,000 (`0.21%`) were at least as ill-conditioned. This uniform-position
comparator is not a sample from the production Markov distribution; it
establishes that the helper is adversarial, not that a typical production
trajectory fails after 21 or 331 updates.

Conditioning does not explain the entire result: at the same failing
configurations, the maintained residual is about `64x` the fresh residual in
FP32 and `78x` in FP64. This isolates update-history error amplified by the
difficult matrix. The audit has not established either a production-energy bias
or the absence of one over arbitrarily long trajectories.

**Why this matters:** Acceptance probabilities, nodes, and derivatives depend on
the maintained determinant and inverse. Accumulated update error can change the
Markov chain and local energy, particularly near nodes. QMCPACK explicitly
documents periodic determinant recomputation because inverse errors accumulate.
The current diagnostic does not justify calling shipped energies untrustworthy,
but a credible long-running implementation still needs a bounded recovery path.

**Proposed fix:**

1. Add a configurable `determinant_rebuild_interval`, expressed in completed
   sweeps or accepted moves and given a conservative precision-dependent default.
2. Recompute the determinant, inverse, and dependent caches from current
   positions at that interval.
3. In debug/validation builds, compute a residual such as
   \(||D D^{-1}-I||_\infty\) and rebuild early when it exceeds a documented
   tolerance.
4. Keep a fixed upper bound on time between rebuilds even if the residual appears
   small; a diagnostic check is not a proof for untested matrices.
5. Implement the same logical cadence in both backends, though the CPU and CUDA
   factorization routines may differ.

**Done when:**

- [ ] Paired rebuild-enabled and rebuild-disabled controls exercise at least
      5,000 accepted updates in FP32 and FP64; the test proves that a rebuild was
      actually triggered and measures its effect rather than passing trivially.
- [ ] Validation includes both production-like Metropolis trajectories and
      deliberately ill-conditioned paths, with maintained-versus-fresh errors
      reported at identical configurations and grouped by condition number.
- [ ] Long-run tests report maximum and quantile residuals, determinant-ratio
      differences, rebuild counts, and the seed/configuration needed to reproduce
      any failure; absence of a trend in one trajectory is not treated as proof.
- [ ] Tests cover a near-singular proposed row, a rejected move, and a forced
      rebuild.
- [ ] Rebuild cadence and residuals are emitted in run metadata.
- [ ] CUDA rebuilding leaves matrices on device; no matrix is copied to the host
      for factorization.

### PHYS-003 — Make accepted moves transactional

**Priority:** P0

**Verified problem:** `Simulation::metropolis_step` updates positions and
`log_psi_current_`, calls `SlaterPlaneWave::accept_move`, and then updates Ewald
caches. `accept_move` silently returns when `1/ratio` is non-finite. If this path
is reached, the logical move can be accepted while determinant and other caches
describe different configurations.

**Why this matters:** The physical walker state is the complete tuple of
positions, wavefunction caches, and energy caches. Partial commits make later
energies neither those of the old nor the new state. RNG state and attempt
counters intentionally advance even for a rejected move and should be handled
separately from the physical-state transaction.

**Proposed fix:**

1. Divide a move into `prepare`, `decide`, and `commit` stages.
2. During `prepare`, compute every delta and validate the determinant ratio and
   proposed numerical state without modifying committed data.
3. During `commit`, update all components cooperatively. If an update cannot be
   represented safely, rebuild from the proposed positions before declaring
   success, or reject without changing any committed component.
4. Return an explicit status (`accepted`, `rejected`, `rebuilt`, or `numerical_failure`)
   instead of silently returning from a component.
5. Use `log_psi_current_` as a periodic consistency check against a full
   recomputation, or remove it if it is not part of the state contract.

**Done when:**

- [ ] Fault-injection tests force each component update to fail and prove that
      the physical state remains old or becomes fully new, never mixed, while
      RNG/counter advancement follows the documented attempt semantics.
- [ ] Full recomputation after randomized accepted/rejected sequences matches all
      incremental caches.
- [ ] Numerical failures are counted and reported, not hidden.

### PHYS-004 — Correct and separate the Perdew–Zunger reference test

**Priority:** P0 for validation

**Verified problem:** `tests/test_known_energy.cpp` combines polarized
`gamma=-0.0843` with unpolarized `beta1=1.0529` and `beta2=0.3334`, then divides
the correlation energy by two. Perdew and Zunger state atomic units and give the
fully polarized large-\(r_s\) parameters `beta1=1.3981` and `beta2=0.2611`.

PZ81 defines its equations with \(\hbar=m=e^2=1\), so these are Hartree atomic
units, not Rydbergs. For example, the polarized fit gives

\[
e_c(1)=\frac{-0.0843}{1+1.3981+0.2611}=-0.031701263538\ \mathrm{Ha}.
\]

Halving this correlation term would produce `1.160896961267 Ha/electron` for
the total at \(r_s=1\); that value is a unit-conversion error. PZ81 Table XIII's
polarized correlation value near `-0.86 eV` at \(r_s=1\) independently agrees
with approximately `-0.0316 Ha`, not half that value.

The correct thermodynamic-limit polarized totals from exact Hartree–Fock kinetic
and exchange plus the PZ81 correlation fit are:

| \(r_s\) | Correct PZ81 total, Ha/electron |
|---:|---:|
| 1 | `1.145046329498` |
| 2 | `0.125784112432` |
| 5 | `-0.060810301537` |
| 10 | `-0.050680495037` |
| 20 | `-0.031235395165` |

**Why this matters:** A passing test against an incorrect published oracle gives
false confidence. Also, a finite-\(N\), Gamma-point VMC calculation is not the
same Hamiltonian and boundary-condition limit as the PZ81 thermodynamic fit, so
their difference is not a direct implementation-error measure.

**Proposed fix:**

1. Add a small, deterministic test of the PZ81 parametrization itself using the
   published polarized coefficients and the values above.
2. Rename the existing end-to-end comparison as a qualitative finite-size sanity
   check unless controlled finite-size and twist corrections are supplied.
3. For quantitative validation, compare against either:
   - a reference code using the same finite cell, twist, Ewald convention, and
     trial wavefunction; or
   - a documented twist-averaged and finite-size-corrected thermodynamic study.
4. Store all reference-code inputs, versions, raw results, and unit conversions
   alongside golden data.

**Done when:**

- [ ] The PZ81 unit test reproduces all five table values above.
- [ ] The unit test independently checks the correlation term before it is added
      to kinetic and exchange energy, and records that the PZ81 fit coefficients
      are in Hartree atomic units.
- [ ] No test labels a finite-cell approximation as an exact thermodynamic
      reference.
- [ ] At least one end-to-end HEG comparison uses demonstrably matching physical
      conventions.

### PHYS-005 — Make the periodic Jastrow differentiable

**Priority:** P0

**Verified mathematical and measured estimator defect:** The current Padé
function is evaluated at a minimum-image distance and has nonzero radial
derivative at the boundary of the Wigner–Seitz cell. The distance is continuous
there, but its derivative changes branch, so the wavefunction gradient is
discontinuous. The distributional Laplacian therefore contains a surface delta
term. Evaluating the ordinary pointwise Laplacian, as the current local-energy
path does, omits that term; the fact that the surface itself has zero sampling
measure does not make its integrated kinetic contribution vanish. CASINO's
documented periodic Jastrow construction instead uses a cutoff at the inscribed
Wigner–Seitz radius and a smoothness choice that keeps the local energy
continuous.

An independent clean-reference FP64 diagnostic compared the gradient-form
kinetic estimator with the implemented regular Laplacian estimator at
\(r_s=2\), \(a=0.25\), and \(b=0.5\). It used 5,000 warmup sweeps, full
derivative evaluation once per measured sweep, determinant rebuilding every ten
sweeps, 200-sweep analysis blocks, and seeds `8675309 + N`:

| \(N\) | Measured sweeps | \((T_{\mathrm{grad}}-T_{\mathrm{regular}})/N\), Ha/electron | Block SE |
|---:|---:|---:|---:|
| 7 | 200,000 | `0.06869` | `0.00524` |
| 19 | 150,000 | `0.10205` | `0.00565` |
| 57 | 50,000 | `0.11944` | `0.00709` |

These values diagnose a material omitted contribution; they are not universal
constants. A quoted correction without \(N\), seed, sampling protocol, rebuild
policy, and uncertainty must not be used as an acceptance value.

**Why this matters:** The local kinetic energy contains first and second
derivatives. Omitting the surface term can lower the reported energy and destroy
the variational interpretation of the value printed by the clean reference CPU
binary. The short-range Padé form also lacks the HEG's conventional long-range
reciprocal-space/RPA correlation behavior. The optimizer's special rejection of
an anomalously low `b_min` result in STAT-003 is consistent with this defect and
must not be used as a substitute for correcting the ansatz.

**Proposed fix:**

1. Select and document a cutoff form satisfying the parallel-spin cusp at the
   origin while bringing the value and the derivatives needed by the local
   energy smoothly to zero at an inscribed-cell cutoff.
2. Derive its value, gradient, Laplacian, and parameter derivatives on paper
   before implementation.
3. Implement the derived value, gradient, Laplacian, incremental-update, and
   parameter-derivative formulas from the same reviewed scalar definition.
4. Add a separate periodic long-range reciprocal-space/RPA Jastrow component if
   the target accuracy requires it; do not claim that Ewald potential evaluation
   supplies wavefunction correlation.
5. Optimize and validate the new ansatz only after its fixed-parameter energy is
   correct.

**Done when:**

- [ ] Value, gradient, and Laplacian tests cover both sides of the cutoff and
      periodic cell faces.
- [ ] The parallel-spin cusp remains \(1/4\) in Hartree atomic units.
- [ ] The local energy remains finite and continuous across the chosen cutoff,
      away from determinant nodes and true Coulomb coalescences.
- [ ] For the corrected smooth ansatz, gradient-form and Laplacian-form kinetic
      estimators agree within a predeclared combined uncertainty on several
      \(N,r_s,b\) cases.
- [ ] The old ansatz's energy impact is measured with a committed, reproducible
      diagnostic that records \(N\), seed, sampling length, rebuild policy, and
      uncertainty; no single correction is generalized to other cases.

### PHYS-006 — Turn Ewald accuracy into a controlled parameter

**Priority:** P1

**Verified result:** The implemented coefficients are algebraically correct for
`alpha=6/L`: self term \(-\alpha N/\sqrt\pi\), neutralizing-background term
\(-\pi N^2/(2\alpha^2V)\), and reciprocal half-space weighting. The current real
sum uses only the minimum image and the reciprocal selection uses a fixed
`1e-6` exponential tolerance.

An independent calculation using additional real lattice images and a tighter
reciprocal cutoff observed approximately:

| \(N\) | \(r_s\) | Difference in total energy | Difference per electron |
|---:|---:|---:|---:|
| 3 | 2 | `4.36e-7` Ha RMS | `1.45e-7` Ha RMS |
| 19 | 1 | `2.10e-5` Ha RMS | `1.10e-6` Ha RMS |
| 19 | 5 | `4.19e-6` Ha RMS | `2.21e-7` Ha RMS |
| 57 | 5 | `2.75e-5` Ha RMS | `4.82e-7` Ha RMS |
| 485 | 5 | about `9.5e-4` Ha observed | about `1.97e-6` Ha |

The existing test named an exact Ewald reference repeats the same minimum-image
real term, alpha, half-space convention, and reciprocal cutoff as production.
It tests internal agreement but is not an independent convergence oracle.

**Why this matters:** Whether microhartree-per-electron error is acceptable
depends on the scientific target. A hard-coded cutoff cannot substantiate an
accuracy claim, and a self-referential test cannot measure truncation error.

**Proposed fix:**

1. Expose an Ewald accuracy target or explicit real/reciprocal cutoffs.
2. Derive real and reciprocal bounds from that target, or converge them
   independently during initialization.
3. Retain the current fast minimum-image mode only if its measured error satisfies
   a named tolerance; label it as such.
4. Add a slow reference implementation that includes explicit lattice images
   and a much tighter reciprocal cutoff and does not call production helpers.
5. Validate full initialization and many incremental accepted-move updates
   against the reference.

**Done when:**

- [ ] Tightening the requested tolerance produces convergent energy.
- [ ] Tests cover several cell sizes, densities, random configurations, and
      accepted/rejected update sequences.
- [ ] Output records alpha, cutoffs, number of G vectors, and requested accuracy.

### PHYS-007 — Handle coalescence and singular limits explicitly

**Priority:** P1

**Verified problem:** Several paths replace `r < 1e-12` behavior with a skip or
an arbitrary finite inverse distance. In the real Ewald update, for example,
`1/r` becomes `1` below the threshold. This is not the Coulomb limit and hides
duplicate-particle states.

**Why this matters:** The same-spin cusp and antisymmetric determinant govern the
physical coalescence limit. Arbitrarily regularizing only selected terms breaks
cusp cancellation and can conceal corrupted positions or invalid proposals.

**Proposed fix:**

1. Distinguish the moved particle's exact self-index, which should be excluded by
   index before evaluating its distance, from two distinct particles that happen
   to coincide.
2. Reject non-finite or physically singular distinct-particle proposals before
   committing them.
3. Where a removable analytic limit exists, implement the derived limit
   consistently in value, gradient, and Laplacian paths.
4. Count and report singular-proposal and rebuild events in debug/validation
   output.

**Done when:**

- [ ] Tests distinguish self terms, near-coalescence pairs, exact duplicate
      positions, and ordinary small distances.
- [ ] No production Coulomb path substitutes an arbitrary finite value for a
      true `1/r` singularity.

### PHYS-008 — Preserve already verified formulas through independent tests

**Priority:** P1

**Verified correct in the CPU reference path:**

- the real plane-wave basis consisting of the constant orbital and cosine/sine
  representatives of nonzero \(\pm\mathbf{k}\) pairs;
- the log-space acceptance calculation for an exactly symmetric single-electron
  proposal, \(\min(1,|D'/D|^2e^{2\Delta J})\), subject to correcting the actual
  finite RNG support in GPU-005;
- the parallel-spin Padé cusp parameter \(a=1/4\);
- CPU Padé value, gradient, and Laplacian formulas;
- the local kinetic-energy identity used by `EnergyTracker`;
- short-horizon determinant ratios and Sherman–Morrison updates; and
- Ewald self/background/reciprocal coefficients described in PHYS-006.

**Why this matters:** Backend redesign should not replace known-correct formulas
without a reference. Preserving behavior by duplicating code is insufficient;
the CUDA Jastrow bug shows that copies can diverge while CPU tests continue to
pass.

**Proposed fix:**

1. Move small scalar formulas into backend-neutral, side-effect-free helpers
   compilable by both C++ and nvcc.
2. Keep independent tests derived from equations, high-precision calculations,
   or external-code golden data rather than reusing those helpers in the oracle.
3. Test complete fixed configurations, not only isolated pair functions.

**Done when:**

- [ ] Every preserved formula has an independent analytic or numerical test.
- [ ] CPU and CUDA invoke the same scalar formula where practical.
- [ ] The reference test implementation shares no reduction, cutoff-selection,
      or cache-update code with production.

---

## Phase 2 — Sampling, statistics, and optimization

### STAT-001 — Replace fixed blocking with a demonstrated correlation analysis

**Priority:** P0 for reported uncertainties

**Measured problem:** `BlockingAnalysis` correctly calculates the variance of
its completed block means, but it never establishes that blocks are independent.
For three independent \(N=19,r_s=5\) runs with 200,000 measured single-particle
moves, block size 100 produced SE/electron values `0.00020956`, `0.00022279`, and
`0.00021531`. Reblocking plateaus over block sizes 500–2000 were approximately
`0.0002400–0.0002429`, `0.0002691–0.0002708`, and
`0.0002595–0.0002664`. The short blocks underestimated SE by roughly 13–20%.

The default `N=485, Block_Size=100` configuration has only 0.206 completed
sweeps per block. Its underestimate was not quantified by this audit and must
not be inferred from the \(N=19\) percentages.

**Why this matters:** An energy without a defensible uncertainty cannot support
comparisons, optimizer decisions, or publication claims. Correlated Markov-chain
samples do not obey the independent-sample SE formula.

**Proposed fix:**

1. Define one Monte Carlo sweep as \(N\) attempted particle moves and make
   measurement cadence explicit in sweeps or substeps.
2. Implement dyadic reblocking or an integrated-autocorrelation estimator that
   exposes all levels and identifies a stable plateau with sufficient blocks.
3. If no stable estimate exists, return `insufficient_statistics`; never report
   a smaller provisional SE as final.
4. On CUDA, accumulate base blocks and reblocking levels on device. Returning the
   compact final level statistics does not violate device residency.
5. Preserve block-level results in output so an independent analysis can be run.
6. Define how an incomplete final block is handled. The current implementation
   silently omits its samples once two complete blocks exist; either include it
   with a statistically correct unequal-block treatment or explicitly discard
   it from both the reported mean and sample count.

**Done when:**

- [ ] Synthetic IID and correlated AR(1) series recover their known uncertainty.
- [ ] CPU and CUDA statistics agree for identical supplied energy sequences.
- [ ] Representative VMC runs demonstrate a plateau with a documented minimum
      number of blocks.
- [ ] Output distinguishes raw sample count, attempted moves, sweeps, base block
      size, chosen blocking level, and effective sample size or autocorrelation.

### STAT-002 — Correct multi-walker aggregation

**Priority:** P0 for reported uncertainties

**Verified problem:** `src/main.cu` sums available per-thread variances but always
divides the resulting SE by the total thread count. If only some walkers provide
an SE, this underestimates uncertainty. Equal averaging is only appropriate when
walkers contribute equal valid sample weights.

**Why this matters:** The requested GPU design will increase the number of
walkers substantially. Aggregation rules must remain correct for failed walkers,
unequal run lengths, and different effective sample counts.

**Proposed fix:**

1. For equal-target chains, define the mean from pooled sample/block sums with
   explicit measurement-count weights. Combine its uncertainty using each
   chain's autocorrelation-aware variance. Use inverse-variance weighting only
   if it is selected and justified as a separate documented estimator; do not
   casually substitute “effective sample count” for a derivation.
2. Require a valid uncertainty contribution from every included walker, or
   exclude and report that walker from both the mean and uncertainty.
3. Report within-walker and between-walker consistency statistics. Treat a large
   between-walker discrepancy as a convergence failure rather than averaging it
   away.
4. Give each walker a unique, recorded RNG stream.

**Done when:**

- [ ] Unit tests cover equal and unequal lengths, missing SEs, failed walkers,
      and analytically known independent Gaussian inputs.
- [ ] The number of contributing walkers and their weights appear in output.

### STAT-003 — Rebuild the Jastrow optimizer as a scientific workflow

**Priority:** P0 because it currently runs by default

**Verified problem:** `JastrowOptimizer` documents a two-phase grid plus
golden-section procedure and a variance-penalized objective. The implementation
performs only a short grid scan, has two grid points when `num_threads=1`, uses
50–500 sweeps per point, applies an ad hoc boundary rejection, and returns
`standard_error=0`. It sets blocking length from `measure_sweeps`, even though
`BlockingAnalysis` counts individual measured moves; at `N=485` the resulting
minimum block of 50 moves is only about 0.103 sweeps. Its comment that Ewald
should handle long-range Jastrow correlation confuses the Hamiltonian with the
trial wavefunction.

The default `config.cfg` has `N=485`, `L=190`, and eight CPU workers, giving
\(r_s=15.0018\), `b_min=0.0666587`, and 16 scan points. The next point after
`b_min` is about `0.39555`, so the special boundary branch compares the
longest-range candidate against a very coarse interior grid. PHYS-005 provides
a concrete mechanism for making that boundary candidate's regular local-energy
estimate artificially low. Whether the branch fires for a particular run must
still be demonstrated with its revision, complete configuration, seed, raw
candidate results, and uncertainties; the source relationship alone is not a
runtime observation.

**Why this matters:** Main always replaces the configured `b` with this result.
Selecting the minimum of noisy short estimates produces selection bias and can
make a run less reproducible or variationally worse. Correcting PHYS-005 removes
one source of false low energies, but does not fix this workflow's independent
sampling, uncertainty, selection-bias, and documentation defects.

**Proposed fix:**

1. Disable automatic optimization by default immediately. A configured,
   recorded `b` must be used unless optimization is requested explicitly.
2. First implement a statistically adequate one-parameter workflow:
   bracket candidates, use sufficiently equilibrated samples, propagate valid
   uncertainties, and confirm the winner with an independent production-length
   run.
3. Prefer correlated sampling on a fixed configuration population when the
   weight stability is demonstrably adequate.
4. If optimization expands beyond one parameter, implement a documented linear
   method or stochastic-reconfiguration method rather than generalizing the
   current noisy grid heuristic.
5. On the GPU, candidate groups may use separate walkers while samples and
   reductions remain on device. Small candidate summaries may be returned only
   at optimization phase boundaries.

**Done when:**

- [ ] Header documentation exactly describes the implemented algorithm.
- [ ] The optimizer never returns an invented zero uncertainty.
- [ ] Repeated independent optimization runs choose statistically consistent
      parameters and the confirmation run verifies the claimed improvement.
- [ ] Fixed-`b` runs bypass all optimizer work and are the default.

### STAT-004 — Validate measurement cadence and warmup policy

**Priority:** P1

**Verified status:** Adapting proposal size only during warmup is valid, and the
current target acceptance of about 50% is a conventional efficiency choice. The
current implementation evaluates a full local energy after every one-particle
attempt. This remains an unbiased stationary estimator but creates strongly
correlated, expensive samples.

**Why this matters:** On GPU, full-energy evaluation after every attempted move
can dominate runtime and conceal the benefit of parallel walkers. Statistical
meaning must not depend on ambiguous uses of “step” and “sweep.”

**Proposed fix:**

1. Define `attempt`, `sweep`, `substep`, `measurement`, and `block` in code and
   documentation.
2. Permit multiple full sweeps/substeps between measurements.
3. Keep adaptation strictly in warmup and freeze all transition parameters
   before measurements begin.
4. Bound the adapted proposal scale away from numerical underflow and at a
   documented maximum such as the existing `L/2`; report when a bound is reached.
5. Validate that chosen cadence reduces cost without changing the mean outside
   combined confidence intervals.

**Done when:**

- [ ] Output reports all five counts unambiguously.
- [ ] Tests prove no adaptation occurs during measurement.
- [ ] Cadence studies report energy, SE, autocorrelation, and cost per effective
      sample.

### STAT-005 — Make proposal and particle-selection distributions exact

**Priority:** P1

**Verified problem:** The Metropolis ratio in PHYS-008 is correct for an exactly
symmetric proposal, but the finite RNG mappings do not establish that condition.
The CPU uses `std::uniform_real_distribution(-step_size, step_size)`, whose
specified interval is half-open. CUDA maps `curand_uniform[_double]`, whose
documented support is `(0,1]`, to `(-step_size,+step_size]`. The CUDA mapping
therefore has an unmatched positive endpoint. The C++ half-open contract permits
the opposite unmatched endpoint but does not specify a universal raw-bit-to-real
mapping, so it fails to prove exact sign symmetry without an implementation-level
support/mass check. CUDA also chooses a particle with `curand() % N`, which has
modulo bias whenever `N` does not divide the generator's integer range.

The magnitude is bounded well enough to set priority honestly. The documented
FP64 cuRAND mapping uses 53 random bits, so one unmatched grid endpoint is an
order-\(2^{-53}\) event; the FP32 mapping is materially coarser and must have its
exact mass verified against the selected implementation. For a uniform 32-bit
integer mapped with `% 7`, four particles receive one extra raw integer: the
absolute heavy-versus-light probability difference is exactly \(2^{-32}\), and
the relative imbalance is about `1.63e-9`. The C++ standard specifies the CPU
interval but not one universal raw-bit mapping, so its endpoint mass must be
measured for the supported standard-library implementation rather than assumed.

**Why this matters:** The acceptance formula omits a Hastings proposal-ratio term
and therefore requires `q(delta)=q(-delta)`. The endpoint probability is small,
but an exact correctness claim should not depend on never drawing it. A fixed,
state-independent nonuniform mixture of individually valid particle-move kernels
can still preserve the target distribution, so the modulo issue is not by itself
proof of an energy bias; it is nevertheless unnecessary, contradicts the
intended random-scan distribution, and gives particles unequal update rates.
These defects should be removed for an exact and reviewable kernel, but the
verified magnitudes do not justify classifying current energies as untrustworthy
ahead of the much larger PHYS-005 defect.

**Proposed fix:**

1. Construct displacement samples from raw random bits with a support and mass
   function exactly invariant under `delta -> -delta` (for example, a paired
   midpoint grid), then scale by `step_size`. Merely swapping which side of an
   interval is open is not sufficient.
2. Use a proved rejection method for an unbiased particle integer in `[0,N)`;
   a multiply-high implementation is acceptable only with the rejection step
   required to remove its remainder bias.
3. Accept unconditionally when the log ratio is nonnegative. Otherwise define
   the acceptance uniform's raw-bit mapping and discrete interval convention,
   test its exactly-zero and maximum representable values against the log-space
   comparison, and document the finite-grid acceptance-probability error bound.
4. Specify these distributional contracts once for CPU and CUDA even if their
   engines and trajectories differ.

**Done when:**

- [ ] Exhaustive reduced-bit tests prove proposal sign symmetry and integer
      uniformity, including the endpoint cases.
- [ ] Statistical tests cover displacement moments and particle frequencies for
      awkward values such as `N=7`, `19`, and `485`.
- [ ] A small enumerated transition test using acceptance ratios exactly
      representable on the chosen random grid fails for the old endpoint mapping
      and satisfies detailed balance with the replacement; nonrepresentable
      ratios satisfy the documented discretization bound.

---

## Phase 3 — Configuration, errors, and reproducibility

### CFG-001 — Make invalid scientific inputs impossible to run silently

**Priority:** P0

**Verified problem:** Configuration validates particle count, block size, box
length, and measurement sweeps, but not `num_threads >= 1`, finite Jastrow
parameters, a safe domain for `b`, or overflow in `N * sweeps`. A missing file
silently uses defaults. Unsigned parsing does not reject syntactically negative
text explicitly, conversions do not verify that the complete input string was
consumed, and the generic parser does not check the destination range. Malformed
lines without `=` are skipped, duplicate keys overwrite earlier values, and
unknown keys are only warned about. The Boolean helper would interpret any text
other than `true` or `1` as false.

Particle count alone also does not state an occupation policy. The initializer
fills real orbitals in `(n^2,n_x,n_y,n_z)` order and can stop partway through a
degenerate shell. That is a deterministic Slater determinant, but it is not a
closed-shell Gamma-point HEG state and should not be presented as one implicitly.

**Why this matters:** A typo can launch a valid-looking simulation with different
physics. `num_threads=0` results in no chains and invalid final arithmetic. A
negative `b` can place a pole in `1+b*r`.

**Proposed fix:**

1. Parse into checked intermediate types using complete-string conversion.
2. Require finite `a`, finite `b`, a documented allowed `b` domain, positive
   finite box length and proposal scale, at least one walker/thread, and enough
   measurements for the requested statistical estimator.
   For the current fully polarized model, either enforce the cusp value
   `a=0.25` for production or require an explicit, recorded override for a
   non-cusp trial state.
3. Check multiplication and allocation-size overflow before computing derived
   counts.
4. Treat missing configuration, unknown keys, malformed lines, and invalid
   Booleans as errors by default; reject duplicate keys rather than silently
   taking the last. Provide an explicit `--allow-defaults` mode if desired for
   examples.
5. Make derived values immutable after validation or recompute them through one
   validated constructor rather than writing fields manually in the optimizer.
6. If the supported scientific model is closed-shell Gamma-point HEG, reject
   particle counts that do not fill a complete momentum shell. If open shells
   are intentionally supported, require/document their occupation and emit the
   exact ordered orbital list so the state is reproducible.

**Done when:**

- [ ] Boundary and malformed-input tests cover every field and overflow path.
- [ ] Invalid input exits nonzero before allocating simulation state.
- [ ] The final normalized configuration is written to output.
- [ ] Closed/open-shell tests verify occupation, degeneracy tie-breaking, and
      the model label written to results.

### CFG-002 — Add a real command-line contract

**Priority:** P1

**Verified problem:** `main` ignores command-line arguments, hard-codes
`config.cfg`, always runs optimization, and always targets `output/vmc.bin`.

**Why this matters:** The requested runtime backend selection cannot exist
without a real CLI. Hard-coded paths and mandatory optimization also make it too
easy to run different physics than intended while producing a familiar-looking
output file.

**Proposed fix:**

Support, at minimum:

```text
vmc [--cuda | --backend cpu|cuda]
    --config PATH
    [--output PATH]
    [--seed UINT64]
    [--walkers COUNT]
    [--optimize-jastrow]
    [--checkpoint-interval BLOCKS]
```

No backend fallback should be silent. `--cuda` with a missing CUDA artifact,
missing device, insufficient memory, or incompatible device must fail clearly.

**Done when:**

- [ ] CLI integration tests cover defaults, aliases, conflicts, missing files,
      backend failures, and exit codes.
- [ ] `vmc --help` states exactly which computation and transfers occur in each
      backend.

### CFG-003 — Record complete provenance

**Priority:** P1

**Verified problem:** Current output omits backend, precision, device, source
revision, xpu revision, compiler, optimizer details, Ewald cutoffs, determinant
rebuild cadence, and walker seed mapping.

**Why this matters:** Energies and error bars cannot be reproduced or audited if
the executable, numerical policy, hardware backend, and random streams that
created them are unknown.

**Proposed fix:**

Record a versioned run manifest containing:

- normalized input and derived \(r_s\);
- physical model and units;
- VMC and xpu revisions and dirty-tree indicator;
- CPU/CUDA backend, precision, compiler, build type, and fast-math status;
- GPU name, compute capability, runtime/driver versions, and walker count;
- every walker seed/subsequence;
- ansatz parameters and whether/how they were optimized;
- exact orbital occupations, shell closure, and twist;
- Ewald parameters and convergence target;
- determinant rebuild policy and observed early rebuilds;
- measurement/reblocking configuration and contributing walkers; and
- start/end times and completion/failure status.

**Done when:**

- [ ] A result file alone is sufficient to reconstruct the run configuration.
- [ ] CPU and CUDA emit the same schema and units.

---

## Phase 4 — Separate CPU and CUDA products

### ARCH-001 — Use separate backend artifacts behind one launcher

**Priority:** P0 for the requested architecture

**Verified constraint:** xpu documents that `XPU_CUDA` selects either aligned
host allocation or `cudaMalloc`, must have the same value in every translation
unit linked into a program, and cannot mix CPU and CUDA backends in one linked
program.

**Why this matters:** A single ordinary executable cannot safely instantiate
both current xpu modes and select between them at runtime. Trying to do so would
weaken the type/memory boundary that is supposed to prevent host dereferences of
device pointers.

**Proposed fix:**

```text
vmc                         small CPU-only launcher
├── default  -> vmc-cpu     xpu compiled with XPU_CUDA=OFF
└── --cuda  -> vmc-cuda    xpu compiled with XPU_CUDA=ON
```

1. Build CPU and CUDA implementations in separate CMake sub-builds or a
   superbuild, because the one `xpu::xpu` target cannot carry both definitions in
   one configure/link context.
2. Have the launcher validate arguments and replace itself with the chosen
   executable, preserving standard I/O and exit status.
3. Share source-level request/result schemas and physics helpers, not storage
   objects or linked xpu instantiations.
4. Prefer this process boundary over dynamically loaded libraries. A one-process
   design would require a narrow C ABI, hidden visibility, and proof that inline
   xpu symbols cannot interpose between contradictory definitions.

**Done when:**

- [ ] `vmc` without a backend flag executes `vmc-cpu`.
- [ ] `vmc --cuda` executes `vmc-cuda` and never silently falls back.
- [ ] `vmc-cpu` has no CUDA runtime dependency and runs on a CUDA-free host.
- [ ] `vmc-cuda` fails clearly when no compatible device is available.
- [ ] Both artifacts consume the same validated configuration schema and emit
      the same result schema.

### ARCH-002 — Define explicit backend-neutral boundaries

**Priority:** P1

**Verified problem:** Current domain classes combine model equations, memory
ownership, host orchestration, CUDA launches, output, and cache mutation behind
`#ifdef XPU_CUDA`. This makes address space and execution location implicit.

**Why this matters:** Implicit ownership permits host code to compile around
device pointers and encourages divergent copies of physics formulas. A narrow
boundary is also necessary for the two separately linked backends to consume the
same scientific request and produce comparable results.

**Proposed fix:**

Define narrow value-level interfaces such as:

```cpp
struct RunRequest;        // validated immutable values, no owning raw pointers
struct RunResult;         // estimates, diagnostics, provenance, no backend data
struct ModelParameters;   // N, L, spin/model, Jastrow and Ewald definitions
```

Then maintain backend-specific ownership:

- `CpuSimulation` and `CpuWalkerState` own host xpu buffers and CPU RNGs.
- `CudaSimulation` owns device allocations through host-side RAII handles as
  well as launch metadata, but never dereferences the allocation contents.
- `DeviceWalkerState` is a device-appropriate view/POD passed to kernels and is
  never dereferenced by host code.
- Small scalar physics functions may be host/device callable; orchestration and
  owning classes must not pretend to be portable through one macro.

**Done when:**

- [ ] Backend-neutral headers contain no CUDA runtime types or xpu-owned storage.
- [ ] Address-space ownership is evident from types and naming.
- [ ] No class changes its data-member layout under `__CUDA_ARCH__`.

### ARCH-003 — Complete a pure CPU backend first

**Priority:** P1

**Why this matters:** The CPU backend is the inspectable scientific reference for
the CUDA implementation. It must remain independently usable and must not depend
on nvcc, curand, cudart, or a GPU.

**Proposed fix:**

1. Compile CPU implementation files as C++ rather than treating `.cu` files as
   C++ with `-x c++`.
2. Use host xpu buffers, checked views, and a dedicated CPU RNG type.
3. Preserve parallelism across independent walkers with `std::jthread`, a bounded
   task pool, or OpenMP; do not conflate CPU workers with CUDA blocks.
4. Make the CPU backend pass all reference physics and statistics gates before it
   becomes the CUDA oracle.

**Done when:**

- [ ] CPU configure/build/test succeeds on a machine without a CUDA toolkit.
- [ ] Link inspection confirms no CUDA libraries.
- [ ] CPU production runs pass all release-quality scientific tests in Phase 6.

---

## Phase 5 — Device-resident multi-walker CUDA engine

### GPU-001 — Store all mutable walker state on device

**Priority:** P0 for CUDA correctness

**Verified problem:** xpu CUDA buffers use `cudaMalloc` and cannot be
host-dereferenced, but current VMC code does so in several places: `Simulation`
indexes positions on the host, `EnergyTracker` uses `std::copy_n` into xpu CUDA
storage, Slater initialization reads a device counter and sorts device arrays
with host algorithms, and the Jastrow path reads moved positions on the host.

**Why this matters:** These operations are invalid for ordinary `cudaMalloc`
pointers and are the direct architectural contradiction the xpu migration is
supposed to remove. Even replacing them with a scalar copy at each call would
still violate the requested device-resident execution model.

**Proposed fix:**

1. Build immutable orbital/G-vector tables in checked host staging storage, then
   perform one explicit host-to-device transfer, or generate them with device
   kernels. Do not apply host algorithms to device pointers.
2. Allocate walker-major device storage once. Avoid allocation inside individual
   moves or local-energy evaluations.
3. Initialize positions, RNGs, determinants, Jastrow state, Ewald structure
   factors, and counters on device.
4. Keep only counts, launch configuration, stream/event state, and output staging
   handles on the host.

**Done when:**

- [ ] CUDA API tracing shows no `cudaMallocManaged`.
- [ ] Compute-sanitizer reports no invalid host/device accesses.
- [ ] A device-residency trace shows no per-move H2D/D2H copies.
- [ ] Debug code can query pointer attributes and assert that every walker array
      has the expected address space.

### GPU-002 — Implement one CUDA block per walker as a benchmarkable policy

**Priority:** P1

**Design requirement:** Independent Markov chains provide natural coarse-grain
parallelism. One block per walker is viable, but production QMC codes batch and
tune walkers/crowds; the best mapping cannot be assumed without measurement.

**Why this matters:** A single Markov chain is sequential, so useful GPU
parallelism must come from independent walkers and cooperative work within each
move. A mapping that exhausts registers/shared memory or leaves lanes idle can be
correct yet substantially slower than another batching policy.

**Proposed execution model:**

1. `blockIdx.x` identifies a walker.
2. One lane owns that walker's Philox state, selects a particle, creates a
   proposal, and broadcasts it through shared memory.
3. Threads cooperate on Jastrow deltas, determinant ratios and updates, Ewald
   deltas, derivatives, and energy reductions.
4. One lane performs the accept/reject comparison after cooperative calculations
   complete; the block then commits or discards all prepared state.
5. Each block advances its chain sequentially for a chunk of many attempted moves
   before returning control to the host.
6. Different blocks/walkers never share mutable Markov state.

**Why chunked kernels:** A host launch at a chunk boundary does not move walker
data off device. Chunking permits determinant rebuild kernels, checkpoints, error
inspection, and watchdog-safe execution without a host round trip for every move.

**Done when:**

- [ ] Walker independence is tested by changing grid order/count without changing
      aggregate statistical results beyond expected uncertainty.
- [ ] Benchmarks compare one-block-per-walker with at least one batched/crowd
      alternative for representative small and large `N`.
- [ ] The selected mapping is justified by effective samples per second, not raw
      attempted moves per second.

### GPU-003 — Eliminate per-operation synchronization and scalar copies

**Priority:** P1

**Verified problem:** Current CUDA branches perform `cudaDeviceSynchronize` and
copy scalar values to the host in determinant ratios, Jastrow values/deltas,
kinetic energy, and Ewald updates. That is host-driven accelerator offload, not a
device-resident simulation.

**Why this matters:** Per-operation synchronization serializes work across
walkers, while scalar round trips make the host part of every Markov transition.
Both defeat the requested backend contract and severely limit GPU throughput.

**Proposed fix:**

1. Fuse the sequence needed for one move or a move chunk into walker kernels.
2. Store ratios, deltas, energy components, and decisions in registers, shared
   memory, or walker device state.
3. Use block reductions with a deterministic tree for testability. Avoid global
   floating-point atomics when the result belongs to one walker.
4. Use streams/events only at chunk, rebuild, checkpoint, or final boundaries.
5. Check launch errors asynchronously at safe boundaries rather than forcing a
   device-wide synchronization after every small kernel.

**Done when:**

- [ ] Nsight Systems shows no device-wide synchronization inside the per-move
      steady-state loop.
- [ ] No per-move scalar `cudaMemcpy` remains.
- [ ] Component and end-to-end tests remain correct after fusion.

### GPU-004 — Implement device determinant rebuilds

**Priority:** P1, dependent on PHYS-002's rebuild policy

**Why this matters:** PHYS-002 establishes an adversarial update-history error
and a credible-release requirement for bounded rebuilding; it does not establish
typical production failure at a particular update count. Once rebuilding is
required, copying every walker matrix to the host would restore numerical
accuracy at the cost of abandoning device residency.

**Proposed fix:** First implement one correct, supported baseline with matrices
remaining on device. A conventional cuSOLVER factorization/solve launched at a
walker-chunk boundary is an acceptable candidate when the installed API supports
the required matrix sizes and layout: host launch control does not violate the
residency contract when walker matrices, pivots, status, and device work arrays
remain device-resident. If the selected API requires host workspace, document
its specified role and trace the call: accept it only if no VMC-owned walker
state or numerical result is staged there and the factorization executes on the
GPU. Do not infer that from the pointer name. Select and record the actual API
only after compiling a minimal supported prototype against the pinned CUDA
toolkit.

Treat comparison with cuSolverDx, a custom tiled block/cooperative
factorization, or another batched route as a bounded design/performance spike.
It is not a prerequisite for the first correct backend. Require a custom path
only if the baseline fails a declared capability or performance target while
preserving the same numerical acceptance tests.

The conventional cuSOLVER dense API takes device matrix, pivot, status, and
device-workspace pointers while being launched through a host handle; the generic
`Xgetrf` API also requires a host workspace. Device matrix pointers establish
storage residency but, by themselves, do not establish the host workspace's
role or the absence of hidden transfers; the prototype and trace are the
acceptance evidence. cuSolverDx documents device-side batched LU operators, but
it is a MathDx component and therefore adds a distribution/build/licensing
decision beyond toolkit cuSOLVER. Exact library choice must be confirmed against
the installed CUDA version and benchmarked; this plan does not assume that every
vendor routine supports every desired layout or batch shape.

**Done when:**

- [ ] Rebuilt device determinants/inverses agree with the CPU reference on fixed
      matrices and randomized walker states.
- [ ] Pivot/singularity status is stored and handled per walker.
- [ ] FP32 and FP64 residuals satisfy documented tolerances over long runs.
- [ ] No determinant matrix is copied to host for routine rebuilding.
- [ ] Profiling confirms that the selected vendor path performs no hidden
      walker-state/result transfer and satisfies the documented residency
      interpretation for any host workspace it requires.
- [ ] The first correct baseline has a reproducible capability and timing record;
      any decision to pursue a second factorization implementation names the
      failed target that justifies it.

### GPU-005 — Use separate, unbiased RNG implementations

**Priority:** P1

**Verified problem:** `WalkerRNG` changes member layout under `__CUDA_ARCH__`.
On the audited GCC 15/libstdc++ and CUDA 13.3 FP64 toolchain, a compiled
host/device probe measured `sizeof(WalkerRNG)==2568` in host code and `80` in
device code on the RTX 3070. The exact byte counts are ABI-specific, but the
layout mismatch is inherent in the conditional member list. The distributional
defects in its proposal and particle mappings are addressed in STAT-005. The CPU
path seeds each `std::mt19937_64` with
`master_seed + walker_id`; this gives distinct seed integers but is not a
documented stream-splitting or non-overlap guarantee. CUDA 13.3 documentation
also warns that the legacy cuRAND device API used by `curand_kernel.h` will be
deprecated in a future release. The documented newer device-side extension,
cuRANDDx, is distributed through MathDx rather than as part of the CUDA Toolkit.

**Why this matters:** A type whose host and device layouts differ cannot safely
serve as shared walker storage. Independent walkers also require explicit stream
ownership and a stable dependency boundary for reproducibility and statistical
validation.

**Proposed fix:**

1. Create separate `CpuWalkerRng` and `DeviceWalkerRngState` types behind a
   narrow CUDA RNG adapter; do not expose cuRAND-specific state in shared model
   interfaces. Select legacy cuRAND, cuRANDDx, or another reviewed counter-based
   implementation only after recording its support, build, and licensing cost.
2. Specify stream derivation separately for each backend. For CUDA Philox, key
   subsequences by `(master_seed, walker_id)`, retain checkpointable offsets, and
   bound each walker's consumption below the documented subsequence length. For
   the selected CPU engine, use a reviewed stream-splitting or seed-expansion
   procedure; do not call `master_seed + walker_id` a non-overlap proof.
3. Implement the exact distribution mappings specified by STAT-005 inside each
   backend RNG adapter.
4. Do not require CPU and GPU to produce identical trajectories. Require each to
   be reproducible within its backend and statistically equivalent across
   backends.

**Done when:**

- [ ] The CUDA adapter passes every distributional contract test in STAT-005.
- [ ] Changing walker count does not reuse RNG streams.
- [ ] Cross-stream tests and the documented maximum draws per walker support the
      chosen separation policy; the plan does not infer independence merely from
      distinct integer seeds.
- [ ] Checkpoint/restart reproduces the continuation of every walker.
- [ ] The selected CUDA RNG provider/version and stream mapping are recorded in
      the run manifest, and dependency upgrades have deterministic replay and
      statistical regression tests.

### GPU-006 — Add memory planning and graceful capacity checks

**Priority:** P1

**Verified sizing concern:** At `N=485`, three FP64 `N x N` arrays alone occupy
about 5.4 MiB per walker. Two `N x 243` FP64 trigonometric arrays add about
1.8 MiB. With 128 walkers, those arrays alone approach 0.9 GiB before Ewald state,
scratch, alignment, output, or library workspace.

**Why this matters:** Walker count determines both statistical throughput and
memory use. Discovering an over-capacity request through an aborting allocator
loses diagnostics and can prevent the launcher from recommending a viable run.

**Proposed fix:**

1. Define a function that computes exact bytes per walker and shared immutable
   bytes for the selected precision and algorithm.
2. Include factorization workspace and peak temporary allocations, not only
   persistent arrays.
3. Query available device memory and choose or reject walker count before
   allocation. Never rely on xpu's allocation abort as normal capacity control.
4. Print the estimate and selected occupancy-relevant launch configuration.

**Done when:**

- [ ] Unit tests detect size arithmetic overflow.
- [ ] CUDA integration tests cover a valid allocation and a deliberately
      over-capacity request with a clear nonzero exit.
- [ ] Benchmarks report memory per walker and maximum resident walkers.

### GPU-007 — Design output and checkpointing around residency

**Priority:** P1

**Verified problem:** The current writer flattens and writes all positions after
every measured particle move. Doing that in CUDA mode would require continuous
device-to-host traffic and prevent the requested residency model.

**Why this matters:** Output policy is part of the backend architecture, not a
late formatting detail. A default trajectory stream could silently reintroduce
the transfers that the CUDA redesign is intended to eliminate.

**Proposed fix:**

1. Make block summaries the default output.
2. Make trajectories explicitly optional and sampled at a user-specified block
   interval.
3. For optional snapshots, copy complete checkpoint batches through pinned
   staging storage at chunk boundaries; never read individual device positions
   from host code.
4. Retain block summaries on device until the final result transfer unless the
   user explicitly requests live output or checkpointing. For live output, copy
   summaries in batches only at declared chunk/checkpoint boundaries.
5. A checkpoint must include positions, all state needed for a mathematically
   identical continuation, RNG states, counters, rebuild cadence, and statistics.

**Done when:**

- [ ] A no-trajectory, no-checkpoint CUDA run transfers only validated immutable
      setup data before sampling and one compact result after sampling; walker
      initialization itself remains on device.
- [ ] Live output and checkpoint modes transfer only documented bulk batches at
      chunk/checkpoint boundaries, and report their transfer volume.
- [ ] Restarted CPU and CUDA runs reproduce their respective uninterrupted runs.
- [ ] Snapshot cost and transfer volume are visible in run metadata.

---

## Phase 6 — Verification suite and external validation

### TEST-001 — Build an independent fixed-configuration oracle suite

**Priority:** P0

**Why this matters:** End-to-end stochastic energies can pass broad tolerances
while individual formulas are wrong. Fixed configurations remove sampling noise
and localize discrepancies.

**Required golden quantities:**

- orbital values and `k`/shell assignments;
- determinant, log absolute determinant, inverse residual, and move ratio;
- Jastrow value, move delta, gradient, Laplacian, and cutoff behavior;
- real, reciprocal, self, and background Ewald terms;
- kinetic and total local energy;
- all incremental-cache states after accepted and rejected moves.

**Proposed fix:**

1. Generate small deterministic configurations away from nodes and singularities.
2. Compute golden values with high precision or an implementation that does not
   share production helpers.
3. Store configuration, units, equations, generator version, and expected values
   in a human-readable versioned format.
4. Run the exact same dataset through CPU and CUDA backends.

**Done when:**

- [ ] Every quantity above has precision-specific tolerances justified by
      conditioning and reduction order.
- [ ] Golden-data generation is reproducible and reviewed independently of the
      production implementation.

### TEST-002 — Add CPU/CUDA parity tests at component and state levels

**Priority:** P0

**Why this matters:** The CUDA Jastrow error survived because CPU tests did not
exercise the independent CUDA formula. Parity at only the final noisy energy is
too weak to identify which component or state transition diverged.

**Proposed fix:**

1. Compare fixed inputs, not backend RNG trajectories.
2. Test initialization, one proposal, rejection, acceptance, long move sequences,
   determinant rebuilding, Ewald updates, reblocking, and checkpoint/restart.
3. Require statistical agreement of independent full runs using predeclared
   confidence criteria; do not require bitwise equality between different RNGs
   or floating-point reduction orders.

**Done when:**

- [ ] CUDA CI executes these tests on a real GPU.
- [ ] A deliberately perturbed CUDA formula is detected by the suite.
- [ ] Full-run comparison criteria are documented before looking at results.

### TEST-003 — Cross-check with established QMC codes

**Priority:** P1

**Why this matters:** CPU/CUDA agreement can prove that two local implementations
match while both share the same conceptual error. An established independent code
provides a separate implementation lineage, provided the Hamiltonian, boundary
conditions, and trial state genuinely match.

**Proposed fix:**

1. Start with a bounded feasibility spike and decision record. For QMCPACK and
   CASINO, enumerate which cell, spin, twist, Ewald, orbital, and Jastrow
   conventions can actually be made identical; identify the smallest exposed
   fixed-configuration quantity; and record the code version plus build/run
   requirements. Select one code only on that evidence.
2. If no identical case is representable or the required components are not
   exposed, stop the spike, record the specific mismatch, and do not force a
   total-energy comparison into service as an oracle.
3. For a feasible case, first compare fixed-configuration wavefunction ratios
   and local-energy components. Treat this as the deliverable independent of any
   later stochastic campaign.
4. Only after all physical conventions and fixed-configuration quantities agree,
   run a separately scoped statistically converged VMC comparison with
   predeclared uncertainty and agreement criteria.
5. Store external input decks, code release/commit, execution commands, raw
   outputs, and conversion scripts in `tests/reference/` or an equivalent
   provenance directory.

**Done when:**

- [ ] The feasibility decision is reproducible and states every relevant matched
      convention, exposed quantity, dependency, and known mismatch.
- [ ] If feasible, at least one fixed-configuration comparison is reproducible
      by another developer; otherwise the plan records why the external code
      cannot be an exact oracle for this ansatz.
- [ ] A controlled stochastic comparison is a separate follow-up and is accepted
      only when finite-size and twist conventions are exactly matched or their
      effect is independently quantified rather than assumed away.

### TEST-004 — Add property, stress, and failure tests

**Priority:** P1

**Why this matters:** Example-based unit tests cover only selected configurations.
Periodic symmetry, long update sequences, singular cases, and failure paths span
a much larger state space and are where several audited defects appear.

**Proposed fix:** Add reproducible property and stress tests for the following
behaviors, with precision-specific tolerances and failure seeds/configurations
printed in the test report.

Required properties include:

- periodic translation invariance;
- particle-permutation behavior of determinant sign and `|Psi|`;
- detailed-balance acceptance formula for symmetric proposals;
- free-gas zero variance for supported shells;
- energy/cache agreement after long randomized move sequences;
- singular and near-singular determinant behavior;
- invalid configurations and allocation overflow;
- zero/one/many walkers and unequal walker completion;
- FP32/FP64 determinant drift and Ewald convergence; and
- deterministic replay within each backend.

Property tests must use reproducible seeds and print the seed/configuration on
failure.

**Done when:**

- [ ] Long-running stress tests are separated from fast unit tests but run on a
      scheduled CI cadence.
- [ ] Every previously confirmed defect in this document has a regression test.

### TEST-005 — Verify CUDA residency and safety with tooling

**Priority:** P0 for the CUDA release

**Why this matters:** Ordinary numerical tests do not reliably detect invalid
address-space accesses, uninitialized device data, shared-memory races, divergent
barriers, or accidental host/device transfers. These are defining risks of the
new backend.

**Proposed fix:** Add self-hosted or otherwise GPU-capable CI stages for:

- CUDA runtime tests;
- compute-sanitizer memcheck;
- racecheck where supported and applicable;
- initcheck/synccheck where supported;
- an Nsight Systems residency/transfer smoke trace; and
- an automated source check rejecting `cudaMallocManaged`.

**Done when:**

- [ ] The steady-state move loop has no host/device transfer or device-wide
      synchronization according to the trace.
- [ ] All enabled compute-sanitizer tools report zero errors.
- [ ] Tool versions and GPU model are recorded in CI artifacts.

---

## Phase 7 — Code quality, CI, documentation, and release polish

### QUAL-001 — Remove CPU/CUDA formula duplication and clarify ownership

**Priority:** P1

**Verified problem:** Large functions interleave `#ifdef XPU_CUDA` branches with
separate implementations of the same physics. Raw pointers do not communicate
host versus device address space. This contributed directly to the CUDA Jastrow
error.

**Why this matters:** Scientific formulas must have one reviewable definition
where practical, while ownership and scheduling must remain explicit. The current
mixture makes both goals harder and lets one backend regress unnoticed.

**Proposed fix:**

1. Put the project in a `vmc` namespace.
2. Share small pure formulas, but keep scheduling, storage, RNG, and error
   handling in backend-specific code.
3. Use named component enums/views instead of anonymous positional pointers at
   call sites.
4. Keep ownership RAII-based and make non-owning views explicit.
5. Restrict `noexcept` to functions whose complete failure behavior satisfies it;
   do not hide numerical or CUDA failures behind silent returns.

**Done when:**

- [ ] Physics source contains no large CPU/CUDA preprocessor forks.
- [ ] Code review can identify address space and ownership from types.
- [ ] Public headers expose the minimal stable API needed by callers/tests.

### QUAL-002 — Make output formats durable and complete

**Priority:** P1

**Verified problem:** CSV methods throw intentionally. The binary output has no
magic number, schema version, endian/precision tag, checksums, or complete
metadata, and `output/` is not created automatically. `render.py` always decodes
a little-endian FP64 header (`<QdQ`) and FP64 frame values, while the C++ writer
serializes native `real_t` values. Consequently, an FP32 result has a different
header/frame width and cannot be decoded correctly by the current renderer.
Only the first CPU worker receives a writer, `BinOutputWriter::write_done`
discards every field of `DoneData`, and the cross-worker aggregate is printed
only to standard output rather than stored in the result file.

**Why this matters:** An unversioned native binary cannot be interpreted safely
after precision, platform, or record layout changes. Silent or late output
failure can also waste a long simulation or leave a file that looks complete.

**Proposed fix:**

1. Choose a versioned portable result/checkpoint format. If the raw binary format
   remains, specify magic, version, byte order, scalar representation, record
   lengths, and failure detection.
2. Implement CSV for block summaries or remove it from the public API until it is
   supported.
3. Create parent directories safely, use temporary files plus atomic rename for
   completed summaries/checkpoints, and report I/O failures. Preflight a required
   output path before simulation rather than warning and continuing without the
   requested artifact.
4. Do not encode “missing SE” as numeric zero.
5. Store per-walker diagnostics and the final multi-walker aggregate, including
   its estimator, weights, contributing walkers, and uncertainty status.

**Done when:**

- [ ] Round-trip and truncated/corrupt-file tests exist.
- [ ] FP32 and FP64 files are both decoded according to an explicit schema and
      covered by renderer/reader integration tests.
- [ ] Readers reject unsupported versions clearly.
- [ ] Output includes the provenance manifest from CFG-003.
- [ ] The durable result contains the same final aggregate reported by the CLI;
      a run does not depend on captured terminal text for its primary result.

### QUAL-003 — Build a meaningful CI matrix

**Priority:** P1

**Why this matters:** The current single Ubuntu Release job cannot protect FP32,
CPU-only, CUDA, sanitizer, or documentation behavior. A backend may therefore be
broken or scientifically divergent while the only required check remains green.

**Proposed fix:** Replace the starter workflow with the following minimum matrix,
using a real/self-hosted GPU for runtime and sanitizer jobs where hosted runners
cannot provide one.

The minimum CI coverage should include:

| Job | Required coverage |
|---|---|
| CPU Debug FP64 | complete unit/integration suite |
| CPU Debug FP32 | complete suite and drift thresholds |
| CPU Release FP64 | production flags and smoke run |
| CPU ASan+UBSan | complete compatible suite |
| CPU LeakSanitizer | separate compatible run |
| CUDA compile | supported nvcc/host compiler and architectures |
| CUDA real GPU | parity, runtime, residency, and sanitizer tests |
| Formatting | `clang-format --dry-run --Werror` |
| Static analysis | selected `clang-tidy` checks and documented suppressions |
| Documentation | links, examples, and schema/reference-data checks |

The clean audited reference revision already compiles on CPU with zero project
warnings under `-Wall`, `-Wextra`, `-Wpedantic`, `-Wshadow`, `-Wconversion`, and
`-Wsign-conversion`; enable warnings-as-errors for that CPU baseline immediately
and preserve it while resolving the current migration build failure. Track CUDA
diagnostics separately: the clean CUDA 13.3 build emits 16 warnings from applying
`CUDA_CALLABLE` to `metropolis_step` while its callees remain host-only. Remove
that invalid annotation boundary or make the complete call graph genuinely
device-callable; do not hide these diagnostics with a blanket suppression. Pin
major toolchain versions or publish the exact tested range.

**Done when:**

- [ ] Required checks protect the main branch.
- [ ] CPU strict-warning builds use warnings-as-errors, and CUDA builds have no
      unexpected host/device call diagnostics or an unexplained warning allowlist.
- [ ] CI artifacts retain test logs, configuration manifests, sanitizer output,
      and reference comparison results.

### QUAL-004 — Align documentation with actual science and support

**Priority:** P1

**Verified problem:** README requirements contradict CMake. `docs/paper.tex` is
an uncompleted template. There is no repository-level explanation of the
Hamiltonian, wavefunction, Ewald convention, finite-size limitations, estimator,
or uncertainty method. Local LaTeX build products are present but are correctly
ignored and absent from `git ls-files`; they should not be described as tracked
repository content.

**Why this matters:** Users cannot reproduce or correctly interpret a simulation
whose model, units, numerical approximations, and support matrix are undocumented.
Template text and contradictory prerequisites make the repository look less
finished than its substantive tests deserve.

**Proposed fix:**

1. Write a model document containing equations, units, spin assumptions,
   boundary conditions, shell/twist convention, Jastrow form, Ewald convention,
   local energy, sampling transition, estimator, and statistical analysis.
2. Explain what is and is not comparable to thermodynamic PZ81 values.
3. Document the exact CPU/CUDA residency contract and CLI behavior.
4. Replace or complete the paper template and keep generated `.aux`, `.log`,
   `.fls`, `.fdb_latexmk`, and `.synctex.gz` files out of version control.
5. Add contribution, citation, benchmark, and release-validation instructions.
6. Preserve the existing MIT license and audit every redistributed dependency,
   reference input deck, golden dataset, and generated artifact for attribution
   and redistribution terms. Add a third-party notice only where the audit shows
   one is required; do not infer compatibility from a repository URL.

**Done when:**

- [ ] A new user can reproduce a documented CPU example and, with suitable
      hardware, its CUDA statistical equivalent.
- [ ] Every published number points to a configuration and provenance manifest.
- [ ] Documentation claims are checked in review against code and CI.
- [ ] A recorded license/attribution audit covers code dependencies and external
      validation data included in a release.

### QUAL-005 — Establish performance baselines early; optimize after correctness gates

**Priority:** P1 for measurement and baseline infrastructure; P2 for optimization

**Verified baseline gap:** The audited current tree and clean reference each
contain 20 explicit `cudaDeviceSynchronize` call sites. This static count
documents a synchronization-heavy structure but is not a runtime slowdown factor.
The current migration does not build end-to-end, so no current CPU/CUDA timing is
yet defensible. A reported “126x slower” experimental run was not accompanied by
the revision, hardware, configuration, timing boundary, or completed-run output
needed to reproduce it and therefore is not an acceptance baseline.

**Why this matters:** Attempted moves per second can improve while statistical
efficiency worsens, and aggressive fusion, reduced precision, or fast math can
change numerical results. The relevant product metric is time to a defensible
uncertainty at unchanged physics. An early baseline is also needed to measure
whether the backend redesign removes the cost it is intended to remove.

**Proposed fix:**

1. Immediately after BUILD-001, add a reproducible benchmark harness and record
   a CPU baseline plus a clearly labelled pre-redesign CUDA baseline for the same
   fixed-parameter workload. If CUDA cannot complete correctly, preserve the
   failure/profile evidence and do not publish a speed ratio.
2. Record revision and tree status, CPU/GPU model, driver, compiler/toolkit,
   build flags, precision, particle and walker counts, warmup and measurement
   lengths, optimizer state, initialization inclusion, timing boundaries, and
   repeated-sample summary. Separate compile/startup time from simulation time.
3. Use Nsight Systems or equivalent evidence to count runtime synchronizations,
   launches, and transfer volume by phase rather than inferring cost from source
   call-site counts.
4. Report effective independent samples per second and time-to-target-error, not
   only kernel throughput. Break down initialization, warmup, proposals,
   determinant work, Jastrow, Ewald, local energy, statistics, transfers,
   rebuilds, and output.

Benchmark at representative combinations such as small closed shells and the
current `N=485` target across walker counts and precision. Determine rather than
assume the best CUDA block size, walkers per device, kernel chunk length,
factorization method, and measurement cadence.

**Done when:**

- [ ] A versioned baseline artifact contains the complete workload, environment,
      timing protocol, raw measurements, and summary; no speed ratio is accepted
      without those fields and a completed correct run.
- [ ] A runtime trace, not merely a static call-site count, establishes the
      synchronization and transfer baseline used by GPU-003.
- [ ] Every optimization is accompanied by unchanged correctness/parity results.
- [ ] Benchmark inputs, hardware, compiler flags, and revisions are stored.
- [ ] Fast-math remains off until each affected function has an error budget and
      scientific validation.

### QUAL-006 — Remove stale, unreachable, and ambiguous operational behavior

**Priority:** P2

**Verified examples:** The progress condition checks `i == measure_steps` inside
a loop whose condition is `i < measure_steps`, so that final branch is
unreachable. The program calls CPU execution units “threads” even where the
scientific concept is an independent walker. Exception output is written to
standard output and exits through `std::exit(-1)`; configuration is loaded before
the top-level `try`, so parsing/validation exceptions bypass that handler anyway.
`git diff --check` currently finds trailing whitespace in
`src/energy_tracking/energy_tracking.cu`.
`VMC_FAST_MATH` is explicitly described as inert and has no consumer outside its
CMake definition, yet enabling it forces FP64 off. Thus it changes precision
without enabling the optimization named by the option. There are 24
`#pragma omp simd` directives in the audited source, 14 immediately preceded by
comments claiming “Not vectorized” (four spell it “Not vectorzied”). A GCC 15.2
release build with `-fopt-info-vec-all` nevertheless reports six of the seven
annotated loops in `energy_tracking.cu` as vectorized; the first loop in the CPU
`update_structure_factors`, which calls `sincos` twice per iteration, is the
exception. The comments therefore are not valid evidence of generated code and
are already stale on the audited toolchain. The clean CPU reference produces
zero project warnings under the strict set recorded in QUAL-003, while the clean
CUDA build produces 16 `CUDA_CALLABLE` host/device diagnostics rather than 17.

**Why this matters:** These issues do not change the Hamiltonian, but collectively
make failures, logs, scripts, and reviews less reliable. Ambiguous terminology is
especially costly when mapping CPU workers and CUDA blocks to scientific walkers.

**Proposed fix:**

1. Separate `walker_count`, CPU worker/thread count, and CUDA block/grid settings
   in types, CLI, logs, and documentation.
2. Fix progress using `i + 1`, a time/rate-limited reporter, and a guaranteed final
   completion event outside the loop.
3. Report fatal diagnostics on standard error, return standard nonzero exit codes,
   and allow RAII cleanup rather than calling `std::exit` from ordinary error
   handling.
4. Remove unreachable code, obsolete comments, and whitespace/style violations
   as independently reviewable cleanup; keep generated artifacts untracked.
5. Replace vectorization claims in comments with compiler-report and benchmark
   evidence for the supported release toolchains. Retain pragmas only where they
   are legal, useful, and covered by numerical parity tests.
6. Remove `CUDA_CALLABLE` from host-only orchestration such as the current
   `metropolis_step`, or make its complete call graph device-safe as part of the
   backend redesign. Fix the source of the warnings rather than suppressing the
   category globally.
7. Track intentionally unfinished features as issues rather than public methods
   whose only implementation is to throw.
8. Remove the inert fast-math option now, or implement it only after QUAL-005's
   error-budget gate. Keep precision selection orthogonal: an optimization flag
   must not silently change FP64 to FP32.

**Done when:**

- [ ] Compiler warnings, `git diff --check`, format, and selected static-analysis
      jobs are clean.
- [ ] Stored vectorization reports and focused benchmarks, not source comments,
      support claims about intended hot loops on each supported toolchain.
- [ ] CUDA compilation emits no unexpected host/device call diagnostics.
- [ ] CLI/log terminology distinguishes scientific walkers from execution
      resources.
- [ ] Progress and failure-path integration tests check output and exit status.
- [ ] Every public CMake option has an observable, tested effect matching its
      name, and conflicting options fail rather than silently rewriting one
      another.

### QUAL-007 — Consolidate build scripts and visualization tooling

**Priority:** P2

**Verified problem:** CMake rejects every native Windows build and requires GCC,
but the repository still ships PowerShell scripts that locate Visual Studio and
attempt MSVC builds. Those scripts use inconsistent hard-coded Visual Studio
paths (`...\18\...` in build/profile and `...\2022\...` in tests). The shell
scripts duplicate option parsing, accept an unrecognized argument as a build
type or directory, and the test script always disables CUDA. README does not
document the complete set of flags implemented by those scripts. The renderer is
also coupled directly to the unversioned binary layout described in QUAL-002.

**Why this matters:** Multiple contradictory entry points make it unclear which
build and test paths are supported, and can give users a false impression that a
platform or CUDA test has been exercised. A visualization tool that independently
guesses file layout will drift from the writer again.

**Proposed fix:**

1. Make checked CMake presets the single source of build/test configuration.
   Retain scripts only as thin, strict wrappers around named presets.
2. Remove native-Windows scripts if native Windows remains unsupported, or add
   actual MSVC support and CI before advertising them. WSL instructions should
   call the Linux presets.
3. Reject unknown wrapper arguments and provide `--help`; do not reinterpret a
   typo as a build type or directory.
4. Give CPU and CUDA tests explicit wrapper/preset names so running CPU tests
   cannot be mistaken for CUDA coverage.
5. Provide one tested result-reader module for the versioned format and have
   `render.py` consume it. Keep Python dependencies pinned and document the
   supported Python environment.

**Done when:**

- [ ] Every supported script maps visibly to a CI-tested preset.
- [ ] Unsupported platforms have no apparently functional build wrapper.
- [ ] Unknown options fail nonzero with a useful message.
- [ ] The renderer uses the same documented schema tested by the C++ writer and
      supports every production precision intentionally retained.

---

## Execution tracking and near-term milestones

The audit found 41 separately actionable items across four release gates. They
are not all prerequisites for starting the GPU backend, and the release gates
are not a useful day-to-day work queue. Before implementation begins, track each
item with an owner, status, dependency list, review/evidence artifact, and a
coarse relative size (`S`, `M`, or `L`) assigned by the implementers. Do not
invent calendar estimates before the current tree builds. Mark uncertain work,
especially GPU-004 alternatives and TEST-003, as a bounded spike with a named
decision output rather than an open-ended implementation task.

### Milestone M0 — Restore an observable baseline

BUILD-001 is the only global implementation gate: the current source cannot
validate any change until it compiles. Complete the mechanical xpu migration and
pin its contract (`BUILD-001`, `DEP-001`), then establish the minimum CPU/CUDA
build presets from BUILD-002. Preserve the clean-reference test and warning
results recorded above so a migration change cannot be mistaken for a preexisting
behavior. As soon as both current targets compile, run the QUAL-005 benchmark
harness and capture profiles even though the CUDA implementation is still
labelled scientifically untrusted and non-resident.

**M0 is complete when:**

- [ ] Current CPU and CUDA configurations compile and run their applicable stock
      tests without reintroducing managed allocation or per-move readbacks as a
      purported migration fix.
- [ ] The pinned xpu revision and backend selection are visible in build output.
- [ ] Strict CPU and CUDA warning baselines and reproducible CPU/CUDA timing and
      trace artifacts exist; an incomplete CUDA run is recorded as a failure,
      not converted into a speed ratio.

### Milestone M1 — Defensible fixed-parameter CPU result

This is the first near-term scientific deliverable, not a public-release claim.
Work on these items may proceed in parallel after M0:

1. Make fixed parameters the default immediately; put the existing automatic
   optimizer behind an explicit experimental opt-in (`STAT-003`). Its complete
   redesign follows PHYS-005 and STAT-001.
2. Repair transactional state, the PZ81 validation oracle, and the periodic
   Jastrow (`PHYS-003`, `PHYS-004`, `PHYS-005`).
3. Make Ewald convergence controlled and reject invalid scientific inputs
   (`PHYS-006`, `CFG-001`).
4. Correct correlation analysis and multi-walker aggregation, and validate
   warmup/measurement cadence (`STAT-001`, `STAT-002`, `STAT-004`).
5. Build the relevant TEST-001 fixed-configuration oracle data, including
   independent finite-difference derivatives, local-energy components, and an
   Ewald convergence study. Preserve the already checked identities through the
   focused tests in PHYS-008.
6. Add PHYS-002's fixed-maximum plus residual-triggered determinant rebuilding
   and its condition-stratified tests before calling a long production run
   credible. This is a robustness safeguard, not evidence that ordinary chains
   fail after 21 or 331 accepted moves.

**M1 is complete when:**

- [ ] A fixed-parameter CPU run passes the independent component/local-energy
      oracles and reports the Ewald controls actually used.
- [ ] Its aggregate and uncertainty follow predeclared estimators and demonstrate
      an autocorrelation/blocking plateau; the default path does not select its
      Jastrow parameter from the same production samples.
- [ ] Periodic-Jastrow regularity and determinant-rebuild tests prove that the
      intended correction/rebuild paths actually execute.
- [ ] The result is accompanied by its validated configuration and enough
      provenance to reproduce the M1 evidence, even if the final release schema
      and documentation are not yet complete.

### Parallel tracks after M0

- **CPU science/reference:** execute M1, then PHYS-007, the full STAT-003
  optimizer redesign, STAT-005, ARCH-003, and remaining property tests.
- **CUDA architecture:** in parallel with M1, define the separate artifacts and
  ownership boundary (`ARCH-001`, `ARCH-002`), fix PHYS-001 using the shared
  derivative definition, and implement/prototype GPU-001 through GPU-003 plus
  the split RNG state in GPU-005. Use QUAL-005 traces to choose chunk and walker
  mappings. Stochastic CPU/CUDA parity cannot be accepted until the M1 CPU
  reference is trusted.
- **Verification/infrastructure:** grow TEST-001 and the reproducible build
  matrix, configuration/provenance, result schema, CI, and documentation without
  blocking independent GPU ownership work. Start TEST-003 only as its bounded
  feasibility/decision spike.

The non-negotiable cross-track dependencies are:

1. PHYS-005 and STAT-001 precede validation of the redesigned optimizer; the
   current boundary-rejection rule is removed, not calibrated around the broken
   ansatz.
2. PHYS-002's rebuild policy precedes the final GPU-004 design. GPU-004 first
   supplies one supported device-resident baseline; a second factorization path
   requires a recorded capability or performance failure.
3. Shared scalar formulas and TEST-001 goldens precede acceptance of TEST-002
   CPU/CUDA parity. TEST-002 may be scaffolded earlier.
4. TEST-003 fixed-configuration feasibility precedes any external stochastic
   campaign; a nonidentical ansatz is not used as a total-energy oracle.
5. The QUAL-005 measurement baseline precedes optimization, while all
   performance-changing optimization remains behind the relevant correctness
   and residency gates.

After M1, Gate A, the GPU implementation, and release infrastructure can advance
concurrently wherever those edges permit. Gates A through D below define what
may be claimed at release, not an instruction to postpone all GPU work until
every CPU, documentation, and external-validation task is finished.

## Release gates

### Gate A — Trusted CPU scientific reference

- [ ] CPU-only clean build has no CUDA dependency.
- [ ] All CPU-applicable P0 physics and statistics defects are fixed with
      regression tests; CUDA-only P0 items remain release-blocking at Gate C.
- [ ] Determinant drift remains controlled for long FP32/FP64 sequences.
- [ ] Ewald and fixed-configuration independent oracles pass.
- [ ] Reported uncertainties demonstrate a blocking/autocorrelation plateau.
- [ ] Fixed-parameter runs are the default; optimizer results are independently
      confirmed when requested.

### Gate B — Backend separation

- [ ] `vmc`, `vmc-cpu`, and `vmc-cuda` satisfy the artifact contract.
- [ ] Backend selection is explicit, tested, and never silently falls back.
- [ ] CPU and CUDA use the same validated model/request/result schemas.
- [ ] The pinned xpu revision and its tested allocation/address-space behavior
      are recorded in build and run metadata.

### Gate C — Trusted CUDA scientific backend

- [ ] CUDA initialization, sampling, energy, and statistics stay device-resident.
- [ ] Component and state parity passes against the CPU reference.
- [ ] Independent full runs agree under predeclared statistical criteria.
- [ ] Device determinant rebuild, RNG stream separation, and restart pass.
- [ ] Compute-sanitizer and residency traces are clean.
- [ ] Capacity checks and representative `N=485` runs succeed on target hardware.

### Gate D — Public/research release

- [ ] CI matrix, format/static checks, sanitizers, documentation, and durable
      output are complete.
- [ ] At least one exact-convention fixed-configuration QMCPACK or CASINO
      comparison is reproducible if TEST-003's feasibility spike establishes
      that one is representable; any stochastic comparison is separately scoped
      and predeclared. If neither code can represent the ansatz exactly, the
      recorded negative result and an independently reviewed replacement
      validation strategy are release requirements.
- [ ] Finite-size/twist limitations and PZ81 comparison rules are documented.
- [ ] Every reported result carries full provenance and a defensible uncertainty.
- [ ] Performance claims use effective samples per second or time-to-error and
      include reproducible benchmark inputs.

## Primary references used to verify this plan

1. J. P. Perdew and A. Zunger, “Self-interaction correction to density-functional
   approximations for many-electron systems,” *Physical Review B* 23, 5048
   (1981): <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.23.5048>.
   The parameter table is also visible in this copy:
   <https://www.colorado.edu/faculty/zunger-matter-by-design/sites/default/files/attached-files/55.pdf>.
2. R. J. Needs et al., *CASINO User's Guide*, for parallel-spin cusp conditions,
   periodic Jastrow cutoff behavior, and the three-dimensional Ewald expression:
   <https://casinoqmc.net/casino_manual_dir/casino_manual.pdf>.
3. QMCPACK manual, VMC methods, for particle-by-particle moves, walker crowds,
   measurement/substep semantics, and periodic determinant recomputation:
   <https://qmcpack.readthedocs.io/en/develop/methods.html>.
4. QMCPACK manual, wavefunction introduction, for periodic and long-range/RPA
   Jastrow constructions:
   <https://qmcpack.readthedocs.io/en/develop/intro_wavefunction.html>.
5. C. Lin, F. H. Zong, and D. M. Ceperley, “Twist-averaged boundary conditions in
   continuum quantum Monte Carlo algorithms,” *Physical Review E* 64, 016702
   (2001):
   <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.64.016702>.
6. H. Flyvbjerg and H. G. Petersen, “Error estimates on averages of correlated
   data,” *The Journal of Chemical Physics* 91, 461–466 (1989):
   <https://doi.org/10.1063/1.457480>.
7. The pinned [UWHPC/xpu](https://github.com/UWHPC/xpu) backend and memory-model
   documentation and its `cudaMalloc` implementation, audited at commit
   `5ab6884` (the sibling checkout on the audit machine contained the same clean
   revision).
8. NVIDIA cuSOLVER documentation, for host-launched dense LU using device matrix,
   pivot, and status storage plus its documented device and host workspaces:
   <https://docs.nvidia.com/cuda/cusolver/>.
9. NVIDIA cuSolverDx LU documentation, for device-side batched LU operators with
   and without pivoting:
   <https://docs.nvidia.com/cuda/cusolverdx/get_started/getrf.html>.
10. NVIDIA cuRAND documentation, for Philox seed, subsequence, offset, state,
    reproducibility semantics, and the legacy-device-API deprecation notice:
    <https://docs.nvidia.com/cuda/curand/>.
11. NVIDIA Compute Sanitizer documentation, for memcheck, shared-memory
    racecheck, global-memory initcheck, and synchronization checking:
    <https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html>.
12. N. D. Drummond et al., “Finite-size errors in continuum quantum Monte Carlo
    calculations,” *Physical Review B* 78, 125106 (2008), including its 2014
    erratum, for interaction/kinetic finite-size corrections and twist-averaging
    limitations: <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.78.125106>.
13. S. Chiesa et al., “Finite-Size Error in Many-Body Simulations with Long-Range
    Interactions,” *Physical Review Letters* 97, 076404 (2006), for long-wavelength
    kinetic and potential finite-size corrections:
    <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.076404>.
