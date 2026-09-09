# Hosted CI Frozen-Publication Correction

Date: 2026-09-09
Base candidate: `705201a0230d8548f416793a3a0491041be6a5e1`
Hosted run: `34393522107`, test job `102607451542`
Status: implemented and validated; unstaged and uncommitted

## Cause and correction

The two hosted failures were independent. The slot-drainage fixture has exact
zero energy, but compared against one signed macOS cancellation value with a
`2e-32` absolute tolerance. It now uses the established `256N*eps(Float64)`
roundoff scale multiplied by the one-body operator scale. The target-projector
and drainage/lifecycle assertions remain unchanged, and the returned energy is
also compared with an independent energy reconstruction under that envelope.

The substantive UHF failure was in historical
`frozen_occupied=:roundoff_exact`, introduced at `313149b`, rather than either
compact engine. `_frozen_local!` evaluated the local energy with its temporary
damped density but published the pure `Uup`/`Udn` determinant. The return value
therefore did not certify the returned orbitals whenever the final evaluation
started with damping below one.

The correction records the damping used before each density evaluation. When
the final evaluated density was mixed, it uses the existing frozen scalar,
effective one-body fields, and density-backend Fock operation once more to
evaluate the pure final active projectors. The mixed local-energy history,
rise detection, damping update, local convergence, frozen schedule, and user
controls are unchanged. If a rise lowers damping only after the final energy
evaluation, the recorded evaluation damping remains one and the undamped
energy path is retained.

## Downstream behavior

The corrected value remains the window diagnostic `final_energy`, the
end-of-sweep convergence energy, the observer energy, and the returned tuple
energy. Consequently, a sweep whose terminal window used a mixed density can
make a different and more truthful convergence decision: it now compares the
energy of the pure determinant actually carried between windows. This is the
only demonstrated solver behavior change.

The necessary pure-publication certificate adds one density-backend Fock build
only for affected frozen windows. The existing Fock-count regression now
derives that count from `damping_before` and rises preceding the final local
iteration. No additional work occurs when the final evaluated density was
undamped. Frozen occupied remains a correctness mode with no performance
claim.

## Bounded reproduction and validation

Runtime: Julia 1.12.6 on macOS, ILP64 OpenBLAS
`libopenblas64_.dylib`, four BLAS threads.

The focused ordinary and checked-bounds runner used deterministic N=12 RHF
and UHF cases whose terminal windows entered at damping `0.5`. The results
were:

| Mode | Mixed trial energy | Published pure energy | Reconstruction error |
| --- | ---: | ---: | ---: |
| RHF | -5.037598953815567 | -5.037603542333162 | 8.88e-16 ordinary; 0 checked |
| UHF | -5.1902565053298195 | -5.1902567127380115 | 8.88e-16 ordinary; 3.55e-15 checked |

- Focused ordinary: 11/11 assertions, 11.3 seconds.
- Focused checked-bounds: 11/11 assertions, 16.8 seconds.
- Fresh isolated package test with `--check-bounds=yes`: 2231/2231
  assertions passed. The first attempt exposed only the obsolete Fock-count
  expectation; after accounting for the required conditional publication
  Fock, the complete rerun passed.

Exact isolated command:

```sh
/Users/srw/.juliaup/bin/julialauncher +release --startup-file=no -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path=pwd()); Pkg.test("HFDMRG"; julia_args=["--check-bounds=yes"])'
```

Post-validation SHA-256:

- `src/frozen_occupied.jl`:
  `086ade07ac86abe6f3450e7bd97561641e4670027ea82a9839eb098b7d962dab`
- `test/runtests.jl`:
  `7819d6407d1a81b11d6140df4b4b00b92beb4d2c8a0e718a40508f904a336d69`
- machine-local focused runner:
  `2d13e7dc72144f6427b473c6cdb34ec88e8384fdcf31e060c0152565483f8299`

No compact source, public API, scientific control, dependency, fixture,
external repository, remote, or visibility setting changed.

-- hfdmrg-manager@Mac-Studio
