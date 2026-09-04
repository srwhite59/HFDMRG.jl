# Packet 10A5R2 H1000 extraction qualification

Packet 10A5R2 stops at the mandatory compact-UHF performance gate. The
accepted Packet 10A5R1B source remains unchanged and no production repair is
folded into this qualification.

## RHF qualification

The compact typed RHF route was run from the block-native dimer-product start
with Julia 1.12.6, four BLAS threads, `state_maxdim=16`,
`state_cutoff=1e-12`, and `energy_tolerance=1e-11`. It converged in nine
half-sweeps to represented energy `-2227.128386923008 Ha`, differing from the
accepted value by `4.55e-13 Ha`.

The nonlinear phase took `88.9003 s`, or `1.21796` times the frozen PPP
dimer-start sweep time of `72.9910 s`; this passes the `1.25` target. The
8,992 visited centers each used one Fock construction, eigensolution,
post-factorization true-energy evaluation, and publication. The additional
nine preparations and Focks are the half-sweep invariant evaluations.

The reflected and physical update-disabled replay pair changed represented
energy by `4.55e-13 Ha` net. A four-half-sweep stationary nonlinear replay
changed energy by `1.36e-12 Ha`. Initialization ingested 500 dimers once and
later dense-row reads were zero. The one collection occupied 103,615,656
summary bytes, the nonlinear workspace 14,275,200 bytes, and the reusable
lifecycle workspace 13,978,664 bytes. Peak RSS was 2,023,112,704 bytes.

One streamed producer-energy evaluation returned
`-425.61104061526885 Ha` in `11.7050 s` with 160 allocated bytes. Its
`6.96e-9 Ha` rounding envelope contains the `1.91e-11 Ha` difference from the
accepted producer record, and serialized compact-state hashes were identical
before and after evaluation. No producer residual was evaluated. Represented
and producer energies remain distinct contracts.

## UHF hard stop

The compact typed UHF route used the accepted block-native atomic-Neel start
and the same controls. It reached represented energy
`-2282.2570732618165 Ha` in 12 half-sweeps, differing from the accepted value
by `3.18e-12 Ha`. Its nonlinear phase took `289.4157 s`, however, or
`2.33861` times frozen PPP's `123.7554 s`. This exceeds the mandatory `1.5`
hard-stop ratio.

The runner stopped immediately after classifying that timing. It did not run
UHF stationary, spin-swapped, or producer-energy replays, and it did not
serialize UHF allocation, reusable-ownership, call-count, or peak-RSS fields.
Those missing fields are not inferred from the earlier H100 gate, and no
second H1000 run is authorized merely to recover them.

## Static surface and inventory

The source and tests are byte-clean at R1B commit
`3e94fd707f5980f8c522d67351f174445b942155`. The compact production slice is
5,091 physical lines. The historical dense, sliced, cached,
target-residual, frozen-occupied, and history slice is 5,692 physical lines;
the module/timing entry layer is 291 lines. The public surface remains exactly
eight exported names, `solve_hfdmrg` has 19 methods, and Julia's ambiguity
check reports zero ambiguities. The removed `PreparedUHFCenter`,
`_prepare_center_cross!`, `_append_component_cross!`, and obsolete
`interaction_centers.jl` path remain absent. Historical block-local
`pair_field` names are intentionally retained and are not compact total-center
fallbacks.

## Readiness decision

The corrected compact RHF extraction is H1000-qualified. The combined compact
RHF/UHF public transfer is **not ready** because UHF fails the explicit H1000
performance ratio despite numerical same-basin agreement. No claim is made
about producer-Hamiltonian stationarity or a global HF minimum. The next work,
if authorized, must be a separate bounded UHF performance localization or
source-repair packet, not an amendment to R2 evidence.
