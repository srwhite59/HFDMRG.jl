# Packet 10A7D compact extensive-energy stopping policy

## Scope

Packet 10A7D is based on
`c5c21cf5aca72d9650387141b989a9bb8feaa489` plus the accepted unstaged
10A7B occupation-default and 10A7C RHF transfer corrections. The six recovered
10A6D documentation deltas remain separate and unchanged.

This packet changes compact RHF/UHF stopping and result diagnostics only. It
does not change local optimizer acceptance, schedules, state or interaction
selection, cutoffs, history, historical backends, residual policy, or
projector convergence requirements.

## Frozen PPP policy

Frozen donor commit:
`0bcd4f11831a98ae4bf4fb3d09eac5968df9d96c`.

The generic PPP RHF and independent-spin UHF stored solvers compare terminal
represented energies against the energy from the previous half-sweep in the
same orientation. Convergence requires two consecutive quiet half-sweeps,
therefore covering opposite orientations. PPP's accepted extensive H34,000
qualification corrected the bare absolute tolerance source-free to

```text
max(requested_energy_tolerance,
    256 * eps(Float64) * max(abs(represented_energy), 1) *
        sqrt(realized_state_rank))
```

The reference is terminal represented energy, not a local-center acceptance
energy and not producer energy. System-size dependence enters through the
extensive energy magnitude; the rank factor is the maximum durable realized
state rank, not center dimension. Finite replay and discarded-weight
certificates remain separate post-convergence safeguards. No producer
residual or mandatory projector threshold participates.

Donor slices, all adapted rather than byte-identical:

- `_ControlledMonotoneLocalUpdates.jl`, lines 1--10, energy-envelope formula,
  SHA-256 `10a3d71bfa83f96c265c02676b1c98cfe20d5ea8ce138798423c37a997cdcd6b`;
- `_ControlledStoredNonlinearSweeps.jl`, lines 616--676, orientation-specific
  energy references and consecutive-quiet semantics, SHA-256
  `c8b1f061e1be289b9105177890fdbc76eb564b4faf903252cee853ddce3c4f8a`;
- `_IndependentSpinStoredNonlinear.jl`, lines 234--285, independent-spin
  orientation and quiet-state orchestration, SHA-256
  `841aa59d4cec806af519fa725d9db960df01ce15f45ae67081b8676cd1199e61`;
- `validation/packet8f8c1/run_official_h34000_uhf.jl`, lines 125--198,
  accepted extensive-energy correction and recorded effective envelope,
  SHA-256 `dc0e7502056f185a260cee6c5cae3a8cebb6def9d37f2994fe65784329950343`.

## Recipient correction

`_energy_stopping_threshold` implements that global policy once for both
compact spins. The solver loops continue comparing same-orientation terminal
energies and require two consecutive quiet half-sweeps. Local monotonicity
continues using its existing center-dimension envelope and is unchanged.

`HFDMRGResult` and `UHFDMRGResult` now expose:

- `requested_energy_tolerance`;
- `effective_energy_tolerance` at the terminal represented energy/rank;
- `stopping_tolerance_reason`, either `:requested_energy_tolerance` or
  `:float64_extensive_energy_floor`.

Existing result fields and `reason` values are preserved. Compatibility
constructors keep existing internal fixed-state/evidence consumers working;
their stopping fields are `NaN`, `NaN`, and `:not_applicable` because no
nonlinear stopping decision occurred.

Exact small-system replay remains gated by its roundoff envelope. Finite
replay remains gated separately by the requested tolerance plus its existing
roundoff envelope, while discarded-weight-derived projector envelopes remain
reported safeguards. No roundoff projector convergence was added.

## Validation

The new permanent tests cover low/normal/high requested values
`1e-7/1e-9/1e-11`, ordinary request-controlled thresholds, an extensive-energy
Float64 floor, exact boundary quiet/nonquiet classification, invalid inputs,
both-spin public diagnostics, genuine two-half-sweep nonconvergence, and exact
small-system replay.

- Focused ordinary B/C/D and compact RHF/UHF nonlinear suite: 1,210 passed.
- Focused checked-bounds B/C/D and compact RHF/UHF nonlinear suite: 1,210
  passed.
- One fresh temporary-environment checked-bounds `Pkg.test()`: 2,195 passed.

The isolated suite used Julia 1.12.6. Its Project and Manifest were created
under the OS temporary directory; no repository `Manifest.toml` exists.
No H100/H1000, physical reproduction, broad timing campaign, producer
residual, or external integration was run. The bounded policy arithmetic,
both-spin solver gates, and complete suite made a physical rerun unnecessary.

## Packet-owned paths

- 10A7D hunks in `src/nonlinear_rhf.jl` and `src/nonlinear_uhf.jl` (both also
  contain accepted 10A7B hunks);
- one 10A7D include in `test/runtests.jl` beside accepted B/C includes;
- `test/compact_stopping_policy.jl`;
- this report;
- `validation/packet10a7d/evidence.toml` and `SHA256SUMS`.

All B/C/D work remains unstaged and uncommitted. Remaining release work is
separate: review/landing of the cohesive corrections, the pending 10A6D guide,
and a distributable pinned physical H-chain fixture/input. Main fast-forward,
remote publication, and clean public-checkout/CI confirmation have not been
performed.
