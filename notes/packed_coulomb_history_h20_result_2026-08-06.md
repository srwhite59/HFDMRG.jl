# Packed-Coulomb sliced history result

Date: 2026-08-06

Status: private implementation checkpoint; suitable for review, not yet merged.

## Outcome

The compact physical-Coulomb sliced backend can now build history-space RHF and
UHF models in weighted canonical-pair coordinates. This removes resource
failures caused only by the generic ordered-pair representation. The generic
unsymmetrized sliced history path is unchanged.

On the authenticated H20 sweep history, packing admitted all three prescribed
history sets, including ranks 50 and 58 that the ordered representation had
rejected. It did not demonstrate earlier convergence than ordinary compact
HFDMRG: ordinary sweep 6 was already stationary enough that none of the tested
offline proposals established an acceleration advantage. No history set was
tuned, and no further cleanup or end-to-end timing was run.

## Representation and numerical contracts

For history rank `r`, the physical-Coulomb model stores the interaction in
weighted canonical pairs, with `Kpair = r(r+1)/2`, rather than ordered `r^2`
pairs. Slice pair maps use the same square-root-of-two convention as the compact
cache. The contracted direct and exchange Focks remain algebraically equivalent
to the projection-sliced oracle.

The accepted private UHF stopping rule is unchanged:

```text
max(R_alpha, R_beta) <= policy.target
```

where each `R_spin` is the occupied--virtual orbital-residual norm. The driver
also reports the diagnostics

```text
K_alpha = sqrt(2) R_alpha
K_beta  = sqrt(2) R_beta
K       = hypot(K_alpha, K_beta).
```

In the contracted Pulay solve, `commutator_rms = K/sqrt(2)` is the RMS of the
alpha and beta commutator norms. The norm after horizontal concatenation is
`K`. The RMS is a DIIS ranking quantity, not the physical stopping target.

## H20 offline evidence

The source was the authenticated six-sweep compact H20 capture. All proposals
passed energy, basin, route, and selected-window parity audits.

| Captures | Rank | Setup (s) | Live model bytes | max spin residual | combined K | Result |
|---|---:|---:|---:|---:|---:|---|
| 1,2,3 | 50 | 2.745603583 | 13,025,120 | 2.4605301893541886e-6 | 4.788192608130532e-6 | accepted, above target |
| 1,3,5 | 42 | 1.607720667 | 6,537,504 | 5.163371411399857e-7 | 1.0192444305100262e-6 | accepted, target by existing rule |
| 2,4,6 | 58 | 3.369253459 | 23,447,200 | 2.824366278658122e-7 | 5.529279582260403e-7 | accepted, target |

At rank 42, packing reduced setup time from 6.8733 s to 1.6077 s (4.27x)
and live model storage from 24.91 MB to 6.54 MB (73.8%). The process high-water
RSS was 2,168,897,536 bytes.

For comparison under the separately audited common `K` convention, ordinary
sweep 6 had `K = 7.819e-7`; the sweep-5 and `(1,3,5)` values were `1.242e-6`
and `1.019e-6`. Thus the packed model proves backend completeness and removes
artificial resource failures, but this H20 ladder is negative acceleration
evidence rather than a timing success.

## Validation

- Package suite: 667/667 passed.
- Cached/projection verifier: fixed and ragged RHF/UHF window parity passed;
  maximum Fock difference `5.10702591327572e-15 Ha`, maximum sweep-energy
  difference `1.0658141036401503e-14 Ha`.
- Compact tests cover fixed and ragged contraction parity, RHF/UHF workflows,
  common and split one-body Hamiltonians, an empty spin, resource planning, and
  exact preservation of maximum-spin-residual stopping semantics.
- `git diff --check`: passed.

## Immutable scratch evidence

- Root: `/Users/srw/dmrgtmp/hfdmrg_h20_offline_packed_history_8a42618_20260806_v1`
- Runner SHA-256: `cb5219bccbd1c1a989ddaa64f910710a892f079a9eeb5b81749cbfdd1edf13f6`
- Output manifest SHA-256: `5f3c0702aa0352c745561ef0600bbb9c4ac016bd832dec62025466c889ca5e5a`
- Ready marker SHA-256: `860077ec6b93686ac6da49e4132d2ae12476baa1055781788a005e28aefb84a5`
- Log SHA-256: `66b1fc90ae48130e072889d2493045c031063ed93deb04397d898788bbc7889a`

## Limits

This is private backend completeness, not an H20 speed claim. It does not alter
the public API, history policy, sweep core, cached-sliced arithmetic, or generic
unsymmetrized sliced model. Changing UHF stopping from maximum-spin residual to
combined `K` would require a separate design and cross-backend validation.
