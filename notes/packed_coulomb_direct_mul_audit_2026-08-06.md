# Packed Coulomb direct-`mul!` audit

Date: 2026-08-06

Status: scratch attribution accepted and implemented as a separately guarded
kernel milestone after packed-exchange promotion.

## Outcome

The two hand-written packed matrix-vector contractions in the compact Coulomb
cross-direct kernel are profitably replaceable by two batched `mul!` calls. On
the authenticated H20 middle window, the complete direct phase fell from
0.65369165 ms to 0.355291675 ms at one BLAS thread, a 45.6% reduction. The
candidate remained zero-allocation after its workspace was preallocated.

This audit followed, and is separate from, exchange implementation commit
`1ba3b87d55d0fc05e17f940227025166bbcf0a44`. No HFDMRG solve was run.

## Candidate algebra

For the exterior slices represented by one block channel, concatenate the
weighted canonical-pair densities into `q`:

```text
q_s[(a,b)] = s(a,b) D_s[a,b].
```

With the active exterior columns of the stored packed channel denoted `W`, the
two current scalar contractions become

```text
block_pair_field = W q
exterior_pair_field = W' block_pair_density.
```

The first result is unpacked into the block Fock field. The second is unpacked
slice by slice and transformed through each exterior slice basis. This is the
same square-root-of-two convention as the current scalar loops.

Only one exterior-pair vector is needed: after `W*q` completes, the same vector
can be overwritten by `W'*block_pair_density`. For the audited 300 width-four
slices it contains 3,000 values, or 24,000 bytes. This is reusable window
scratch, not persistent block/cache storage.

## Numerical and timing evidence

All measurements used Julia 1.12.6 and committed HEAD `1ba3b87`. Times are
medians of five warmed batches of 40 calls.

| Shape | BLAS threads | Current (ms) | Batched (ms) | Ratio | Extra workspace |
|---|---:|---:|---:|---:|---:|
| H20 `kL=20,kR=15,d=4` | 1 | 0.653692 | 0.355292 | 0.5435 | 24,000 B |
| H10-like `kL=10,kR=9,d=4` | 1 | 0.259074 | 0.207511 | 0.8010 | 24,000 B |
| Small `kL=6,kR=5,d=4` | 1 | 0.174806 | 0.161558 | 0.9242 | 24,000 B |
| H20 `kL=20,kR=15,d=4` | 8 | 0.764898 | 0.296693 | 0.3879 | 24,000 B |

The maximum absolute Fock difference was `4.440892098500626e-16`; the maximum
relative error was `3.278906630744624e-16`. Both current and candidate kernels
allocated zero bytes in steady state.

After the accepted exchange change, direct work is no longer dominant. At the
H20 shape it is about 0.65 ms versus about 2.30 ms for fused exchange. The
direct candidate is nevertheless useful and has no observed small-rank
regression.

## Required design work before implementation

1. Add one reusable vector sized to the active exterior canonical-pair columns;
   do not change persistent `G/W` cache storage.
2. Prove the active columns are contiguous for fixed and ragged aligned windows,
   and handle an empty exterior range without relying on invalid endpoints.
3. Cover both sweep directions and unequal ranks; either block can supply the
   packed channel in a future generalized call.
4. Preserve restricted and unrestricted total-density conventions.
5. Compare complete cross time after combining the direct candidate with the
   committed exchange kernel; do not accept an isolated direct improvement that
   regresses the full cross application.
6. Require fixed/ragged projection parity, changing-rank absorption, zero steady
   allocation, small-rank no-regression, and the same authenticated H20 witness.

No H20/H12/H16 solve, direct-kernel source edit, public API change, or cache
lifecycle change is part of this audit.

## Implemented guard and validation

The complete-cross measurement exposed a small-rank boundary that the isolated
direct timing could not establish. The implementation therefore keeps the
ordered direct kernel for left ranks below 12 and uses the two-`mul!` kernel at
rank 12 and above. This is a kernel-level selection inside one compact backend;
it does not construct two cache states or fall back per window. Guarded windows
allocate no exterior-pair scratch.

On the authenticated H20 middle window at one BLAS thread, the implemented
kernel reduced direct time from 0.717 to 0.366 ms (49.0%) and complete cross
time from 3.094 to 2.709 ms (12.4%). At the host's eight-thread BLAS setting,
direct time fell from 0.837 to 0.301 ms (64.1%) and complete cross time from
3.119 to 2.962 ms (5.0%). The maximum direct and complete-cross differences
were `4.440892098500626e-16`, and all compared calls allocated zero bytes in
steady state. The active H20 scratch is 24,000 bytes.

Rank-10 and rank-6 controls take the unchanged ordered arithmetic and allocate
no new window scratch. Package tests explicitly exercise both sides of the
rank-12 guard, unequal left/right ranks, projection parity, and zero steady
allocation. The complete package suite and cached/projection verifier passed.

## Scratch evidence

Root:
`/Users/srw/dmrgtmp/hfdmrg_packed_exchange_audit_a97fdd3_20260806`

- Direct runner SHA-256:
  `37e1aa4107786becaf09110215e83d748a2d245e6f96ddd719b2338c9875c131`
- Direct log SHA-256:
  `50cad1639dc14359f790252ba6b57c5abd95b5d7d6466fde501138c5719a7f07`
- Witness SHA-256:
  `2a64d04123330b4e4420ba2e59434ff3842c04c305d7ab44997ab3565141fd20`
