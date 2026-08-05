# History acceleration

## Status

**Private and provisional.** History acceleration is implemented for
density-density and fixed/ragged sliced RHF/UHF research workflows, but no
history driver or policy is exported. This page records the numerical design;
it deliberately provides no calling syntax and makes no universal performance
claim.

## Purpose and spaces

Slow sweep convergence can reflect one or a few global occupied--virtual
directions that are transmitted only gradually by local windows. History
acceleration intermittently solves HF in a small one-particle model assembled
from recent full-physical-basis determinants, then resumes ordinary HFDMRG.

For current occupied orbitals ``C`` and earlier snapshots ``C_j``, RHF builds
an anchored orthonormal basis

```math
Q=\operatorname{orth}\left[C,\;\mathcal{T}
\left((I-CC^T)[C_1,\ldots,C_{h}] /\sqrt{h}\right)\right],
```

where ``\mathcal{T}`` is a rank-revealing absolute-plus-relative-tail
truncation. The current occupied space is retained exactly. For UHF, the anchor
is the common spatial union of current alpha and beta occupied spaces, and past
spin histories are normalized by the square root of the number of populated
spin channels and snapshots. This removes trivial history-count/spin-count
scaling while retaining genuine differences between alpha and beta histories.

The model rank is ``r=\operatorname{cols}(Q)``.

## Inputs and outputs

The private lifecycle consumes full-basis sweep-end determinants, one supported
interaction representation, and a fixed internal history policy. A proposed
RHF determinant, or paired UHF determinant, is lifted from ``Q`` to the
physical basis and audited before it can seed cleanup sweeps.

The UHF queue stores paired ``(C_\alpha,C_\beta)`` snapshots. One common
spatial `Q` is used, but the contracted alpha and beta densities, Focks, and
eigensolves remain distinct. Acceptance and fallback are atomic across spins.

## Contracted interactions

### Density-density

The pair map

```math
P_{p,(ab)}=Q_{pa}Q_{pb}
```

gives a contracted four-index model

```math
G=\operatorname{reshape}\left(P^TVP\right).
```

The workspace guard bounds the dominant pair-map/model allocation before it is
made; exceeding it cleanly falls back to ordinary sweeps.

### Sliced

For each slice ``s``, define the ordered pair map

```math
R_s[(a,b),(i,j)]=Q_s[a,i]Q_s[b,j].
```

The exact history interaction is

```math
G_Q=\sum_{s,t}R_s^T V_{st}R_t.
```

It is accumulated in bounded row batches. With
``S_2=\sum_s d_s^2``, it stores ``O(r^4)`` model data and at most
``O(S_2r^2)`` pair maps, without constructing a physical ``N^4`` tensor.
Projection and cached sliced backends share this contraction and the same
full-space audit; they differ only in their ordinary cleanup sweeps.

## Contracted HF and DIIS

The contracted solver builds Focks from the contracted density and interaction.
Pulay errors are gauge-invariant density commutators ``FD-DF``. A
rank-revealing, scale-relative affine solve produces coefficients whose sum is
normalized to one. The incoming determinant seeds a tie-aware best-iterate
ledger, so the contracted solve always retains a no-change candidate.

For UHF, one coefficient vector is obtained from the summed alpha/beta error
inner products. The resulting mixed Focks are diagonalized separately. Two
consecutive invalid Pulay solves reject the proposal rather than returning a
partially updated spin pair.

## Physical acceptance and cleanup

Every proposal is recomputed with the original full physical interaction.
Acceptance requires finite energy and orbitals, small Gram errors, nonincreasing
energy up to a scale-aware roundoff tie, and sufficient occupied-projector
overlap. The full physical commutator is not formed in the driver; the audit
uses orbital residuals ``(I-CC^T)FC`` with ``O(N^2n)`` work.

Accepted, rejected, no-op, and resource-fallback proposals all continue through
ordinary HFDMRG cleanup sweeps. Rejection preserves the complete incoming
determinant before cleanup.

**Sliced cleanup always starts a fresh solver initialization.** A global
history proposal can contain occupied directions absent from compressed
environments. Projecting it into stale block bases could silently destroy the
correction, so stale-environment projection is forbidden. Fresh initialization
is specific to these between-segment global proposals; it is not an ordinary
in-sweep global cache rebuild.

## Invariants

- The current occupied projector is contained in `Q` to roundoff.
- `Q` is orthonormal and its rank passes an overflow-safe resource guard.
- UHF spin proposals are accepted or rejected together.
- Density and sliced model Focks/energies agree with their physical
  representations under projection.
- Physical energy is not allowed to rise beyond the scale-aware tie.
- Cleanup never projects a history proposal into stale environments.

## Convergence semantics

The private driver targets a full-space orbital residual, not merely sweep
energy change. Raw finite-tolerance projector displacement is diagnostic; a
separate, tighter stationary closure is used when validating basin identity.
History acceleration helps approach a supplied basin. It is not a
metastability search or proof of the global HF minimum.

## Cost and scaling

For density interactions, model construction is dominated by ``P``, ``VP``,
and ``G`` and is guarded against excessive ``O(Nr^2+r^4)`` workspace. The
sliced model uses bounded ``O(S_2r^2+r^4)`` persistent/model workspace and
batched interaction rows. Contracted eigensolves cost ``O(r^3)``.

Cleanup sweeps usually dominate once ``r\ll N``. For cached sliced data, every
accepted or rejected cycle also pays a fresh block/cache initialization. That
cost is accepted for correctness but must be measured on any proposed
production fixture. Repeated-restart validation has shown cache families become
unreachable rather than accumulating, though allocator high-water RSS can
remain large.

## Code map

- `src/history_accelerated_rhf.jl`: shared policy, normalized RHF basis,
  density model, Pulay solve, physical audit, and lifecycle.
- `src/history_accelerated_uhf.jl`: paired UHF history and coupled Pulay logic.
- `src/history_accelerated_sliced.jl`: sliced contraction, physical audit, and
  sliced cleanup dispatch.

## Current deviations and non-goals

- The policy and drivers are private and may change without compatibility
  guarantees.
- There is no target-residual history model.
- There is no live cache continuation for global proposals.
- History does not expose a public observer schema or checkpoint format.
- The feature is not advertised as an acceleration for every system; a
  naturally slow calculation must establish a performance benefit.
