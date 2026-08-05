# Algorithms

These pages describe the numerical algorithms that the current implementation
executes. They are organized around one contract so that mathematical meaning,
runtime behavior, and code ownership remain visible together.

## Algorithm-page contract

Each algorithm page records:

1. **Status** -- public, private/provisional, or developer-only.
2. **Spaces and dimensions** -- physical, retained, center, and model spaces.
3. **Inputs and outputs** -- including ownership and basis conventions.
4. **Pseudocode** -- the actual lifecycle, not only the final formula.
5. **Linear algebra** -- factorizations, projections, and solve conventions.
6. **Invariants** -- orthogonality, determinant integrity, cache consistency,
   and energy conventions.
7. **Convergence semantics** -- what a local or global stopping test proves and
   what it does not prove.
8. **Cost and scaling** -- work, temporary storage, and retained state.
9. **Code map** -- the implementation files that own the operation.
10. **Current deviations and non-goals** -- provisional behavior and excluded
    generalizations.

The first four algorithm pages cover the moving sweep, local RHF/UHF solves,
interaction backends, and private history acceleration.
