# Solver controls

## Historical dense and sliced routes

The historical matrix-based `solve_hfdmrg` routes use these moving-window controls.

| Keyword | Default | Meaning |
|---|---:|---|
| `blocksize` | `200` | Requested number of physical functions in an interior chunk. This is not a DMRG bond dimension. Small systems may use a smaller effective size to form enough chunks. |
| `block_partition` | `nothing` | Optional `SliceLayout` giving one physical slice per chunk. It overrides `blocksize`; when the backend has a layout, both layouts must agree. |
| `nblockcenter` | `1` | Number of consecutive chunks kept explicit in the moving center. Larger centers enlarge the local variational space but increase local work. |
| `maxiter` | `1000` | Positive integer maximum number of complete left-to-right plus right-to-left sweeps. |
| `cutoff` | `1e-11` | End-of-sweep energy stopping tolerance. Common-H routes use an absolute comparison; genuinely split-H routes scale it by `max(1, abs(E))`. |
| `scf_cutoff` | `nothing` | Local SCF energy tolerance. `nothing` means `cutoff` for common H and `cutoff/10` for split H. Each window allows at most four microiterations. |
| `environment_cutoff` | `1e-10` | Absolute singular-value cutoff for grown occupied-orbital environment bases. Every singular value strictly above the cutoff is retained, with at least one vector. |
| `frozen_occupied` | `:off` | The correctness-oriented density-only `:roundoff_exact` mode is available, but is not a recommended general performance optimization. |
| `observer` | `nothing` | Called after each complete sweep. Return `true` to stop cleanly; return `false` or `nothing` to continue. |
| `verbose` | `false` | Print partition and sweep-energy progress. |

`environment_cutoff` must be finite and nonnegative. It controls environment
compression, not sweep convergence. Lower values retain more restricted
occupied-orbital directions. `cutoff`, `scf_cutoff`, and
`environment_cutoff` should therefore be varied and interpreted separately.

## Center geometry

`nblockcenter=1` is the current default. With more than one explicit chunk,
adjacent windows overlap. That overlap can transmit orbital corrections more
quickly per nominal sweep, but it is not universally faster: larger local Fock
and eigensolve work, and the number of windows, determine total cost. Treat
`nblockcenter` as a scientific control and compare both stationarity and work.

Cached sliced calculations should pass `block_partition=layout` so windows
remain aligned with interaction slices. A split-slice cached window delegates
its complete interaction work to the projection route.

## Convergence and stationarity

HFDMRG's built-in convergence flag is an end-of-sweep energy-change test.
The first completed sweep supplies the baseline for the second.
Reaching `maxiter` without satisfying it is not an exception: the solver
returns the last completed-sweep determinant. Conversely, a small energy
change alone is not proof of Fock stationarity, Aufbau ordering, uniqueness, or
the global Hartree--Fock minimum.

For important calculations, also check:

- alpha and beta Gram errors;
- a full-space occupied--virtual Fock residual;
- occupied--virtual gaps or an appropriate stability diagnostic;
- sensitivity to the initial determinant, center geometry, and tolerances; and
- problem-specific spin or symmetry quantities.

See [Local RHF and UHF](@ref) for the distinction between local energy stopping
and physical stationarity.

## Complete-sweep observation

An observer receives a qualified `HFDMRG.SweepInfo` after each complete sweep.
The orbitals correspond to the determinant associated with the reported
end-of-sweep energy evaluated on those orbitals.

```julia
last_info = Ref{Any}(nothing)
observer = info -> begin
    last_info[] = (sweep = info.sweep, energy = info.energy,
        converged = info.converged, psiup = copy(info.psiup),
        psidn = copy(info.psidn))
    false
end

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0, psidn0;
    maxiter = 20, observer)
```

Observer orbital arrays are read-only and must be copied before retaining or
mutating them. Observer exceptions propagate. Checkpoint formats, scientific
measurements, and user stop rules belong in the caller; their time is included
in solver wall time.

## Troubleshooting

| Symptom | First checks |
|---|---|
| Energy oscillates or rises | End-of-sweep energy is not guaranteed monotone. Inspect the sweep history, initial determinant, and center size. Local energy rises reduce that window's later density-mixing weight. |
| `maxiter` is reached | Examine both energy change and a Fock residual before extending the run. A flat residual may need a different initial determinant or geometry. |
| Initial guesses reach different energies | HF can have multiple stationary determinants. Match occupations and tolerances, then compare stationary energies and physical diagnostics. |
| Returned orbitals lose orthogonality | Check the input Gram matrices and that the observer did not mutate its arrays. Large errors indicate invalid input or numerical failure. |
| Recomputed and returned energies differ | During damped, unconverged local SCF the internal density can be non-idempotent. Require final stationarity and use the same interaction and scalar-offset convention. |
| A sliced calculation is unexpectedly expensive | Use fixed `V6` for uniform slices, construct the cached backend explicitly, align `block_partition`, and account for observer work. |

## Compact finite state selection

On each compact cut, `state_maxdim` caps the retained active rank. Let
``P_\sigma`` denote the product of retained probabilities in one spin
sector. A mode made filled contributes its occupation, a mode made empty
contributes one minus its occupation, and active modes contribute one.
Finite UHF applies `state_cutoff` to ``1-P_\sigma`` independently in
each spin. Finite RHF applies the cutoff to the closed-shell loss
``1-P^2``. When the cap prevents reaching the cutoff, the selection uses
the lowest loss among permitted whole-cluster choices.

In v0.1.0 finite UHF used ``1-P_\sigma^2`` in each spin. The
`UHFDMRGResult.maximum_alpha_discarded_squared_weight` and
`maximum_beta_discarded_squared_weight` field names remain the same.
They report the largest published nonlinear loss in their respective
spins using ``1-P_\sigma``. The associated per-spin certificates during
initialization and replay use the same probability convention. The same
numeric cutoff can therefore yield different UHF ranks and energies.
`exact=true` retains the existing exact-selection behavior.
