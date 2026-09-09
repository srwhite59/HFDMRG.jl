# Packet 10A7B: safe compact occupation defaults

Packet 10A7B prevents the compact dimer and atomic-Neel initializers from
interpreting traversal fragments as physical atoms.  The compact RHF and UHF
public solvers and both direct lifecycle initializers now share one preflight.
When `atom_intervals` is omitted, the preflight accepts only phase-one
geometries whose traversal consists entirely of complete 18-row cells.

A partial prefix, partial suffix, or nontrivial first phase is ambiguous
without producer-owned physical atom boundaries.  Such calls fail before
nonlinear workspace construction or lifecycle mutation and direct the caller
to pass ordered, contiguous physical-atom ranges covering the complete
one-body interval.  Explicit partitions remain authoritative: terminal rows
must be assigned to their physical atom, never discarded or counted as a new
electron-bearing atom.

The change does not add an electron-count API or alter selection, nonlinear
updates, schedules, tolerances, historical routes, or either compact
interaction engine.  Integral-cell defaults and existing explicit terminal
partitions retain their result contracts.

Focused Julia 1.12.6 validation passed:

- new occupation-default tests, ordinary: 74/74;
- new occupation-default tests, checked bounds: 74/74;
- existing compact RHF/UHF initialization, terminal, reflection, ownership,
  and allocation subset: 761/761;
- compact RHF/UHF public nonlinear regression subset: 689/689.

No complete package suite or large-chain calculation was run.  The cohesive
isolated suite remains reserved for the later accepted correction sequence.

Recovery boundary: this work is based on recovered commit
`c5c21cf5aca72d9650387141b989a9bb8feaa489` in the persistent worktree.
Six pending Packet 10A6D documentation files remain separate and uncommitted.
The former untracked `validation/packet10a5/` diagnostics were absent from
the orphaned checkout before recovery; their cause is unverified and no
reconstruction was attempted.
