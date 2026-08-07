# Slot-scoped frozen-state port to current main

Date: 2026-08-07

Status: behavioral port acceptance passed on a branch from current main. The
port is not merged or pushed.

Code commit: `8adb4938e66cb9a93cc156fe32a03cd199acc152`

Code tree: `514c9192e38cc82310066489508685e4b9ffa9df`

Base: `main@f344a698c72be2e25a2d6e3c658cdcf6ee1c870d`

Branch: `port/frozen-slot-scoped-main-20260807`

## Outcome

The slot-scoped frozen mechanism now coexists with current-main history
acceleration, packed Coulomb kernels and compact cache code in one tree. The
port preserves the experimental branch's numerical behavior:

- the compact witness freezes the initial, left-grown and right-grown block
  families and drains to the baseline within roundoff;
- the package suite and cached/projection verifier pass;
- every alpha and beta H20 checkpoint hash through the fixed eleven-sweep
  ordinary and frozen budgets is identical to the experimental branch;
- H20 frozen-family rank statistics are identical to the experimental branch;
- H40 reproduces the known fail-closed defect in sweep 2 at the same semantic
  site, `_protected_record`, with the same message,
  `protected subspace partially overlaps a frozen record`.

This is behavioral port evidence. It is not publication timing and does not
repair the H40 defect.

## Corrected fixed-budget protocol

The accepted H20 protocol is bounded, not gate-driven:

- frozen: eight finite-rung sweeps, one `:off` drainage sweep and two ordinary
  cleanup sweeps;
- ordinary: eleven sweeps from the same inputs;
- conventional UHF: reuse the accepted reference at `51.707964959 s`; do not
  rerun it for a code-port check.

Stationarity is recorded after every sweep but is not a stopping rule. The
accepted conventional reference has physical energy
`-36.726609287132206 Ha`, combined commutator norm
`6.147722252165476e-9`, minimum gap `0.36926377325530135 Ha`, and orbital
hashes `e4bdd77f...e6882` / `b8a32f63...085d6`.

An earlier version of the port transaction had already launched a fresh
conventional calculation and allowed ordinary and frozen HFDMRG to continue
beyond sweep 11. Those computations are retained as raw provenance but are
protocol-excluded. No claim below uses them.

## Compact witness and checks

The committed `N=16` witness gives:

| family | records | mean active rank | mean frozen columns | nonzero records |
|---|---:|---:|---:|---:|
| initial right | 6 | 1.0 | 0.5 | 3 |
| sweep left | 6 | 1.0 | 0.5 | 3 |
| sweep right | 6 | 1.0 | 0.5 | 3 |

The drained endpoint is within `2e-32 Ha` of the baseline
`-4.336808689964836e-19 Ha`. It records zero global rebuilds and zero
transition rebuilds.

Validation results:

- package suite: `816/816` passed;
- cached/projection verifier: fixed and ragged RHF/UHF passed, worst Fock
  difference `5.107e-15 Ha` and worst sweep-energy difference `1.066e-14 Ha`;
- `git diff --check`: passed.

## H20 fixed trajectories

The combined physical commutator norm `K` and cumulative solver time are:

| sweep | ordinary K | ordinary time (s) | frozen K | frozen time (s) |
|---:|---:|---:|---:|---:|
| 1 | `7.554e-3` | 24.660 | `2.818e0` | 22.521 |
| 2 | `1.009e-4` | 42.855 | `1.251e-1` | 41.026 |
| 3 | `4.510e-6` | 58.898 | `7.100e-4` | 65.459 |
| 4 | `2.141e-6` | 75.056 | `3.354e-5` | 86.509 |
| 5 | `1.241e-6` | 91.159 | `4.669e-7` | 106.959 |
| 6 | `7.844e-7` | 107.245 | `3.654e-8` | 138.612 |
| 7 | `5.221e-7` | 123.433 | `3.643e-8` | 158.550 |
| 8 | `3.604e-7` | 139.660 | `3.639e-8` | 188.800 |
| 9 | `2.560e-7` | 155.849 | `3.637e-8` | 206.649 |
| 10 | `1.864e-7` | 172.062 | `3.636e-8` | 223.338 |
| 11 | `1.388e-7` | 188.362 | `3.636e-8` | 240.007 |

All 22 alpha/beta checkpoint pairs are byte-identical to the corresponding
experimental-branch checkpoints. Thus the reported-energy sequence, physical
energy, residual trajectory, gaps and Aufbau diagnostics are also identical.
The new times differ because this is a new process on a different code
substrate and were not collected as publication timing.

At sweep 11:

| route | time (s) | allocation (GB) | peak RSS (GB) | K | physical energy (Ha) |
|---|---:|---:|---:|---:|---:|
| ordinary, port | 188.362 | 99.929 | 6.527 | `1.388e-7` | `-36.726609287132370` |
| frozen, port | 240.007 | 206.130 | 6.747 | `3.636e-8` | `-36.726609287132870` |
| ordinary, old branch | 174.983 | 99.909 | 6.535 | `1.388e-7` | identical |
| frozen, old branch | 217.731 | 206.129 | 6.994 | `3.636e-8` | identical |

The port processes were `7.65%` and `10.23%` slower respectively than the old
branch processes. This is recorded as a measurement, not attributed causally
to the packed kernels because there was no warmed order-reversed timing
protocol.

## H20 frozen-family signature

The complete fixed schedule reproduces the experimental-branch family table:

| phase/family | records | mean active rank | mean frozen columns | nonzero records |
|---|---:|---:|---:|---:|
| initial right, `1e-4` | 598 | 6.4498 | 8.0468 | 426 |
| first-rung left | 1,192 | 7.9614 | 4.0990 | 449 |
| first-rung right | 1,192 | 10.2122 | 4.1032 | 432 |
| second-rung left/right | 1,192 each | 16.6300 / 16.7257 | 0 / 0 | 0 / 0 |
| third-rung left/right | 1,192 each | 16.7198 / 16.7232 | 0 / 0 | 0 / 0 |
| fourth-rung left/right | 1,192 each | 16.7232 / 16.7232 | 0 / 0 | 0 / 0 |

There are 30 alpha and 30 beta promotions, zero global rebuilds and zero
transition rebuilds. The ninth sweep records `:off, drained=true`; sweeps 10
and 11 are ordinary cleanup.

## H40 fail-closed signature

The port reproduces sweep 1 byte for byte:

- energy: `-64.42539799177939 Ha`;
- alpha SHA-256:
  `7d7233b30f3df0143f1e83e46290bbc4365a9e96c16247379adc5e60d7cf2d47`;
- beta SHA-256:
  `c5ec43b928cfe18bd8118c2fb657f77832a43b669940f4597ee837e3dc8a3953`.

Sweep 2 then fails closed with
`protected subspace partially overlaps a frozen record` in
`_protected_record` (`src/frozen_occupied.jl:221` in the port). The old branch
reports the same function, message and release/growth call chain; its source
line was 280 before the slot-scoped deletion and port reconciliation.

## Port contents

Relative to `main@f344a69`, the code commit changes seven files,
`+1621/-310` overall (`+1471/-250` source and `+150/-60` tests):

- `src/frozen_occupied.jl`: slot-scoped lifecycle and reachability audit;
- `src/core.jl`: frozen determinant transport, growth and drainage hooks;
- `src/backends/sliced_basis_cached.jl`: frozen fields and scratch re-derived
  against current-main generic cached kernels;
- `src/backends/sliced_basis_cached_coulomb.jl`: compact-cache frozen wrappers;
- `src/sliced_fock_audit.jl`: matrix-free physical audit required by frozen
  semantics;
- `src/HFDMRG.jl`: private audit include;
- `test/runtests.jl`: compact lifecycle and exact Fock-call witnesses.

No ordinary cached-sliced arithmetic, history policy, public interface,
threshold, rung schedule or convergence rule was changed. Main's timing,
subset eigensolve, generic tensor, compact cache, packed direct/exchange,
history, documentation and cleanup implementations win over the branch's
independent historical equivalents.

## Branch disposition

The implementation content through `fc295e1` is superseded for active
development by `8adb493`. The branch copies of timing (`9917143`), generic
tensor repair (`53865cd`), compact cache (`c026921`) and packed direct
(`e51e36e`) are superseded by current-main commits. The earlier F2A/F2C
identity, successor-rebasing and transition experiments are superseded by the
slot-scoped semantics.

The experimental branch can be retired as a code-development branch, but it
should be retained as an archive until its unique design and Track-B records
(`70d4090`, `9843463`, `6423c0a` and their predecessors) are copied or judged
unnecessary. No unique shipping implementation remains there.

## Artifacts

Accepted artifact root:

`/Users/srw/dmrgtmp/hfdmrg_frozen_port_acceptance_8adb493`

`SHA256SUMS` SHA-256:

`cee31dd42ebbd767f4affe09fd4af7dfd2e37633ffd9b03025158f2eb3b58c34`

The manifest verifies the runner, fixed-trajectory audit, H20 bounded records,
frozen lifecycle tables and H40 failure log. Inputs retain the accepted H20
transfer-manifest SHA-256
`c35f3db346800db97672308a93ef3826bebf4de6710171469aedc47955dedb40`.
