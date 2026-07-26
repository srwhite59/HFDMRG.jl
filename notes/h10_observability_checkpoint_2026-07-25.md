# H10 observability checkpoint

Date: 2026-07-25

Status: branch-only implementation checkpoint. This commit is blocked from
merge by a pre-existing one-body cache coordinate-consistency defect; it is
not an accepted H10 solver result.

The checkpoint adds a private run-local diagnostics collector for environment
spectra, local iteration and damping history, and eigensolver timing, together
with a scratch-oriented H10 preflight runner. The public solver API,
`SweepInfo`, backend API, and diagnostics-off arithmetic are unchanged.

Validation before checkpointing:

- all eight frozen diagnostics-off routes reproduced their baseline energies
  and orbitals exactly;
- the complete package suite passed 267/267;
- the maintained cached/projection verifier passed, with maximum window Fock
  difference `5.107e-15` and maximum sweep-energy difference
  `1.24345e-14 Ha`;
- diagnostics-off allocation was unchanged, and repeated small-fixture wall
  changes ranged from `-1.36%` to `+1.25%`;
- runner SHA-256
  `8ec8e7a836f828a2dec19921cd3d807398907591a0b0fd97038e82c94dd19655`
  passed its smoke-only live-cache count and memory-traversal check.

The frozen natural-order H10 preflight exposed the blocker at LTR window 200,
entry Fock call 0. The live one-body window differed from direct `B' * H * B`
by `3.7484443993207606e-8` maximum and `3.839487842284572e-8` Frobenius.
The Gram error was only `1.3327019905658136e-13`, and cached, projection, and
independent physical interaction-only Focks agreed within
`1.9790197816609273e-12`. Thus the total-Fock discrepancy is inherited from
the one-body window, not from the cached sliced interaction recurrence.

The observability implementation is retained because it made this defect
measurable, but it must remain branch-only until a separately reviewed repair
restores cache coordinates and the unchanged natural and odd-even H10 window
gates pass.
