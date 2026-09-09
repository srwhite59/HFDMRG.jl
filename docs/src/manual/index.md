# Manual

The manual covers the supported user workflow: preparing inputs, selecting a
backend, choosing solver controls, observing complete sweeps, and interpreting
the returned determinant and energy.

HFDMRG expects an ordered real one-particle basis. A calculation fixes the
alpha and beta occupations from the column counts of the initial orbital
matrices. Different initial determinants can converge to different
Hartree--Fock stationary points, especially for UHF.

The recommended reading order is:

1. [Getting started](@ref)
2. [Solver controls](@ref)
3. [Numerical conventions](@ref)
4. [Paper-release reproducibility](@ref)
5. [Interaction backends](@ref)

For the internal algorithm, continue with [Moving-window sweep](@ref). For
extension work, read [Backend architecture](@ref) before changing an
interaction route.

For a self-contained physical compact calculation, run the distributed H10 or
H20 hydrogen-chain fixture described in [Paper-release reproducibility](@ref).
It requires no PPP or GaussletBases checkout at execution time.
