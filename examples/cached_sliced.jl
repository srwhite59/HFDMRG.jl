using HFDMRG
using LinearAlgebra
using Random

rng = MersenneTwister(31)
ns, nj = 6, 2
layout = HFDMRG.SliceLayout(fill(nj, ns))
N, Nocc = layout.offs[end], 2

onsite = collect(range(-1.0, 1.0; length = N))
hopping = diagm(1 => fill(-0.1, N - 1), -1 => fill(-0.1, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
psi0 = Matrix(qr(randn(rng, N, Nocc)).Q)[:, 1:Nocc]

# A compact, symmetric fixed-slice interaction. V6[a,b,c,d,n,m] stores
# (p q|r s) for p=(n,a), q=(m,c), r=(n,b), and s=(m,d).
pair_shape = [1.0 0.15; 0.15 0.8]
V6 = [pair_shape[a, b] * 0.08 / (1 + abs(n - m)) * pair_shape[c, d]
      for a = 1:nj, b = 1:nj, c = 1:nj, d = 1:nj, n = 1:ns, m = 1:ns]

# The layout convenience overload selects the projection backend. Construct
# the cached backend explicitly for the production sliced route.
backend = HFDMRG.SlicedBasisBackendCached(layout, V6)
psiup, psidn, energy = solve_hfdmrg(H, backend, psi0;
    maxiter = 4, blocksize = 2, cutoff = 1e-9, scf_cutoff = 1e-9)

@assert psiup == psidn
@assert isfinite(energy)
@assert norm(psiup' * psiup - I) < 1e-10
println("cached-sliced RHF energy = ", round(energy; digits = 12))
