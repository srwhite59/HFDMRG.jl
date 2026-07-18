using HFDMRG
using LinearAlgebra
using Random

rng = MersenneTwister(23)
N, Nup, Ndn = 14, 2, 1
onsite = collect(range(-1.2, 1.0; length = N))
hopping = diagm(1 => fill(-0.08, N - 1), -1 => fill(-0.08, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
V = [0.3 / (1 + abs(i - j)) for i = 1:N, j = 1:N]

psiup0 = Matrix(qr(randn(rng, N, Nup)).Q)[:, 1:Nup]
psidn0 = Matrix(qr(randn(rng, N, Ndn)).Q)[:, 1:Ndn]

target_dim = 2
Q = Matrix(qr(randn(rng, N, target_dim)).Q)[:, 1:target_dim]

# Canonical unnormalized pair order: (1,1), (2,1), (2,2).
# Negative entries are allowed: residual_pair is a signed correction.
residual_pair = [0.020 -0.004 0.002;
                -0.004 -0.010 0.003;
                 0.002 0.003 0.015]
backend = HFDMRG.DensityDensityTargetResidualBackend(V, Q, residual_pair)

psiup, psidn, energy = solve_hfdmrg(H, backend, psiup0, psidn0;
    maxiter = 4, blocksize = 2, cutoff = 1e-9, scf_cutoff = 1e-9)

@assert isfinite(energy)
@assert norm(psiup' * psiup - I) < 1e-10
@assert norm(psidn' * psidn - I) < 1e-10
println("target-residual UHF energy = ", round(energy; digits = 12))
