using HFDMRG
using LinearAlgebra
using Random

rng = MersenneTwister(7)
N = 12
Nocc = 2

onsite = collect(range(-1.0, 1.0; length = N))
hopping = diagm(1 => fill(-0.1, N - 1), -1 => fill(-0.1, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
V = [0.4 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
psi0 = Matrix(qr(randn(rng, N, Nocc)).Q)[:, 1:Nocc]

psiup, psidn, energy = solve_hfdmrg(H, V, psi0;
    maxiter = 4, blocksize = 2, cutoff = 1e-9,
    scf_cutoff = 1e-9, verbose = false)

@assert psiup == psidn
@assert norm(psiup' * psiup - I) < 1e-10
println("energy = ", round(energy; digits = 12))
