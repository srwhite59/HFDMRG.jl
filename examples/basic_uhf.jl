using HFDMRG
using LinearAlgebra

# A small, deterministic open-shell model. Alpha and beta have independent
# occupations and start from different occupied subspaces.
N, Nup, Ndn = 14, 3, 2
onsite = collect(range(-1.3, 1.1; length = N))
hopping = diagm(1 => fill(-0.12, N - 1), -1 => fill(-0.12, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
V = [0.28 / (1 + abs(i - j)) for i = 1:N, j = 1:N]

one_body_orbitals = eigen(Symmetric(H)).vectors
for orbital in eachcol(one_body_orbitals)
    pivot = argmax(abs.(orbital))
    orbital[pivot] < 0 && (orbital .*= -1)
end
psiup0 = hcat(one_body_orbitals[:, 1], one_body_orbitals[:, 2],
    cos(0.38) * one_body_orbitals[:, 3] +
    sin(0.38) * one_body_orbitals[:, 7])
psidn0 = hcat(
    cos(0.31) * one_body_orbitals[:, 1] -
    sin(0.31) * one_body_orbitals[:, 6],
    one_body_orbitals[:, 2])

psiup, psidn, energy = solve_hfdmrg(H, V, psiup0, psidn0;
    maxiter = 80, blocksize = 2, cutoff = 1e-13,
    scf_cutoff = 1e-13, environment_cutoff = 1e-12)

# Reconstruct the complete physical-basis UHF quantities. For the
# density-density convention, direct terms use the total diagonal density and
# exchange uses the density of the same spin.
Dup, Ddn = psiup * psiup', psidn * psidn'
direct = V * diag(Dup + Ddn)
Fup = H + Diagonal(direct) - V .* Dup
Fdn = H + Diagonal(direct) - V .* Ddn
recomputed_energy =
    0.5 * dot(Dup, Fup + H) + 0.5 * dot(Ddn, Fdn + H)

function spin_audit(F, C)
    nocc = size(C, 2)
    diagonalization = eigen(Symmetric(F))
    occupied = diagonalization.vectors[:, 1:nocc]
    aufbau_projector = occupied * occupied'
    residual = F * C - C * (C' * F * C)
    (; residual = norm(residual),
       gap = diagonalization.values[nocc + 1] -
             diagonalization.values[nocc],
       aufbau_error = norm(C * C' - aufbau_projector))
end

up = spin_audit(Fup, psiup)
dn = spin_audit(Fdn, psidn)
gram_up = norm(psiup' * psiup - I)
gram_dn = norm(psidn' * psidn - I)

@assert gram_up < 1e-11 && gram_dn < 1e-11
@assert abs(recomputed_energy - energy) < 5e-12
@assert up.residual < 1e-8 && dn.residual < 1e-8
@assert up.gap > 1e-2 && dn.gap > 1e-2
@assert up.aufbau_error < 2e-8 && dn.aufbau_error < 2e-8

println("UHF energy = ", round(energy; digits = 12))
println("residuals (alpha, beta) = ", (up.residual, dn.residual))
println("occupied-virtual gaps = ", (up.gap, dn.gap))
