using HFDMRG
using LinearAlgebra
using Random
using Serialization

function write_checkpoint(path, info)
    payload = (sweep = info.sweep, energy = info.energy,
        converged = info.converged, psiup = copy(info.psiup),
        psidn = copy(info.psidn))
    temporary = path * ".tmp"
    open(temporary, "w") do io
        serialize(io, payload)
    end
    mv(temporary, path; force = true)
    nothing
end

rng = MersenneTwister(17)
N, Nocc = 12, 2
onsite = collect(range(-1.0, 1.0; length = N))
hopping = diagm(1 => fill(-0.1, N - 1), -1 => fill(-0.1, N - 1))
H = Matrix(Diagonal(onsite) + hopping)
V = [0.4 / (1 + abs(i - j)) for i = 1:N, j = 1:N]
psi0 = Matrix(qr(randn(rng, N, Nocc)).Q)[:, 1:Nocc]

# This temporary directory keeps the example self-cleaning. For a long run,
# replace it with machine-local scratch and choose an archival format suited
# to the surrounding project.
mktempdir() do directory
    checkpoint = joinpath(directory, "latest_checkpoint.jls")
    history = NamedTuple[]

    function observer(info)
        push!(history, (sweep = info.sweep, energy = info.energy,
            converged = info.converged))
        write_checkpoint(checkpoint, info)
        info.sweep == 3 # Demonstrate a caller-controlled clean stop.
    end

    psiup, psidn, energy = solve_hfdmrg(H, V, psi0;
        # A zero sweep cutoff keeps the built-in stop false in this example.
        maxiter = 10, blocksize = 2, cutoff = 0.0,
        scf_cutoff = 1e-9, observer)

    saved = open(checkpoint) do io
        deserialize(io)
    end
    @assert saved.sweep == 3 == length(history)
    @assert !saved.converged
    @assert saved.energy == energy
    @assert saved.psiup == psiup && saved.psidn == psidn
    @assert norm(psiup' * psiup - I) < 1e-10

    _, _, restart_energy = solve_hfdmrg(H, V, saved.psiup;
        maxiter = 1, blocksize = 2, cutoff = 0.0, scf_cutoff = 1e-9)
    @assert isfinite(restart_energy)

    println("clean stop after sweep ", saved.sweep)
    println("checkpoint energy = ", round(saved.energy; digits = 12))
    println("one-sweep restart energy = ", round(restart_energy; digits = 12))
end
