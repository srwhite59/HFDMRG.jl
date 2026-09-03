using Dates
using LinearAlgebra
using SHA
using TOML
using GaussletBases

function array_sha256(A::Array{Float64})
    bytes2hex(sha256(reinterpret(UInt8, vec(A))))
end

function projector_sha256(C::Matrix{Float64})
    array_sha256(Matrix{Float64}(C * transpose(C)))
end

function load_hfdmrg(source::String, name::Symbol)
    namespace = Module(name)
    Base.include(namespace, source)
    getfield(namespace, :HFDMRG)
end

function method_record(hf, H, V, psiup0, psidn0)
    method = which(hf.solve_hfdmrg,
        Tuple{typeof(H),typeof(V),typeof(psiup0),typeof(psidn0)})
    location = Base.functionloc(method)
    Dict(
        "signature" => string(method.sig),
        "file" => string(first(location)),
        "line" => last(location),
    )
end

function run_one(hf, adapter, controls)
    result = Base.invokelatest(run_atomic_injected_angular_hfdmrg_hf, adapter;
        hfmod=hf,
        nblockcenter=controls.nblockcenter,
        maxiter=controls.maxiter,
        blocksize=controls.blocksize,
        cutoff=controls.cutoff,
        scf_cutoff=controls.scf_cutoff,
        verbose=controls.verbose)
    record = Dict(
        "energy" => result.energy,
        "psiup_size" => collect(size(result.psiup)),
        "psidn_size" => collect(size(result.psidn)),
        "psiup_sha256" => array_sha256(result.psiup),
        "psidn_sha256" => array_sha256(result.psidn),
        "projector_up_sha256" => projector_sha256(result.psiup),
        "projector_dn_sha256" => projector_sha256(result.psidn),
    )
    (record=record, result=result)
end

function main()
    length(ARGS) == 3 || error(
        "usage: run_matched_gausslet_uhf.jl CLEAN_SOURCE CANDIDATE_SOURCE OUTPUT")
    clean_source, candidate_source, output = abspath.(ARGS)
    clean_hf = load_hfdmrg(clean_source, :Packet10A3R0Clean)
    candidate_hf = load_hfdmrg(candidate_source, :Packet10A3R0Candidate)

    rb = build_basis(RadialBasisSpec(:G10;
        count=6,
        mapping=AsinhMapping(c=0.15, s=0.15),
        reference_spacing=1.0,
        tails=3,
        odd_even_kmax=2,
        xgaussians=[XGaussian(alpha=0.2)]))
    grid = radial_quadrature(rb; accuracy=:medium, quadrature_rmax=12.0)
    radial_ops = atomic_operators(rb, grid; Z=2.0, lmax=2)
    adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(radial_ops;
        nup=2, ndn=1)
    controls = (
        nblockcenter=1,
        maxiter=1,
        blocksize=min(size(adapter.hamiltonian, 1), 64),
        cutoff=0.0,
        scf_cutoff=Inf,
        verbose=false,
    )

    clean_run = run_one(clean_hf, adapter, controls)
    candidate_run = run_one(candidate_hf, adapter, controls)
    clean = clean_run.record
    candidate = candidate_run.record
    clean_up = clean["projector_up_sha256"]
    clean_dn = clean["projector_dn_sha256"]

    clean_result = clean_run.result
    candidate_result = candidate_run.result
    Pup_clean = clean_result.psiup * transpose(clean_result.psiup)
    Pdn_clean = clean_result.psidn * transpose(clean_result.psidn)
    Pup_candidate = candidate_result.psiup * transpose(candidate_result.psiup)
    Pdn_candidate = candidate_result.psidn * transpose(candidate_result.psidn)

    record = Dict(
        "schema" => "hfdmrg-packet10a3r0-gausslet-uhf-v1",
        "generated_at" => string(now()),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "gausslet_source_commit" => "891bf3e84653197b90b930d7d5787be17d8fb998",
        "gausslet_source_file" => "src/angular/angular_atomic_benchmark.jl",
        "gausslet_function" => "run_atomic_injected_angular_hfdmrg_hf",
        "fixture" => Dict(
            "basis_family" => "G10",
            "radial_count" => 6,
            "mapping_c" => 0.15,
            "mapping_s" => 0.15,
            "reference_spacing" => 1.0,
            "tails" => 3,
            "odd_even_kmax" => 2,
            "xgaussian_alpha" => 0.2,
            "quadrature_accuracy" => "medium",
            "quadrature_rmax" => 12.0,
            "nuclear_charge" => 2.0,
            "lmax" => 2,
            "nup" => adapter.nup,
            "ndn" => adapter.ndn,
            "route" => string(adapter.route),
            "solver_mode" => string(adapter.solver_mode),
        ),
        "inputs" => Dict(
            "H_size" => collect(size(adapter.hamiltonian)),
            "V_size" => collect(size(adapter.interaction)),
            "psiup0_size" => collect(size(adapter.psiup0)),
            "psidn0_size" => collect(size(adapter.psidn0)),
            "H_sha256" => array_sha256(adapter.hamiltonian),
            "V_sha256" => array_sha256(adapter.interaction),
            "psiup0_sha256" => array_sha256(adapter.psiup0),
            "psidn0_sha256" => array_sha256(adapter.psidn0),
        ),
        "controls" => Dict(string(k) => v for (k, v) in pairs(controls)),
        "clean" => merge(clean, Dict(
            "source" => clean_source,
            "commit" => "19e11fc2b12a142138dc3324f585d8aa92d46098",
            "method" => method_record(clean_hf, adapter.hamiltonian,
                adapter.interaction, adapter.psiup0, adapter.psidn0),
        )),
        "candidate" => merge(candidate, Dict(
            "source" => candidate_source,
            "method" => method_record(candidate_hf, adapter.hamiltonian,
                adapter.interaction, adapter.psiup0, adapter.psidn0),
        )),
        "comparison" => Dict(
            "energy_difference" => candidate["energy"] - clean["energy"],
            "alpha_projector_frobenius" => norm(Pup_candidate - Pup_clean),
            "beta_projector_frobenius" => norm(Pdn_candidate - Pdn_clean),
            "alpha_projector_hash_equal" => clean_up == candidate["projector_up_sha256"],
            "beta_projector_hash_equal" => clean_dn == candidate["projector_dn_sha256"],
        ),
    )
    open(output, "w") do io
        TOML.print(io, record; sorted=true)
    end
    println("output=", output)
    println("clean_energy=", clean["energy"])
    println("candidate_energy=", candidate["energy"])
    println("energy_difference=", record["comparison"]["energy_difference"])
    println("alpha_projector_frobenius=", record["comparison"]["alpha_projector_frobenius"])
    println("beta_projector_frobenius=", record["comparison"]["beta_projector_frobenius"])
end

main()
