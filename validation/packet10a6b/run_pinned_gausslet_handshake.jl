using Dates
using LinearAlgebra
using SHA
using TOML
using GaussletBases

array_sha256(A::Array{Float64}) =
    bytes2hex(sha256(reinterpret(UInt8, vec(A))))

projector_sha256(C::Matrix{Float64}) =
    array_sha256(Matrix{Float64}(C * transpose(C)))

function load_hfdmrg(source::String, name::Symbol)
    namespace = Module(name)
    Base.include(namespace, source)
    getfield(namespace, :HFDMRG)
end

function method_record(hf, adapter)
    types = adapter.solver_mode == :restricted_closed_shell ?
        Tuple{typeof(adapter.hamiltonian),typeof(adapter.interaction),
            typeof(adapter.psiup0)} :
        Tuple{typeof(adapter.hamiltonian),typeof(adapter.interaction),
            typeof(adapter.psiup0),typeof(adapter.psidn0)}
    method = which(hf.solve_hfdmrg, types)
    file, line = Base.functionloc(method)
    Dict("signature" => string(method.sig), "file" => realpath(file),
        "line" => line)
end

function require_snapshot_method(hf, adapter, source, label)
    record = method_record(hf, adapter)
    expected = string(realpath(dirname(source)), Base.Filesystem.path_separator)
    actual = record["file"]
    println(label, "_method_file=", actual)
    startswith(actual, expected) || error(
        "$label method resolved outside frozen snapshot: $actual; " *
        "expected prefix $expected")
end

function run_one(hf, adapter, controls)
    result = Base.invokelatest(run_atomic_injected_angular_hfdmrg_hf, adapter;
        hfmod=hf, controls...)
    record = Dict(
        "energy" => result.energy,
        "psiup_size" => collect(size(result.psiup)),
        "psidn_size" => collect(size(result.psidn)),
        "psiup_sha256" => array_sha256(result.psiup),
        "psidn_sha256" => array_sha256(result.psidn),
        "projector_up_sha256" => projector_sha256(result.psiup),
        "projector_dn_sha256" => projector_sha256(result.psidn),
        "method" => method_record(hf, adapter),
    )
    (; record, result)
end

function comparison(clean, candidate)
    clean_result = clean.result
    candidate_result = candidate.result
    Pup_clean = clean_result.psiup * transpose(clean_result.psiup)
    Pdn_clean = clean_result.psidn * transpose(clean_result.psidn)
    Pup_candidate = candidate_result.psiup * transpose(candidate_result.psiup)
    Pdn_candidate = candidate_result.psidn * transpose(candidate_result.psidn)
    Dict(
        "energy_difference" => candidate.record["energy"] - clean.record["energy"],
        "alpha_projector_frobenius" => norm(Pup_candidate - Pup_clean),
        "beta_projector_frobenius" => norm(Pdn_candidate - Pdn_clean),
        "alpha_coefficients_equal" =>
            clean.record["psiup_sha256"] == candidate.record["psiup_sha256"],
        "beta_coefficients_equal" =>
            clean.record["psidn_sha256"] == candidate.record["psidn_sha256"],
        "alpha_projectors_equal" => clean.record["projector_up_sha256"] ==
            candidate.record["projector_up_sha256"],
        "beta_projectors_equal" => clean.record["projector_dn_sha256"] ==
            candidate.record["projector_dn_sha256"],
    )
end

function adapter_inputs(adapter)
    Dict(
        "H_size" => collect(size(adapter.hamiltonian)),
        "V_size" => collect(size(adapter.interaction)),
        "psiup0_size" => collect(size(adapter.psiup0)),
        "psidn0_size" => collect(size(adapter.psidn0)),
        "H_sha256" => array_sha256(adapter.hamiltonian),
        "V_sha256" => array_sha256(adapter.interaction),
        "psiup0_sha256" => array_sha256(adapter.psiup0),
        "psidn0_sha256" => array_sha256(adapter.psidn0),
    )
end

function main()
    length(ARGS) == 4 || error(
        "usage: run_pinned_gausslet_handshake.jl CLEAN_SOURCE " *
        "CANDIDATE_SOURCE GAUSSLET_ROOT OUTPUT")
    clean_source, candidate_source, gausslet_root, output = abspath.(ARGS)
    println("pid=", getpid())
    haskey(ENV, "HFDMRG_PID_FILE") &&
        write(ENV["HFDMRG_PID_FILE"], string(getpid(), '\n'))
    realpath(dirname(Base.active_project())) == realpath(gausslet_root) ||
        error("GaussletBases snapshot must be the explicit active project")
    clean_hf = load_hfdmrg(clean_source, :Packet10A6BClean)
    candidate_hf = load_hfdmrg(candidate_source, :Packet10A6BCandidate)

    rb = build_basis(RadialBasisSpec(:G10; count=6,
        mapping=AsinhMapping(c=0.15, s=0.15), reference_spacing=1.0,
        tails=3, odd_even_kmax=2, xgaussians=[XGaussian(alpha=0.2)]))
    grid = radial_quadrature(rb; accuracy=:medium, quadrature_rmax=12.0)
    radial_ops = atomic_operators(rb, grid; Z=2.0, lmax=2)
    cases = Dict{String,Any}()
    for (label, nup, ndn, mode, blocksize, expected_energy) in (
            ("rhf", 1, 1, :restricted_closed_shell, 16,
                13.02625601311897),
            ("uhf", 2, 1, :unrestricted, 64,
                39.38381696539275))
        adapter = build_atomic_injected_angular_hfdmrg_hf_adapter(radial_ops;
            nup, ndn)
        adapter.solver_mode == mode || error("unexpected $label solver mode")
        controls = (nblockcenter=1, maxiter=1, blocksize, cutoff=0.0,
            scf_cutoff=Inf, verbose=false)
        require_snapshot_method(clean_hf, adapter, clean_source,
            "clean_$label")
        require_snapshot_method(candidate_hf, adapter, candidate_source,
            "candidate_$label")
        clean = run_one(clean_hf, adapter, controls)
        candidate = run_one(candidate_hf, adapter, controls)
        cmp = comparison(clean, candidate)
        cmp["energy_difference"] == 0.0 || error("$label energy changed")
        candidate.record["energy"] == expected_energy || error(
            string(label, " accepted energy changed: ",
                candidate.record["energy"], " != ", expected_energy))
        all(cmp[key] for key in ("alpha_coefficients_equal",
                "beta_coefficients_equal", "alpha_projectors_equal",
                "beta_projectors_equal")) || error("$label hashes changed")
        cases[label] = Dict(
            "fixture" => Dict("nup" => nup, "ndn" => ndn,
                "solver_mode" => string(mode)),
            "controls" =>
                Dict(string(k) => v for (k, v) in pairs(controls)),
            "inputs" => adapter_inputs(adapter),
            "clean" => clean.record,
            "candidate" => candidate.record,
            "comparison" => cmp,
        )
    end

    record = Dict(
        "schema" => "hfdmrg-packet10a6b-pinned-gausslet-v1",
        "generated_at" => string(now()),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "gausslet_commit" => "891bf3e84653197b90b930d7d5787be17d8fb998",
        "gausslet_root" => realpath(gausslet_root),
        "clean_hfdmrg_commit" =>
            "19e11fc2b12a142138dc3324f585d8aa92d46098",
        "clean_hfdmrg_source" => realpath(clean_source),
        "candidate_hfdmrg_commit" =>
            "f6b2a24d3bbc96713146c2228dd54666ab20a6bc",
        "candidate_hfdmrg_source" => realpath(candidate_source),
        "basis" => Dict("family" => "G10", "radial_count" => 6,
            "mapping_c" => 0.15, "mapping_s" => 0.15,
            "reference_spacing" => 1.0, "tails" => 3,
            "odd_even_kmax" => 2, "xgaussian_alpha" => 0.2,
            "quadrature_accuracy" => "medium", "quadrature_rmax" => 12.0,
            "nuclear_charge" => 2.0, "lmax" => 2),
        "cases" => cases,
    )
    open(output, "w") do io
        TOML.print(io, record; sorted=true)
    end
    println("output=", output)
    println("rhf_energy=", cases["rhf"]["candidate"]["energy"])
    println("uhf_energy=", cases["uhf"]["candidate"]["energy"])
    println("source_paths_verified=true")
end

main()
