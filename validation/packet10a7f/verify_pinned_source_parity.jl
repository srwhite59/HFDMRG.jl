source_path = joinpath(ARGS[1], "validation", "packet8f7a",
    "run_block_native_product.jl")
source = read(source_path, String)
main_start = findfirst("\nfunction product_main()", source)
main_start === nothing && error("frozen PPP product fixture changed")
include_string(Main, source[1:first(main_start)-1], source_path)

pushfirst!(LOAD_PATH, abspath(ARGS[3]))
using HFDMRG
using GaussletBases: sliced_h1, sliced_h1_bandwidth, sliced_vee
include(joinpath(ARGS[3], "examples", "hydrogen_chain_fixture.jl"))
using .HydrogenChainFixtures: load_hydrogen_chain_fixture

const PDM = DensityMPOLab.PeriodicDensityMPOs

function verify_case(atoms, manifest)
    chain, model, periodic, _ = fixture(atoms)
    loaded = load_hydrogen_chain_fixture(Symbol("h$atoms");
        manifest_path=manifest)
    loaded.one_body.bands == model.one_body.bands || error(
        "loaded H1 bands differ from pinned adapter")
    loaded.interaction.phase_at_first == periodic.phase_at_first || error(
        "loaded first phase differs")
    loaded.interaction.tail_transfer == periodic.tail_transfer || error(
        "loaded tail transfer differs")
    loaded.interaction.residual_transfer == periodic.residual_transfer ||
        error("loaded residual transfer differs")
    loaded.interaction.phase_source == periodic.phase_source || error(
        "loaded phase source differs")
    loaded.interaction.phase_sink == periodic.phase_sink || error(
        "loaded phase sink differs")
    loaded.interaction.local_triangle == periodic.local_triangle || error(
        "loaded local interaction differs")
    loaded.atom_intervals == local_atomic_intervals(model.layout) || error(
        "loaded atom intervals differ")
    loaded.producer.constant == model.nuclear_repulsion || error(
        "loaded nuclear constant differs")
    h1 = sliced_h1(chain)
    vee = sliced_vee(chain)
    sliced_h1_bandwidth(chain) == loaded.producer.h1_bandwidth || error(
        "loaded producer H1 bandwidth differs")
    source_workspace = PDM.PeriodicEntryWorkspace(periodic)
    target_rank = size(loaded.interaction.residual_transfer, 1)
    first_scratch = zeros(target_rank)
    second_scratch = similar(first_scratch)
    for second in 1:loaded.interaction.sites, first in 1:loaded.interaction.sites
        loaded.producer.one_body[first, second] == h1[first, second] || error(
            "loaded producer H1 entry differs at ($first,$second)")
        loaded.producer.interaction[first, second] == vee[first, second] ||
            error("loaded producer Vee entry differs at ($first,$second)")
    end
    for second in 1:loaded.interaction.sites, first in 1:second
        represented = HFDMRG._interaction_entry(loaded.interaction, first,
            second, first_scratch, second_scratch)
        source_value = PDM.interaction(periodic, first, second,
            source_workspace)
        represented == source_value || error(
            "loaded compact interaction differs at ($first,$second)")
    end
    println("h", atoms, "_source_parity=true sites=",
        loaded.interaction.sites, " bytes=", loaded.metadata["data_bytes"])
end

function main()
    length(ARGS) == 3 || error(
        "usage: verify_pinned_source_parity.jl PPP_ROOT MANIFEST HFDMRG_ROOT")
    manifest = abspath(ARGS[2])
    verify_case(10, manifest)
    verify_case(20, manifest)
end

main()
