include(joinpath(@__DIR__, "..", "examples", "hydrogen_chain_fixture.jl"))
using .HydrogenChainFixtures: DEFAULT_MANIFEST, load_hydrogen_chain_fixture
using TOML

@testset "public physical hydrogen fixture schema" begin
    for name in (:h10, :h20)
        fixture = load_hydrogen_chain_fixture(name)
        atoms = name === :h10 ? 10 : 20
        sites = 18atoms + 3
        @test fixture.interaction.sites == sites
        @test size(fixture.one_body.bands, 2) == sites
        @test size(fixture.producer.one_body) == (sites, sites)
        @test size(fixture.producer.interaction) == (sites, sites)
        @test length(fixture.atom_intervals) == atoms
        @test first(first(fixture.atom_intervals)) == 1
        @test last(last(fixture.atom_intervals)) == sites
        @test fixture.producer.constant ==
            fixture.metadata["nuclear_repulsion_ha"]
        @test fixture.provenance["energy_unit"] == "hartree"
        @test fixture.provenance["length_unit"] == "bohr"
    end
    @test_throws ArgumentError load_hydrogen_chain_fixture(:h12)
end

function write_fixture_copy(root, manifest; mutate=nothing)
    record = TOML.parsefile(manifest)
    for name in ("h10", "h20")
        filename = record["cases"][name]["data_file"]
        cp(joinpath(dirname(manifest), filename), joinpath(root, filename))
    end
    mutate === nothing || mutate(record)
    path = joinpath(root, "manifest.toml")
    open(path, "w") do io
        TOML.print(io, record; sorted=true)
    end
    path
end

@testset "public physical hydrogen fixture hostile input" begin
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST)
        open(joinpath(root, "h10.f64le"), "a") do io
            write(io, UInt8(0))
        end
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST;
            mutate=record -> record["cases"]["h10"]["data_bytes"] += 8)
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST;
            mutate=record -> record["cases"]["h10"]["arrays"][2][
                "offset_bytes"] += 8)
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST;
            mutate=record -> record["cases"]["h10"]["arrays"][1][
                "element_count"] += 1)
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST;
            mutate=record -> record["cases"]["h10"]["atom_starts"][2] += 1)
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
    mktempdir() do root
        path = write_fixture_copy(root, DEFAULT_MANIFEST;
            mutate=record -> record["cases"]["h10"]["data_file"] =
                "../h10.f64le")
        @test_throws ArgumentError load_hydrogen_chain_fixture(:h10;
            manifest_path=path)
    end
end
