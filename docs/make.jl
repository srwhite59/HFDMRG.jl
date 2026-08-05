using Documenter

pushfirst!(LOAD_PATH, joinpath(@__DIR__, ".."))
using HFDMRG

makedocs(
    modules = [HFDMRG],
    sitename = "HFDMRG.jl",
    source = "src",
    build = "build",
    doctest = true,
    checkdocs = :exports,
    warnonly = false,
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = false,
        edit_link = nothing,
        repolink = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Overview" => "manual/index.md",
            "Getting started" => "manual/getting_started.md",
            "Solver controls" => "manual/solver_controls.md",
            "Numerical conventions" => "manual/numerical_conventions.md",
        ],
        "Algorithms" => [
            "Overview" => "algorithms/index.md",
            "Moving-window sweep" => "algorithms/sweep.md",
            "Local RHF and UHF" => "algorithms/local_scf.md",
            "Interaction backends" => "algorithms/interaction_backends.md",
            "History acceleration" => "algorithms/history_acceleration.md",
        ],
        "Examples" => "examples/index.md",
        "Reference" => [
            "Overview" => "reference/index.md",
            "Solvers and controls" => "reference/solvers.md",
            "Layouts and backends" => "reference/layouts_backends.md",
        ],
        "Developer Notes" => [
            "Overview" => "developer/index.md",
            "Backend architecture" => "developer/backend_architecture.md",
        ],
    ],
)
