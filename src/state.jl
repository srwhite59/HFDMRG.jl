struct StateSelection
    completed::Int
    active::Int
    empty::Int
    discarded_squared_weight::Float64
    limiting_reason::Symbol
    function StateSelection(completed::Int, active::Int, empty::Int,
            discarded_squared_weight::Float64, limiting_reason::Symbol)
        min(completed, active, empty) >= 0 || throw(ArgumentError(
            "state counts must be nonnegative"))
        isfinite(discarded_squared_weight) && discarded_squared_weight >= 0 ||
            throw(ArgumentError("discarded weight must be finite and nonnegative"))
        new(completed, active, empty, discarded_squared_weight,
            limiting_reason)
    end
end

struct RHFStateLink
    old_rank::Int
    physical_width::Int
    active_rank::Int
    close_map::Matrix{Float64}
    completed_map::Matrix{Float64}
    selection::StateSelection
end

function RHFStateLink(old_rank::Int, physical_width::Int,
        close_map::AbstractMatrix{Float64}, selection::StateSelection)
    rows = old_rank + physical_width
    completed = zeros(rows, selection.completed)
    RHFStateLink(old_rank, physical_width, close_map, completed, selection)
end

function RHFStateLink(old_rank::Int, physical_width::Int,
        close_map::AbstractMatrix{Float64},
        completed_map::AbstractMatrix{Float64}, selection::StateSelection)
    min(old_rank, physical_width) >= 0 || throw(ArgumentError(
        "link dimensions must be nonnegative"))
    rows = old_rank + physical_width
    size(close_map, 1) == rows || throw(DimensionMismatch(
        "close-map row count is inconsistent"))
    active = size(close_map, 2)
    active == selection.active || throw(DimensionMismatch(
        "close-map rank and certificate disagree"))
    size(completed_map) == (rows, selection.completed) ||
        throw(DimensionMismatch("completed-map rank and certificate disagree"))
    selection.completed + selection.active + selection.empty == rows ||
        throw(DimensionMismatch("link sectors do not span the parent"))
    all(isfinite, close_map) && all(isfinite, completed_map) ||
        throw(ArgumentError("link maps must be finite"))
    retained = hcat(completed_map, close_map)
    gram = transpose(retained) * retained
    tolerance = 2048eps(Float64) * max(rows, active, 1)
    norm(gram - I, Inf) <= tolerance || throw(ArgumentError(
        "completed and active maps must be mutually isometric"))
    RHFStateLink(old_rank, physical_width, active,
        Matrix{Float64}(close_map), Matrix{Float64}(completed_map), selection)
end

struct RHFRoot
    covariance::Matrix{Float64}
    occupied::Int
end

function RHFRoot(covariance::AbstractMatrix{Float64}; occupied::Int)
    size(covariance, 1) == size(covariance, 2) || throw(DimensionMismatch(
        "root covariance must be square"))
    0 <= occupied <= size(covariance, 1) || throw(DimensionMismatch(
        "root occupation exceeds its active rank"))
    all(isfinite, covariance) || throw(ArgumentError(
        "root covariance must be finite"))
    matrix = Matrix{Float64}(covariance)
    tolerance = 4096eps(Float64) * max(size(matrix, 1), 1)
    norm(matrix - transpose(matrix), Inf) <= tolerance || throw(ArgumentError(
        "root covariance must be symmetric"))
    norm(matrix * matrix - matrix, Inf) <= tolerance || throw(ArgumentError(
        "root covariance must be idempotent"))
    abs(tr(matrix) - occupied) <= tolerance || throw(ArgumentError(
        "root covariance occupation is inconsistent"))
    RHFRoot(matrix, occupied)
end
