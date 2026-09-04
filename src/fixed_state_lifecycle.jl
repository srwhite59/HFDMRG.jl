struct RHFStateControl
    maximum_active::Int
    cutoff::Float64
    exact::Bool
    function RHFStateControl(maximum_active::Int, cutoff::Float64;
            exact::Bool=false)
        maximum_active >= 0 || throw(ArgumentError(
            "state maximum rank must be nonnegative"))
        isfinite(cutoff) && cutoff >= 0.0 || throw(ArgumentError(
            "state cutoff must be finite and nonnegative"))
        new(maximum_active, cutoff, exact)
    end
end

mutable struct RHFMovingRoot
    covariance::Matrix{Float64}
    rank::Int
    occupied::Int
    total_occupied::Int
    closes::Int
    opens::Int
    turnarounds::Int
end

struct RHFRootWorkspace
    first::Matrix{Float64}
    second::Matrix{Float64}
end

mutable struct RHFFixedLifecycle
    collection::Vector{RHFOuterBlock}
    root::RHFMovingRoot
    intervals::Vector{UnitRange{Int}}
    atom_intervals::Vector{UnitRange{Int}}
    cells::Int
    forward::Bool
    next_center::Int
    provenance::UInt
    generation::Int
    half_sweeps::Int
    fragment_ingestions::Int
    original_rows_read::Int
    replacements::Int
    center_visits::Int
end

function _traversal_intervals(operator::UnitCellInteraction)
    sites = operator.sites
    sites >= 2 || throw(ArgumentError("RHF traversal requires at least two rows"))
    first_width = _first_group_count(operator)
    intervals = UnitRange{Int}[1:first_width]
    append!(intervals, [first:min(first + UNIT_CELL_WIDTH - 1, sites)
        for first = first_width + 1:UNIT_CELL_WIDTH:sites])
    intervals
end

function _check_traversal_intervals(intervals, sites::Int)
    owned = Vector{UnitRange{Int}}(intervals)
    length(owned) >= 2 || throw(ArgumentError(
        "RHF traversal requires at least two physical groups"))
    previous = 0
    for interval in owned
        !isempty(interval) && first(interval) == previous + 1 &&
            length(interval) <= UNIT_CELL_WIDTH || throw(ArgumentError(
                "physical groups must be ordered contiguous widths at most 18"))
        previous = last(interval)
    end
    previous == sites || throw(ArgumentError(
        "physical groups must partition the complete one-body range"))
    owned
end

@inline _center_width(lifecycle::RHFFixedLifecycle, cell::Int) =
    length(lifecycle.intervals[cell]) + length(lifecycle.intervals[cell + 1])

mutable struct RHFLifecycleWorkspace
    prepared::FastRHFPreparedCenter
    builder::BlockBuildWorkspace
    root::RHFRootWorkspace
    physical::RHFOuterBlock
    empty::RHFOuterBlock
    maximum_rank::Int
end

function RHFLifecycleWorkspace(one_body::BandedOneBody,
        operator::UnitCellInteraction, maximum_active::Int;
        maximum_rhs::Int=3)
    maximum_active >= 0 || throw(ArgumentError(
        "state capacity must be nonnegative"))
    maximum_rank = 2UNIT_CELL_WIDTH + 2maximum_active
    channels = _maximum_channel_rank(operator)
    prepared = FastRHFPreparedCenter(maximum_rank, operator;
        maximum_rhs=max(maximum_rhs, channels))
    builder = BlockBuildWorkspace(max(UNIT_CELL_WIDTH + maximum_active,
        UNIT_CELL_WIDTH), channels)
    root = RHFRootWorkspace(Matrix{Float64}(undef, maximum_rank, maximum_rank),
        Matrix{Float64}(undef, maximum_rank, maximum_rank))
    physical = seed_cell_block(one_body, operator, 1)
    empty = _empty_outer_block(operator, UInt(0))
    RHFLifecycleWorkspace(prepared, builder, root, physical, empty,
        maximum_rank)
end

function _physical_block(workspace::RHFLifecycleWorkspace, cell::Int,
        rows::UnitRange{Int}, one_body::BandedOneBody,
        operator::UnitCellInteraction, provenance::UInt;
        reflected::Bool=false)
    if length(rows) == UNIT_CELL_WIDTH
        _prepare_physical_block!(workspace.physical, cell, rows, one_body,
            operator, provenance; reflected)
    else
        seed_interval_block(one_body, operator, cell, rows; provenance,
            reflected)
    end
end

function _check_moving_root(root::RHFMovingRoot)
    capacity = size(root.covariance, 1)
    size(root.covariance, 2) == capacity && 0 <= root.rank <= capacity ||
        throw(DimensionMismatch("moving-root capacity is inconsistent"))
    0 <= root.occupied <= root.rank && root.occupied <= root.total_occupied ||
        throw(DimensionMismatch("moving-root occupation is inconsistent"))
    root.closes >= 0 && root.opens >= 0 && root.turnarounds >= 0 ||
        throw(ArgumentError("moving-root counters are invalid"))
    matrix = @view root.covariance[1:root.rank, 1:root.rank]
    all(isfinite, matrix) || throw(ArgumentError(
        "moving-root covariance must be finite"))
    tolerance = 8192eps(Float64) * max(root.rank, 1)
    trace = 0.0
    maximum_symmetry = 0.0
    maximum_idempotency = 0.0
    @inbounds for column = 1:root.rank, row = 1:root.rank
        maximum_symmetry = max(maximum_symmetry,
            abs(matrix[row, column] - matrix[column, row]))
        value = -matrix[row, column]
        for inner = 1:root.rank
            value = muladd(matrix[row, inner], matrix[inner, column], value)
        end
        maximum_idempotency = max(maximum_idempotency, abs(value))
        row == column && (trace += matrix[row, column])
    end
    maximum_symmetry <= tolerance && maximum_idempotency <= tolerance &&
        abs(trace - root.occupied) <= tolerance || throw(ArgumentError(
            "moving-root covariance is not an occupied projector"))
    nothing
end

function _sign_columns!(matrix)
    @inbounds for column in axes(matrix, 2)
        isempty(axes(matrix, 1)) && continue
        pivot = first(axes(matrix, 1))
        for row in axes(matrix, 1)
            abs(matrix[row, column]) > abs(matrix[pivot, column]) &&
                (pivot = row)
        end
        matrix[pivot, column] < 0.0 &&
            (@views matrix[:, column] .*= -1.0)
    end
    matrix
end

@inline function _state_cluster_envelope(rows::Int, count::Int,
        scale::Float64)
    64eps(Float64) * max(1, rows, count) * max(1.0, scale)
end

function _state_spectral_clusters(values, envelope::Float64)
    clusters = UnitRange{Int}[]
    first_index = 1
    @inbounds for index = 1:length(values)-1
        if values[index] - values[index+1] > envelope
            push!(clusters, first_index:index)
            first_index = index + 1
        end
    end
    isempty(values) || push!(clusters, first_index:length(values))
    clusters
end

function _finite_rhf_selection(eigenvalues, rows::Int, maximum_occupied::Int,
        control::RHFStateControl)
    explicit = min(rows, maximum_occupied)
    order = collect(length(eigenvalues):-1:length(eigenvalues)-explicit+1)
    singular = [sqrt(clamp(eigenvalues[index], 0.0, 1.0)) for index in order]
    left_weights = [clamp(eigenvalues[index], 0.0, 1.0) for index in order]
    right_weights = 1.0 .- left_weights
    clusters = _state_spectral_clusters(singular,
        _state_cluster_envelope(rows, explicit,
            isempty(singular) ? 1.0 : singular[1]))
    boundaries = Int[0]
    append!(boundaries, last.(clusters))
    best_meeting = nothing
    best_capped = nothing
    for completed in boundaries, kept_end in boundaries
        kept_end < completed && continue
        active = kept_end - completed
        active <= control.maximum_active || continue
        fidelity = 1.0
        @inbounds for index = 1:completed
            fidelity *= max(0.0, 1.0 - right_weights[index])
        end
        @inbounds for index = kept_end+1:explicit
            fidelity *= max(0.0, 1.0 - left_weights[index])
        end
        discarded = max(0.0, 1.0 - fidelity^2)
        candidate = (; completed, active, kept_end, discarded)
        if discarded <= control.cutoff && (best_meeting === nothing ||
                active < best_meeting.active ||
                (active == best_meeting.active &&
                    discarded < best_meeting.discarded))
            best_meeting = candidate
        end
        if best_capped === nothing || active > best_capped.active ||
                (active == best_capped.active &&
                    discarded < best_capped.discarded)
            best_capped = candidate
        end
    end
    best_capped === nothing && throw(DimensionMismatch(
        "no whole-cluster state selection fits the hard maximum"))
    selected = best_meeting === nothing ? best_capped : best_meeting
    reason = best_meeting === nothing ? :maximum_rank : :cutoff
    completed = order[1:selected.completed]
    active = order[selected.completed+1:selected.kept_end]
    completed, active, selected.discarded, reason
end

function _factor_state_link(local_density, left_rank::Int, right_rank::Int,
        control::RHFStateControl; cross=nothing,
        maximum_occupied::Union{Nothing,Int}=nothing)
    rows = left_rank + right_rank
    size(local_density) == (rows, rows) || throw(DimensionMismatch(
        "local restriction has the wrong dimensions"))
    factors = eigen(Symmetric(Matrix(local_density)))
    tolerance = 8192eps(Float64) * max(rows, 1)
    all(value -> isfinite(value) && -tolerance <= value <= 1.0 + tolerance,
        factors.values) || throw(ArgumentError(
        "local restriction lies outside projector bounds"))
    completed = Int[]
    active = Int[]
    undecided = Int[]
    exact_scores = zeros(length(factors.values))
    @inbounds for index in eachindex(factors.values)
        value = factors.values[index]
        if control.exact
            leakage = 0.0
            if cross !== nothing
                size(cross, 1) == rows || throw(DimensionMismatch(
                    "local leakage certificate has the wrong row count"))
                for outside in axes(cross, 2)
                    component = 0.0
                    for row = 1:rows
                        component = muladd(factors.vectors[row, index],
                            cross[row, outside], component)
                    end
                    leakage = max(leakage, abs(component))
                end
            end
            if leakage > 0.0 || (cross === nothing && 0.0 < value < 1.0)
                push!(active, index)
                # A Float64 eigensolve can give structurally null restriction
                # directions nonzero cross components.  The covariance weight
                # (equivalently the squared Schmidt weight), not the raw cross
                # amplitude, orders the certified occupied range.
                exact_scores[index] = max(value, leakage * leakage)
            elseif value >= 0.5
                push!(completed, index)
            end
        elseif value <= tolerance
            continue
        elseif 1.0 - value <= tolerance
            push!(completed, index)
        else
            push!(undecided, index)
        end
    end
    reason = control.exact ? :exact : :cutoff
    discarded = 0.0
    if control.exact && maximum_occupied !== nothing
        maximum_occupied >= 0 || throw(ArgumentError(
            "exact occupied-rank certificate must be nonnegative"))
        available = maximum_occupied - length(completed)
        available >= 0 || throw(DimensionMismatch(
            "completed restriction exceeds occupied-rank certificate"))
        if length(active) > available
            order = sortperm(active; by=index -> exact_scores[index], rev=true)
            resize!(order, available)
            active = active[order]
        end
    end
    if !control.exact
        nominal = maximum_occupied === nothing ? rows : maximum_occupied
        completed, active, discarded, reason = _finite_rhf_selection(
            factors.values, rows, nominal, control)
    end
    sort!(completed); sort!(active)
    completed_map = Matrix(factors.vectors[:, completed])
    active_map = Matrix(factors.vectors[:, active])
    _sign_columns!(completed_map); _sign_columns!(active_map)
    selection = StateSelection(length(completed), length(active),
        rows - length(completed) - length(active), discarded, reason)
    RHFStateLink(left_rank, right_rank, active_map, completed_map, selection)
end

function _canonicalize_candidate!(candidate, rank::Int, occupied::Int)
    rank == 0 && (occupied == 0 || throw(DimensionMismatch(
        "empty root cannot retain occupation")); return nothing)
    factors = eigen(Symmetric(Matrix(@view(candidate[1:rank, 1:rank]))))
    tolerance = 16384eps(Float64) * max(rank, 1)
    all(value -> isfinite(value) && -tolerance <= value <= 1.0 + tolerance,
        factors.values) || throw(ArgumentError(
        "candidate root lies outside covariance bounds"))
    0 <= occupied <= rank || throw(DimensionMismatch(
        "candidate root occupation exceeds rank"))
    fill!(@view(candidate[1:rank, 1:rank]), 0.0)
    if occupied > 0
        vectors = @view factors.vectors[:, rank-occupied+1:rank]
        mul!(@view(candidate[1:rank, 1:rank]), vectors, transpose(vectors))
    end
    nothing
end

function _close_root!(root::RHFMovingRoot, link::RHFStateLink,
        local_rows::Int, local_first::Bool, workspace::RHFRootWorkspace;
        exact::Bool)
    _check_moving_root(root)
    link.old_rank + link.physical_width == local_rows ||
        throw(DimensionMismatch("closing link has the wrong local rank"))
    local_rows <= root.rank || throw(DimensionMismatch(
        "closing local rank exceeds root rank"))
    future = root.rank - local_rows
    active = link.active_rank
    new_rank = active + future
    new_occupied = root.occupied - link.selection.completed
    0 <= new_occupied <= new_rank || throw(DimensionMismatch(
        "closing certificate loses an occupied direction"))
    output = workspace.first
    fill!(@view(output[1:new_rank, 1:new_rank]), 0.0)
    density = root.covariance
    map = link.close_map
    if local_first
        @inbounds for a = 1:active, b = 1:active, i = 1:local_rows,
                j = 1:local_rows
            output[a, b] = muladd(map[i, a] * density[i, j], map[j, b],
                output[a, b])
        end
        @inbounds for a = 1:active, f = 1:future, i = 1:local_rows
            output[a, active+f] = muladd(map[i, a],
                density[i, local_rows+f], output[a, active+f])
        end
        @views output[active+1:new_rank, active+1:new_rank] .=
            density[local_rows+1:root.rank, local_rows+1:root.rank]
    else
        @views output[1:future, 1:future] .= density[1:future, 1:future]
        @inbounds for f = 1:future, a = 1:active, i = 1:local_rows
            output[f, future+a] = muladd(density[f, future+i], map[i, a],
                output[f, future+a])
        end
        @inbounds for a = 1:active, b = 1:active, i = 1:local_rows,
                j = 1:local_rows
            output[future+a, future+b] = muladd(map[i, a] *
                density[future+i, future+j], map[j, b],
                output[future+a, future+b])
        end
    end
    @inbounds for column = 1:new_rank, row = 1:column
        local_first && row <= active && column > active &&
            (output[column, row] = output[row, column])
        !local_first && row <= future && column > future &&
            (output[column, row] = output[row, column])
    end
    exact || _canonicalize_candidate!(output, new_rank, new_occupied)
    copyto!(@view(root.covariance[1:new_rank, 1:new_rank]),
        @view(output[1:new_rank, 1:new_rank]))
    root.rank = new_rank; root.occupied = new_occupied; root.closes += 1
    _check_moving_root(root)
    root
end

function _open_root!(root::RHFMovingRoot, link::RHFStateLink,
        local_first::Bool, workspace::RHFRootWorkspace)
    _check_moving_root(root)
    active = link.active_rank
    active <= root.rank || throw(DimensionMismatch(
        "opening active rank exceeds root rank"))
    future = root.rank - active
    parent = link.old_rank + link.physical_width
    new_rank = future + parent
    new_occupied = root.occupied + link.selection.completed
    new_rank <= size(root.covariance, 1) || throw(DimensionMismatch(
        "opening exceeds root capacity"))
    output = workspace.first
    fill!(@view(output[1:new_rank, 1:new_rank]), 0.0)
    density = root.covariance
    map = link.close_map
    completed = link.completed_map
    if local_first
        @inbounds for i = 1:parent, j = 1:parent, a = 1:active,
                b = 1:active
            output[i, j] = muladd(map[i, a] * density[a, b], map[j, b],
                output[i, j])
        end
        @inbounds for i = 1:parent, j = 1:parent, mode in axes(completed, 2)
            output[i, j] = muladd(completed[i, mode], completed[j, mode],
                output[i, j])
        end
        @inbounds for i = 1:parent, f = 1:future, a = 1:active
            output[i, parent+f] = muladd(map[i, a], density[a, active+f],
                output[i, parent+f])
        end
        @views output[parent+1:new_rank, parent+1:new_rank] .=
            density[active+1:root.rank, active+1:root.rank]
        @inbounds for i = 1:parent, f = 1:future
            output[parent+f, i] = output[i, parent+f]
        end
    else
        @views output[1:future, 1:future] .= density[1:future, 1:future]
        @inbounds for f = 1:future, i = 1:parent, a = 1:active
            output[f, future+i] = muladd(density[f, future+a], map[i, a],
                output[f, future+i])
        end
        @inbounds for i = 1:parent, j = 1:parent, a = 1:active,
                b = 1:active
            output[future+i, future+j] = muladd(map[i, a] *
                density[future+a, future+b], map[j, b],
                output[future+i, future+j])
        end
        @inbounds for i = 1:parent, j = 1:parent, mode in axes(completed, 2)
            output[future+i, future+j] = muladd(completed[i, mode],
                completed[j, mode], output[future+i, future+j])
        end
        @inbounds for f = 1:future, i = 1:parent
            output[future+i, f] = output[f, future+i]
        end
    end
    copyto!(@view(root.covariance[1:new_rank, 1:new_rank]),
        @view(output[1:new_rank, 1:new_rank]))
    root.rank = new_rank; root.occupied = new_occupied; root.opens += 1
    _check_moving_root(root)
    root
end

function _dimer_orbital(one_body::BandedOneBody, rows::UnitRange{Int})
    dimension = length(rows)
    matrix = Matrix{Float64}(undef, dimension, dimension)
    @inbounds for column = 1:dimension, row = 1:dimension
        matrix[row, column] = h1_entry(one_body, rows[row], rows[column])
    end
    factors = eigen(Symmetric(matrix))
    orbital = Vector(factors.vectors[:, 1])
    _sign_columns!(reshape(orbital, :, 1))
    orbital
end

_dimer_orbital(one_body::BandedOneBody, first_cell::Int) = _dimer_orbital(
    one_body, (first_cell - 1) * UNIT_CELL_WIDTH + 1:
        (first_cell + 1) * UNIT_CELL_WIDTH)

struct _RHFLocalFragment
    rows::UnitRange{Int}
    orbital::Vector{Float64}
end

function _rhf_dimer_fragments(one_body::BandedOneBody, atom_intervals)
    atoms = Vector{UnitRange{Int}}(atom_intervals)
    iseven(length(atoms)) && !isempty(atoms) || throw(ArgumentError(
        "RHF dimer initialization requires an even atom count"))
    previous = 0
    for interval in atoms
        !isempty(interval) && first(interval) == previous + 1 ||
            throw(ArgumentError("atom intervals must be ordered and contiguous"))
        previous = last(interval)
    end
    previous == h1_sites(one_body) || throw(ArgumentError(
        "atom intervals must partition the complete one-body range"))
    [_RHFLocalFragment(first(atoms[index]):last(atoms[index + 1]),
        _dimer_orbital(one_body,
            first(atoms[index]):last(atoms[index + 1])))
        for index = 1:2:length(atoms)]
end

mutable struct _RHFFragmentCursor
    fragments::Vector{_RHFLocalFragment}
    active_ids::Vector{Int}
    active_coordinates::Matrix{Float64}
    completed::BitVector
    ingested::BitVector
    ingestions::Int
end

function _RHFFragmentCursor(fragments)
    owned = Vector{_RHFLocalFragment}(fragments)
    count = length(owned)
    _RHFFragmentCursor(owned, Int[], zeros(0, 0), falses(count),
        falses(count), 0)
end

@inline function _fragment_intersection(first::UnitRange{Int},
        second::UnitRange{Int})
    lo, hi = max(Base.first(first), Base.first(second)),
        min(last(first), last(second))
    lo <= hi ? (lo:hi) : (1:0)
end

function _fragment_parent(cursor::_RHFFragmentCursor, rows::UnitRange{Int},
        old_rank::Int, physical_first::Bool)
    ids = copy(cursor.active_ids)
    for id in eachindex(cursor.fragments)
        cursor.completed[id] && continue
        overlap = _fragment_intersection(cursor.fragments[id].rows, rows)
        isempty(overlap) && continue
        id in ids || push!(ids, id)
        if !cursor.ingested[id]
            cursor.ingested[id] = true
            cursor.ingestions += 1
        end
    end
    width = length(rows)
    parent = zeros(width + old_rank, length(ids))
    cross = zeros(width + old_rank, length(ids))
    old_offset = physical_first ? width : 0
    physical_offset = physical_first ? 0 : old_rank
    for (column, id) in enumerate(ids)
        old_column = findfirst(==(id), cursor.active_ids)
        old_column === nothing || copyto!(
            @view(parent[old_offset+1:old_offset+old_rank, column]),
            @view(cursor.active_coordinates[:, old_column]))
        fragment = cursor.fragments[id]
        overlap = _fragment_intersection(fragment.rows, rows)
        @inbounds for physical in overlap
            parent[physical_offset + physical - first(rows) + 1, column] =
                fragment.orbital[physical - first(fragment.rows) + 1]
        end
        future = physical_first ?
            (first(fragment.rows):min(last(fragment.rows), first(rows) - 1)) :
            (max(first(fragment.rows), last(rows) + 1):last(fragment.rows))
        future_weight = isempty(future) ? 0.0 : sum(abs2,
            @view(fragment.orbital[future .- first(fragment.rows) .+ 1]))
        future_weight > 0.0 && copyto!(@view(cross[:, column]),
            sqrt(future_weight) .* @view(parent[:, column]))
    end
    ids, parent, cross
end

function _advance_fragment_cursor!(cursor::_RHFFragmentCursor, ids, parent,
        link::RHFStateLink)
    active = transpose(link.close_map) * parent
    completed = transpose(link.completed_map) * parent
    tolerance = 32768eps(Float64) * max(size(parent)..., 1)
    next_ids = Int[]
    columns = Int[]
    for (column, id) in enumerate(ids)
        active_weight = sum(abs2, @view(active[:, column]))
        completed_weight = sum(abs2, @view(completed[:, column]))
        if active_weight > tolerance^2
            push!(next_ids, id); push!(columns, column)
        elseif completed_weight > tolerance^2
            cursor.completed[id] = true
        end
    end
    cursor.active_ids = next_ids
    cursor.active_coordinates = isempty(columns) ? zeros(link.active_rank, 0) :
        Matrix(@view(active[:, columns]))
    cursor
end

function _root_fragment_projector(cursor::_RHFFragmentCursor,
        rows::UnitRange{Int}, old_rank::Int, physical_first::Bool,
        total_occupied::Int)
    ids, parent, _ = _fragment_parent(cursor, rows, old_rank, physical_first)
    covariance = parent * transpose(parent)
    factors = eigen(Symmetric(covariance))
    tolerance = 32768eps(Float64) * max(size(parent)..., 1)
    realized = count(>(tolerance), factors.values)
    expected = total_occupied - count(cursor.completed)
    realized == expected || throw(ArgumentError(
        "dimer initializer lost a live occupied direction"))
    basis = @view factors.vectors[:, end-realized+1:end]
    projector = Matrix{Float64}(undef, size(parent, 1), size(parent, 1))
    mul!(projector, basis, transpose(basis))
    all(cursor.ingested) || throw(ArgumentError(
        "dimer initializer omitted a local occupied fragment"))
    projector, realized
end

function initialize_rhf_dimer_lifecycle(one_body::BandedOneBody,
        operator::UnitCellInteraction, control::RHFStateControl,
        workspace::RHFLifecycleWorkspace; forward::Bool=true,
        atom_intervals=nothing)
    h1_sites(one_body) == operator.sites || throw(DimensionMismatch(
        "H1 and interaction lengths disagree"))
    intervals = _check_traversal_intervals(_traversal_intervals(operator),
        operator.sites)
    atoms = atom_intervals === nothing ? intervals :
        Vector{UnitRange{Int}}(atom_intervals)
    fragments = _rhf_dimer_fragments(one_body, atoms)
    cells = length(intervals)
    provenance = hash((operator.provenance, intervals, atoms,
        control.maximum_active,
        control.cutoff, control.exact))
    empty = _empty_outer_block(operator, provenance)
    collection = Vector{RHFOuterBlock}(undef, cells - 1)
    cursor = _RHFFragmentCursor(fragments)
    current = empty
    if forward
        collection[1] = empty
        for cell = cells:-1:3
            physical = _physical_block(workspace, cell, intervals[cell],
                one_body, operator, provenance)
            ids, parent, cross = _fragment_parent(cursor, intervals[cell],
                current.rank, true)
            density = parent * transpose(parent)
            link = _factor_state_link(density, length(intervals[cell]),
                current.rank, control; cross,
                maximum_occupied=length(ids))
            current = build_outer_block(physical, current, link, one_body,
                operator, workspace.builder)
            collection[cell-1] = current
            _advance_fragment_cursor!(cursor, ids, parent, link)
        end
    else
        collection[cells-1] = empty
        for cell = 1:cells-2
            physical = _physical_block(workspace, cell, intervals[cell],
                one_body, operator, provenance; reflected=true)
            ids, parent, cross = _fragment_parent(cursor, intervals[cell],
                current.rank, false)
            density = parent * transpose(parent)
            link = _factor_state_link(density, current.rank,
                length(intervals[cell]), control; cross,
                maximum_occupied=length(ids))
            current = build_outer_block(current, physical, link, one_body,
                operator, workspace.builder)
            current.reflected = true
            collection[cell] = current
            _advance_fragment_cursor!(cursor, ids, parent, link)
        end
    end
    center_first = first(intervals[forward ? 1 : cells - 1])
    center_last = last(intervals[forward ? 2 : cells])
    center_rows = center_first:center_last
    projector, occupied = _root_fragment_projector(cursor, center_rows,
        current.rank, forward, length(fragments))
    capacity = workspace.maximum_rank
    covariance = Matrix{Float64}(undef, capacity, capacity)
    fill!(covariance, NaN)
    rank = size(projector, 1)
    copyto!(@view(covariance[1:rank, 1:rank]), projector)
    root = RHFMovingRoot(covariance, rank, occupied, length(fragments),
        0, 0, 0)
    lifecycle = RHFFixedLifecycle(collection, root, intervals, atoms, cells,
        forward,
        forward ? 1 : cells-1,
        provenance, 0, 0, cursor.ingestions, 0, 0, 0)
    _check_moving_root(root)
    lifecycle
end

function _lifecycle_blocks(lifecycle::RHFFixedLifecycle,
        workspace::RHFLifecycleWorkspace)
    cell = lifecycle.next_center
    workspace.empty.provenance = lifecycle.provenance
    left = cell == 1 ? workspace.empty : lifecycle.collection[cell-1]
    right = cell == lifecycle.cells-1 ? workspace.empty :
        lifecycle.collection[cell+1]
    expected = left.rank + _center_width(lifecycle, cell) + right.rank
    lifecycle.root.rank == expected || throw(DimensionMismatch(
        "collection and live-root ranks disagree at center $cell"))
    lifecycle.root.occupied + left.completed_count + right.completed_count ==
        lifecycle.root.total_occupied || throw(ArgumentError(
        "collection and live-root occupied counts disagree"))
    left, right
end

function _restore_root!(root::RHFMovingRoot, backup, rank, occupied,
        closes, opens, turnarounds)
    copyto!(@view(root.covariance[1:rank, 1:rank]),
        @view(backup[1:rank, 1:rank]))
    root.rank = rank; root.occupied = occupied
    root.closes = closes; root.opens = opens; root.turnarounds = turnarounds
    root
end

function _build_advance_candidate!(root::RHFMovingRoot,
        lifecycle::RHFFixedLifecycle, left::RHFOuterBlock,
        right::RHFOuterBlock, one_body::BandedOneBody,
        operator::UnitCellInteraction, control::RHFStateControl,
        workspace::RHFLifecycleWorkspace, root_workspace::RHFRootWorkspace)
    cell = lifecycle.next_center
    forward = lifecycle.forward
    terminal = forward ? cell == lifecycle.cells-1 : cell == 1
    if forward && !terminal
        width = length(lifecycle.intervals[cell])
        local_rows = left.rank + width
        density = @view root.covariance[1:local_rows, 1:local_rows]
        cross = @view root.covariance[1:local_rows,
            local_rows+1:root.rank]
        link = _factor_state_link(density, left.rank, width,
            control; cross=cross, maximum_occupied=root.occupied)
        physical = _physical_block(workspace, cell,
            lifecycle.intervals[cell], one_body, operator,
            lifecycle.provenance)
        candidate = build_outer_block(left, physical, link,
            one_body, operator, workspace.builder)
        _close_root!(root, link, local_rows, true, root_workspace;
            exact=control.exact)
        _open_root!(root, right.link, false, root_workspace)
        next_forward, next_center = true, cell + 1
    elseif forward
        width = length(lifecycle.intervals[cell + 1])
        local_rows = width + right.rank
        start = root.rank - local_rows + 1
        density = @view root.covariance[start:root.rank, start:root.rank]
        cross = transpose(@view root.covariance[1:start-1, start:root.rank])
        link = _factor_state_link(density, width, right.rank,
            control; cross=cross, maximum_occupied=root.occupied)
        physical = _physical_block(workspace, cell + 1,
            lifecycle.intervals[cell + 1], one_body, operator,
            lifecycle.provenance; reflected=true)
        candidate = build_outer_block(physical, right, link,
            one_body, operator, workspace.builder)
        candidate.reflected = true
        _close_root!(root, link, local_rows, false, root_workspace;
            exact=control.exact)
        lifecycle.cells == 2 ?
            _open_root!(root, candidate.link, false, root_workspace) :
            _open_root!(root, left.link, true, root_workspace)
        next_forward = false
        next_center = lifecycle.cells == 2 ? 1 : cell - 1
        root.turnarounds += 1
    elseif !terminal
        width = length(lifecycle.intervals[cell + 1])
        local_rows = width + right.rank
        start = root.rank - local_rows + 1
        density = @view root.covariance[start:root.rank, start:root.rank]
        cross = transpose(@view root.covariance[1:start-1, start:root.rank])
        link = _factor_state_link(density, width, right.rank,
            control; cross=cross, maximum_occupied=root.occupied)
        physical = _physical_block(workspace, cell + 1,
            lifecycle.intervals[cell + 1], one_body, operator,
            lifecycle.provenance; reflected=true)
        candidate = build_outer_block(physical, right, link,
            one_body, operator, workspace.builder)
        candidate.reflected = true
        _close_root!(root, link, local_rows, false, root_workspace;
            exact=control.exact)
        _open_root!(root, left.link, true, root_workspace)
        next_forward, next_center = false, cell - 1
    else
        width = length(lifecycle.intervals[cell])
        local_rows = left.rank + width
        density = @view root.covariance[1:local_rows, 1:local_rows]
        cross = @view root.covariance[1:local_rows,
            local_rows+1:root.rank]
        link = _factor_state_link(density, left.rank, width,
            control; cross=cross, maximum_occupied=root.occupied)
        physical = _physical_block(workspace, cell,
            lifecycle.intervals[cell], one_body, operator,
            lifecycle.provenance)
        candidate = build_outer_block(left, physical, link,
            one_body, operator, workspace.builder)
        _close_root!(root, link, local_rows, true, root_workspace;
            exact=control.exact)
        lifecycle.cells == 2 ?
            _open_root!(root, candidate.link, true, root_workspace) :
            _open_root!(root, right.link, false, root_workspace)
        next_forward = true
        next_center = lifecycle.cells == 2 ? 1 : cell + 1
        root.turnarounds += 1
    end
    (block=candidate, root=root, cell=cell, forward=next_forward,
        next_center=next_center, turned=next_forward != forward,
        selection=candidate.link.selection,
        published_bytes=block_storage_bytes(candidate))
end

function _publish_advance_candidate!(lifecycle::RHFFixedLifecycle,
        candidate)
    root = lifecycle.root
    trial = candidate.root
    trial === root || begin
        copyto!(@view(root.covariance[1:trial.rank, 1:trial.rank]),
            @view(trial.covariance[1:trial.rank, 1:trial.rank]))
        root.rank = trial.rank; root.occupied = trial.occupied
        root.total_occupied = trial.total_occupied
        root.closes = trial.closes; root.opens = trial.opens
        root.turnarounds = trial.turnarounds
    end
    lifecycle.collection[candidate.cell] = candidate.block
    lifecycle.forward = candidate.forward
    lifecycle.next_center = candidate.next_center
    candidate.turned && (lifecycle.generation += 1)
    lifecycle.replacements += 1
    (selection=candidate.selection, published_bytes=candidate.published_bytes)
end


function _transactional_advance!(lifecycle::RHFFixedLifecycle,
        left::RHFOuterBlock, right::RHFOuterBlock,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::RHFStateControl, workspace::RHFLifecycleWorkspace)
    root = lifecycle.root
    old_rank, old_occupied = root.rank, root.occupied
    old_closes, old_opens, old_turnarounds = root.closes, root.opens,
        root.turnarounds
    copyto!(@view(workspace.root.second[1:old_rank, 1:old_rank]),
        @view(root.covariance[1:old_rank, 1:old_rank]))
    candidate = try
        _build_advance_candidate!(root, lifecycle, left, right, one_body,
            operator, control, workspace, workspace.root)
    catch
        _restore_root!(root, workspace.root.second, old_rank, old_occupied,
            old_closes, old_opens, old_turnarounds)
        rethrow()
    end
    _publish_advance_candidate!(lifecycle, candidate)
end

function run_fixed_half_sweep!(lifecycle::RHFFixedLifecycle,
        one_body::BandedOneBody, operator::UnitCellInteraction,
        control::RHFStateControl, workspace::RHFLifecycleWorkspace)
    _check_moving_root(lifecycle.root)
    lifecycle.provenance == hash((operator.provenance, lifecycle.intervals,
        lifecycle.atom_intervals,
        control.maximum_active, control.cutoff, control.exact)) ||
        throw(ArgumentError("lifecycle controls or operator provenance changed"))
    starting_orientation = lifecycle.forward
    evidence = NamedTuple[]
    while lifecycle.forward == starting_orientation
        left, right = _lifecycle_blocks(lifecycle, workspace)
        cell = lifecycle.next_center
        _prepare_lifecycle_center!(workspace.prepared, left, right, cell,
            lifecycle.intervals[cell], lifecycle.intervals[cell + 1],
            one_body, operator)
        energy = _center_energy_fock!(workspace.prepared,
            @view(lifecycle.root.covariance[1:lifecycle.root.rank,
                1:lifecycle.root.rank]), lifecycle.root.occupied)
        lifecycle.center_visits += 1
        publication = _transactional_advance!(lifecycle, left, right, one_body,
            operator,
            control, workspace)
        selection = publication.selection
        push!(evidence, (cell=cell, forward=starting_orientation,
            energy=energy.total, root_rank=workspace.prepared.active_rank,
            left_rank=left.rank, right_rank=right.rank,
            left_completed=left.completed_count,
            right_completed=right.completed_count,
            completed=selection.completed, active=selection.active,
            empty=selection.empty,
            discarded_squared_weight=selection.discarded_squared_weight,
            limiting_reason=selection.limiting_reason,
            published_bytes=publication.published_bytes))
    end
    lifecycle.half_sweeps += 1
    evidence
end

function lifecycle_storage_bytes(lifecycle::RHFFixedLifecycle,
        operator::UnitCellInteraction, workspace::RHFLifecycleWorkspace)
    collection = sum(block_storage_bytes, lifecycle.collection)
    root = sizeof(lifecycle.root.covariance)
    static_operator = sizeof(operator.tail_transfer) +
        sizeof(operator.residual_transfer) + sizeof(operator.phase_source) +
        sizeof(operator.phase_sink) + sizeof(operator.local_triangle)
    reusable = Base.summarysize(workspace)
    (collection=collection, root=root, static_operator=static_operator,
        reusable_workspace=reusable)
end

function _expand_prefix_mode!(output, column::Int, initial, slot::Int,
        lifecycle::RHFFixedLifecycle, workspace::RHFRootWorkspace)
    length(initial) == (slot == 0 ? 0 : lifecycle.collection[slot].rank) ||
        throw(DimensionMismatch("prefix reconstruction rank mismatch"))
    current = initial
    first, second = @view(workspace.first[:, 1]), @view(workspace.second[:, 1])
    for index = slot:-1:1
        block = lifecycle.collection[index]
        parent = block.link.old_rank + block.link.physical_width
        mul!(@view(first[1:parent]), block.link.close_map, current)
        physical = block.link.physical_width
        rows = lifecycle.intervals[index]
        physical == length(rows) || throw(DimensionMismatch(
            "prefix reconstruction physical width is inconsistent"))
        copyto!(@view(output[rows, column]),
            @view(first[block.link.old_rank+1:parent]))
        copyto!(@view(second[1:block.link.old_rank]),
            @view(first[1:block.link.old_rank]))
        current = @view second[1:block.link.old_rank]
    end
    nothing
end

function _expand_suffix_mode!(output, column::Int, initial, slot::Int,
        lifecycle::RHFFixedLifecycle, workspace::RHFRootWorkspace)
    length(initial) == (slot >= lifecycle.cells ? 0 :
        lifecycle.collection[slot].rank) || throw(DimensionMismatch(
        "suffix reconstruction rank mismatch"))
    current = initial
    first, second = @view(workspace.first[:, 1]), @view(workspace.second[:, 1])
    for index = slot:lifecycle.cells-1
        block = lifecycle.collection[index]
        parent = block.link.old_rank + block.link.physical_width
        mul!(@view(first[1:parent]), block.link.close_map, current)
        rows = lifecycle.intervals[index + 1]
        width = length(rows)
        block.link.old_rank == width || throw(DimensionMismatch(
            "suffix reconstruction physical width is inconsistent"))
        copyto!(@view(output[rows, column]), @view(first[1:width]))
        child = block.link.physical_width
        copyto!(@view(second[1:child]),
            @view(first[width+1:parent]))
        current = @view second[1:child]
    end
    nothing
end

function _preflight_terminal_reconstruction(output, lifecycle, workspace)
    size(output) == (last(lifecycle.intervals[end]),
        lifecycle.root.total_occupied) || throw(DimensionMismatch(
        "terminal reconstruction output has the wrong dimensions"))
    arrays = (output, lifecycle.root.covariance, workspace.root.first,
        workspace.root.second)
    @inbounds for first = 1:length(arrays)-1, second = first+1:length(arrays)
        Base.mightalias(arrays[first], arrays[second]) && throw(ArgumentError(
            "terminal reconstruction arrays must not alias"))
    end
    _check_moving_root(lifecycle.root)
    cell = lifecycle.next_center
    1 <= cell < lifecycle.cells || throw(BoundsError(1:lifecycle.cells-1,
        cell))
    left_slot = cell - 1
    right_slot = cell + 1
    previous_rank = 0
    completed = lifecycle.root.occupied
    for slot = 1:left_slot
        block = lifecycle.collection[slot]
        block.cells == (1:slot) || throw(ArgumentError(
            "prefix collection interval is inconsistent"))
        block.link.old_rank == previous_rank &&
            block.link.physical_width == length(lifecycle.intervals[slot]) &&
            block.link.active_rank == block.rank || throw(DimensionMismatch(
            "prefix reconstruction link is inconsistent"))
        previous_rank = block.rank
        completed += block.link.selection.completed
    end
    left_rank = previous_rank
    previous_rank = 0
    for slot = lifecycle.cells-1:-1:right_slot
        block = lifecycle.collection[slot]
        block.cells == (slot+1:lifecycle.cells) || throw(ArgumentError(
            "suffix collection interval is inconsistent"))
        block.link.old_rank == length(lifecycle.intervals[slot + 1]) &&
            block.link.physical_width == previous_rank &&
            block.link.active_rank == block.rank || throw(DimensionMismatch(
            "suffix reconstruction link is inconsistent"))
        previous_rank = block.rank
        completed += block.link.selection.completed
    end
    right_rank = previous_rank
    lifecycle.root.rank == left_rank + _center_width(lifecycle, cell) +
        right_rank ||
        throw(DimensionMismatch("terminal root and collection disagree"))
    completed == lifecycle.root.total_occupied || throw(ArgumentError(
        "terminal reconstruction occupation count is incomplete"))
    maximum_parent = max(lifecycle.root.rank,
        maximum((block.link.old_rank + block.link.physical_width
            for block in lifecycle.collection); init=0))
    maximum_parent <= size(workspace.root.first, 1) &&
        maximum_parent <= size(workspace.root.second, 1) ||
        throw(DimensionMismatch("terminal reconstruction capacity is too small"))
    cell, left_slot, right_slot, left_rank, right_rank
end

function reconstruct_terminal!(output::AbstractMatrix{Float64},
        lifecycle::RHFFixedLifecycle, workspace::RHFLifecycleWorkspace)
    cell, left_slot, right_slot, left_rank, right_rank =
        _preflight_terminal_reconstruction(output, lifecycle, workspace)
    factors = eigen(Symmetric(Matrix(@view(lifecycle.root.covariance[
        1:lifecycle.root.rank, 1:lifecycle.root.rank]))))
    fill!(output, 0.0)
    first_root = lifecycle.root.rank - lifecycle.root.occupied + 1
    column = 0
    for mode = first_root:lifecycle.root.rank
        column += 1
        vector = @view factors.vectors[:, mode]
        left_rank > 0 && _expand_prefix_mode!(output, column,
            @view(vector[1:left_rank]), left_slot, lifecycle, workspace.root)
        center_start = left_rank + 1
        rows = first(lifecycle.intervals[cell]):last(
            lifecycle.intervals[cell + 1])
        width = length(rows)
        copyto!(@view(output[rows, column]),
            @view(vector[center_start:center_start+width-1]))
        right_rank > 0 && _expand_suffix_mode!(output, column,
            @view(vector[end-right_rank+1:end]), right_slot, lifecycle,
            workspace.root)
    end
    for slot = 1:left_slot
        link = lifecycle.collection[slot].link
        for mode in axes(link.completed_map, 2)
            column += 1
            rows = lifecycle.intervals[slot]
            copyto!(@view(output[rows, column]),
                @view(link.completed_map[link.old_rank+1:end, mode]))
            link.old_rank > 0 && _expand_prefix_mode!(output, column,
                @view(link.completed_map[1:link.old_rank, mode]), slot-1,
                lifecycle, workspace.root)
        end
    end
    for slot = right_slot:lifecycle.cells-1
        link = lifecycle.collection[slot].link
        for mode in axes(link.completed_map, 2)
            column += 1
            rows = lifecycle.intervals[slot + 1]
            width = length(rows)
            copyto!(@view(output[rows, column]),
                @view(link.completed_map[1:width, mode]))
            child = link.physical_width
            child > 0 && _expand_suffix_mode!(output, column,
                @view(link.completed_map[width+1:end, mode]),
                slot+1, lifecycle, workspace.root)
        end
    end
    column == lifecycle.root.total_occupied || throw(ArgumentError(
        "terminal reconstruction occupation count is incomplete"))
    output
end
