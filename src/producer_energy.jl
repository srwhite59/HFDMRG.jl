struct ProducerHamiltonian{H,V}
    one_body::H
    interaction::V
    constant::Float64
    h1_bandwidth::Int
    sites::Int
end

function ProducerHamiltonian(one_body, interaction; constant::Real=0.0,
        h1_bandwidth::Int)
    size(one_body, 1) == size(one_body, 2) || throw(DimensionMismatch(
        "producer one-body operator must be square"))
    size(interaction) == size(one_body) || throw(DimensionMismatch(
        "producer interaction dimensions disagree"))
    sites = size(one_body, 1)
    sites > 0 || throw(ArgumentError("producer must contain physical sites"))
    0 <= h1_bandwidth < sites || throw(ArgumentError(
        "producer H1 bandwidth is invalid"))
    value = Float64(constant)
    isfinite(value) || throw(ArgumentError("producer constant must be finite"))
    ProducerHamiltonian(one_body, interaction, value, h1_bandwidth, sites)
end

struct ProducerEnergy
    one_body_alpha::Float64
    one_body_beta::Float64
    hartree::Float64
    exchange_alpha::Float64
    exchange_beta::Float64
    constant::Float64
    total::Float64
    error_envelope::Float64
    projector_block_reads::Int
end


mutable struct _ProducerSum
    value::Float64
    correction::Float64
    absolute_sum::Float64
end

_ProducerSum() = _ProducerSum(0.0, 0.0, 0.0)

@inline function _producer_add!(sum::_ProducerSum, value::Float64)
    isfinite(value) || throw(ArgumentError(
        "producer-energy contribution must be finite"))
    next = sum.value + value
    sum.correction += abs(sum.value) >= abs(value) ?
        (sum.value-next)+value : (value-next)+sum.value
    sum.value = next
    sum.absolute_sum += abs(value)
    nothing
end

@inline _producer_value(sum::_ProducerSum) = sum.value + sum.correction

mutable struct ProducerEnergyWorkspace
    restricted::Bool
    sites::Int
    next_center::Int
    generation::Int
    replacements::Int
    provenance::UInt
    center_interval::UnitRange{Int}
    left_intervals::Vector{UnitRange{Int}}
    right_intervals::Vector{UnitRange{Int}}
    alpha_left::Vector{RHFStateLink}
    alpha_right::Vector{RHFStateLink}
    beta_left::Vector{RHFStateLink}
    beta_right::Vector{RHFStateLink}
    covariance::Matrix{Float64}
    parent::Matrix{Float64}
    factor::Matrix{Float64}
    cross_first::Matrix{Float64}
    cross_second::Matrix{Float64}
    projector_block::Matrix{Float64}
    density_alpha::Vector{Float64}
    density_beta::Vector{Float64}
    projector_block_reads::Int
end

@inline _producer_lifecycle(result::HFDMRGResult) = result.state.lifecycle
@inline _producer_lifecycle(result::UHFDMRGResult) = result.state.lifecycle

@inline _producer_spin_block(lifecycle::RHFFixedLifecycle, slot, alpha) =
    lifecycle.collection[slot]
@inline _producer_spin_block(lifecycle::UHFFixedLifecycle, slot, alpha) =
    alpha ? lifecycle.collection[slot].alpha : lifecycle.collection[slot].beta

function _producer_link_inventory(lifecycle, alpha::Bool)
    cell = lifecycle.next_center
    left = RHFStateLink[_producer_spin_block(lifecycle, slot, alpha).link
        for slot = cell-1:-1:1]
    right = RHFStateLink[_producer_spin_block(lifecycle, slot, alpha).link
        for slot = cell+1:lifecycle.cells-1]
    left_intervals = UnitRange{Int}[lifecycle.intervals[slot]
        for slot = cell-1:-1:1]
    right_intervals = UnitRange{Int}[lifecycle.intervals[slot+1]
        for slot = cell+1:lifecycle.cells-1]
    left, right, left_intervals, right_intervals
end

function ProducerEnergyWorkspace(result::Union{HFDMRGResult,UHFDMRGResult},
        producer::ProducerHamiltonian)
    lifecycle = _producer_lifecycle(result)
    producer.sites == last(lifecycle.intervals[end]) ||
        throw(DimensionMismatch("producer and compact state lengths disagree"))
    restricted = lifecycle isa RHFFixedLifecycle
    alpha_left, alpha_right, left_intervals, right_intervals =
        _producer_link_inventory(lifecycle, true)
    beta_left, beta_right, beta_left_intervals, beta_right_intervals =
        restricted ? (alpha_left, alpha_right, left_intervals,
            right_intervals) : _producer_link_inventory(lifecycle, false)
    left_intervals == beta_left_intervals &&
        right_intervals == beta_right_intervals || throw(ArgumentError(
        "alpha and beta producer intervals disagree"))
    center = first(lifecycle.intervals[lifecycle.next_center]):last(
        lifecycle.intervals[lifecycle.next_center+1])
    maximum_coordinate = max(
        restricted ? lifecycle.root.rank : lifecycle.root.alpha.rank,
        restricted ? lifecycle.root.rank : lifecycle.root.beta.rank, 1)
    for links in (alpha_left, alpha_right, beta_left, beta_right), link in links
        maximum_coordinate = max(maximum_coordinate,
            link.old_rank + link.physical_width)
    end
    maximum_width = length(center)
    for intervals in (left_intervals, right_intervals), interval in intervals
        maximum_width = max(maximum_width, length(interval))
    end
    square = () -> Matrix{Float64}(undef, maximum_coordinate,
        maximum_coordinate)
    ProducerEnergyWorkspace(restricted, producer.sites,
        lifecycle.next_center, lifecycle.generation, lifecycle.replacements,
        lifecycle.provenance, center, left_intervals, right_intervals,
        alpha_left, alpha_right, beta_left, beta_right,
        square(), square(), square(),
        Matrix{Float64}(undef, maximum_width, maximum_coordinate),
        Matrix{Float64}(undef, maximum_width, maximum_coordinate),
        Matrix{Float64}(undef, maximum_width, maximum_width),
        Vector{Float64}(undef, producer.sites),
        Vector{Float64}(undef, producer.sites), 0)
end

@inline _producer_root(lifecycle::RHFFixedLifecycle, alpha) = lifecycle.root
@inline _producer_root(lifecycle::UHFFixedLifecycle, alpha) =
    alpha ? lifecycle.root.alpha : lifecycle.root.beta

function _check_producer_workspace(workspace::ProducerEnergyWorkspace,
        lifecycle, producer::ProducerHamiltonian)
    producer.sites == workspace.sites == last(lifecycle.intervals[end]) ||
        throw(DimensionMismatch("producer workspace length is stale"))
    workspace.next_center == lifecycle.next_center &&
        workspace.generation == lifecycle.generation &&
        workspace.replacements == lifecycle.replacements &&
        workspace.provenance == lifecycle.provenance || throw(ArgumentError(
        "producer workspace compact-state provenance is stale"))
    workspace.restricted == (lifecycle isa RHFFixedLifecycle) ||
        throw(ArgumentError("producer workspace spin mode is stale"))
    _check_moving_root(_producer_root(lifecycle, true))
    workspace.restricted || _check_moving_root(_producer_root(lifecycle,
        false))
    cell = lifecycle.next_center
    left_count, right_count = cell-1, lifecycle.cells-cell-1
    length(workspace.alpha_left) == left_count &&
        length(workspace.alpha_right) == right_count &&
        length(workspace.beta_left) == left_count &&
        length(workspace.beta_right) == right_count &&
        length(workspace.left_intervals) == left_count &&
        length(workspace.right_intervals) == right_count ||
        throw(ArgumentError("producer workspace link inventory is stale"))
    for alpha in (workspace.restricted ? (true,) : (true, false))
        left = alpha ? workspace.alpha_left : workspace.beta_left
        right = alpha ? workspace.alpha_right : workspace.beta_right
        for (ordinal, slot) in enumerate(cell-1:-1:1)
            actual = _producer_spin_block(lifecycle, slot, alpha).link
            left[ordinal].close_map === actual.close_map &&
                left[ordinal].completed_map === actual.completed_map ||
                throw(ArgumentError("producer workspace left links are stale"))
            workspace.left_intervals[ordinal] == lifecycle.intervals[slot] ||
                throw(ArgumentError("producer workspace left interval is stale"))
        end
        for (ordinal, slot) in enumerate(cell+1:lifecycle.cells-1)
            actual = _producer_spin_block(lifecycle, slot, alpha).link
            right[ordinal].close_map === actual.close_map &&
                right[ordinal].completed_map === actual.completed_map ||
                throw(ArgumentError("producer workspace right links are stale"))
            workspace.right_intervals[ordinal] ==
                lifecycle.intervals[slot+1] || throw(ArgumentError(
                "producer workspace right interval is stale"))
        end
    end
    maximum_width = length(workspace.center_interval)
    maximum_coordinate = max(_producer_root(lifecycle, true).rank,
        workspace.restricted ? 0 : _producer_root(lifecycle, false).rank)
    for intervals in (workspace.left_intervals, workspace.right_intervals),
            interval in intervals
        maximum_width = max(maximum_width, length(interval))
    end
    for links in (workspace.alpha_left, workspace.alpha_right,
            workspace.beta_left, workspace.beta_right), link in links
        maximum_coordinate = max(maximum_coordinate,
            link.old_rank + link.physical_width)
    end
    for array in (workspace.covariance, workspace.parent, workspace.factor)
        size(array, 1) >= maximum_coordinate &&
            size(array, 2) >= maximum_coordinate || throw(DimensionMismatch(
            "producer coordinate capacity is too small"))
    end
    for array in (workspace.cross_first, workspace.cross_second)
        size(array, 1) >= maximum_width &&
            size(array, 2) >= maximum_coordinate || throw(DimensionMismatch(
            "producer cross-block capacity is too small"))
    end
    size(workspace.projector_block, 1) >= maximum_width &&
        size(workspace.projector_block, 2) >= maximum_width ||
        throw(DimensionMismatch("producer projector-block capacity is too small"))
    length(workspace.density_alpha) == workspace.sites &&
        length(workspace.density_beta) == workspace.sites ||
        throw(DimensionMismatch("producer density capacity is too small"))
    nothing
end

@inline function _link_ranges(link::RHFStateLink, left::Bool)
    left ? (1:link.old_rank,
        link.old_rank+1:link.old_rank+link.physical_width) :
        (link.old_rank+1:link.old_rank+link.physical_width,
            1:link.old_rank)
end

function _open_producer_covariance!(workspace::ProducerEnergyWorkspace,
        active::Int, link::RHFStateLink, left::Bool)
    link.active_rank == active || throw(DimensionMismatch(
        "producer links do not compose"))
    parent_rank = link.old_rank + link.physical_width
    parent = @view workspace.parent[1:parent_rank, 1:parent_rank]
    fill!(parent, 0.0)
    if active > 0
        intermediate = @view workspace.factor[1:parent_rank, 1:active]
        mul!(intermediate, link.close_map,
            @view(workspace.covariance[1:active, 1:active]))
        mul!(parent, intermediate, transpose(link.close_map))
    end
    if !isempty(link.completed_map)
        mul!(parent, link.completed_map, transpose(link.completed_map),
            1.0, 1.0)
    end
    old, physical = _link_ranges(link, left)
    parent, old, physical
end

function _producer_add_block!(one_body::_ProducerSum,
        exchange::_ProducerSum, first_interval, second_interval, block,
        producer::ProducerHamiltonian, same::Bool)
    @inbounds for local_j in axes(block, 2), local_i in axes(block, 1)
        i = first(first_interval) + local_i - 1
        j = first(second_interval) + local_j - 1
        p = block[local_i, local_j]
        if same
            abs(i-j) <= producer.h1_bandwidth && _producer_add!(one_body,
                Float64(producer.one_body[i, j]) * p)
            _producer_add!(exchange,
                -0.5 * Float64(producer.interaction[i, j]) * p^2)
        else
            abs(i-j) <= producer.h1_bandwidth && _producer_add!(one_body,
                (Float64(producer.one_body[i, j]) +
                    Float64(producer.one_body[j, i])) * p)
            _producer_add!(exchange, -0.5 *
                (Float64(producer.interaction[i, j]) +
                    Float64(producer.interaction[j, i])) * p^2)
        end
    end
    nothing
end

function _producer_branch!(one_body, exchange, density, root_covariance,
        boundary, links, intervals, left::Bool, producer, workspace)
    for target in eachindex(links)
        active = length(boundary)
        active > 0 && copyto!(@view(workspace.covariance[1:active, 1:active]),
            @view(root_covariance[boundary, boundary]))
        for index = 1:target
            parent, old, physical = _open_producer_covariance!(workspace,
                active, links[index], left)
            if index == target
                interval = intervals[index]
                block = @view parent[physical, physical]
                _producer_add_block!(one_body, exchange, interval, interval,
                    block, producer, true)
                @inbounds for local_row in eachindex(interval)
                    density[first(interval)+local_row-1] =
                        block[local_row, local_row]
                end
                width = length(interval); next_active = length(old)
                next_active > 0 && copyto!(@view(workspace.cross_first[
                    1:width, 1:next_active]), @view(parent[physical, old]))
                active = next_active
            else
                active = length(old)
                active > 0 && copyto!(@view(workspace.covariance[
                    1:active, 1:active]), @view(parent[old, old]))
            end
        end
        width = length(intervals[target])
        for source = target+1:length(links)
            link = links[source]
            link.active_rank == active || throw(DimensionMismatch(
                "producer branch cross covariance does not compose"))
            old, physical = _link_ranges(link, left)
            pair = @view workspace.projector_block[1:width,
                1:length(physical)]
            if active == 0
                fill!(pair, 0.0)
            else
                mul!(pair, @view(workspace.cross_first[1:width, 1:active]),
                    transpose(@view(link.close_map[physical, 1:active])))
            end
            _producer_add_block!(one_body, exchange, intervals[target],
                intervals[source], pair, producer, false)
            next_active = length(old)
            if next_active > 0
                next = @view workspace.cross_second[1:width, 1:next_active]
                active == 0 ? fill!(next, 0.0) : mul!(next,
                    @view(workspace.cross_first[1:width, 1:active]),
                    transpose(@view(link.close_map[old, 1:active])))
                copyto!(@view(workspace.cross_first[1:width,
                    1:next_active]), next)
            end
            active = next_active
        end
    end
    nothing
end

function _producer_center_branch!(one_body, exchange, root_covariance,
        center_range, boundary, links, intervals, left, producer, workspace)
    width = length(center_range); active = length(boundary)
    active > 0 && copyto!(@view(workspace.cross_first[1:width, 1:active]),
        @view(root_covariance[center_range, boundary]))
    for index in eachindex(links)
        link = links[index]
        link.active_rank == active || throw(DimensionMismatch(
            "producer center cross covariance does not compose"))
        old, physical = _link_ranges(link, left)
        pair = @view workspace.projector_block[1:width, 1:length(physical)]
        active == 0 ? fill!(pair, 0.0) : mul!(pair,
            @view(workspace.cross_first[1:width, 1:active]),
            transpose(@view(link.close_map[physical, 1:active])))
        _producer_add_block!(one_body, exchange, workspace.center_interval,
            intervals[index], pair, producer, false)
        next_active = length(old)
        if next_active > 0
            next = @view workspace.cross_second[1:width, 1:next_active]
            active == 0 ? fill!(next, 0.0) : mul!(next,
                @view(workspace.cross_first[1:width, 1:active]),
                transpose(@view(link.close_map[old, 1:active])))
            copyto!(@view(workspace.cross_first[1:width, 1:next_active]),
                next)
        end
        active = next_active
    end
    nothing
end

function _producer_cross_branches!(one_body, exchange, root_covariance,
        left_boundary, right_boundary, left_links, right_links,
        left_intervals, right_intervals, producer, workspace)
    for target in eachindex(left_links)
        rows = length(left_intervals[target])
        left_active, right_active = length(left_boundary),
            length(right_boundary)
        left_active > 0 && right_active > 0 && copyto!(
            @view(workspace.factor[1:left_active, 1:right_active]),
            @view(root_covariance[left_boundary, right_boundary]))
        for index = 1:target
            link = left_links[index]
            old, physical = _link_ranges(link, true)
            if index == target
                output = @view workspace.cross_first[1:rows, 1:right_active]
                left_active == 0 || right_active == 0 ? fill!(output, 0.0) :
                    mul!(output,
                        @view(link.close_map[physical, 1:left_active]),
                        @view(workspace.factor[1:left_active,
                            1:right_active]))
            else
                output = @view workspace.parent[1:length(old),
                    1:right_active]
                left_active == 0 || right_active == 0 ? fill!(output, 0.0) :
                    mul!(output, @view(link.close_map[old, 1:left_active]),
                        @view(workspace.factor[1:left_active,
                            1:right_active]))
                copyto!(@view(workspace.factor[1:length(old),
                    1:right_active]), output)
                left_active = length(old)
            end
        end
        active = right_active
        for source in eachindex(right_links)
            link = right_links[source]
            old, physical = _link_ranges(link, false)
            pair = @view workspace.projector_block[1:rows,
                1:length(physical)]
            active == 0 ? fill!(pair, 0.0) : mul!(pair,
                @view(workspace.cross_first[1:rows, 1:active]),
                transpose(@view(link.close_map[physical, 1:active])))
            _producer_add_block!(one_body, exchange,
                left_intervals[target], right_intervals[source], pair,
                producer, false)
            next_active = length(old)
            if next_active > 0
                next = @view workspace.cross_second[1:rows, 1:next_active]
                active == 0 ? fill!(next, 0.0) : mul!(next,
                    @view(workspace.cross_first[1:rows, 1:active]),
                    transpose(@view(link.close_map[old, 1:active])))
                copyto!(@view(workspace.cross_first[1:rows,
                    1:next_active]), next)
            end
            active = next_active
        end
    end
    nothing
end

function _producer_spin!(one_body, exchange, density, lifecycle,
        alpha::Bool, producer, workspace)
    root = _producer_root(lifecycle, alpha)
    covariance = @view root.covariance[1:root.rank, 1:root.rank]
    left_links = alpha ? workspace.alpha_left : workspace.beta_left
    right_links = alpha ? workspace.alpha_right : workspace.beta_right
    left_rank = isempty(left_links) ? 0 : first(left_links).active_rank
    right_rank = isempty(right_links) ? 0 : first(right_links).active_rank
    center_width = length(workspace.center_interval)
    root.rank == left_rank + center_width + right_rank ||
        throw(DimensionMismatch("producer root and center dimensions disagree"))
    left_boundary = 1:left_rank
    center_range = left_rank+1:left_rank+center_width
    right_boundary = left_rank+center_width+1:root.rank
    center = @view covariance[center_range, center_range]
    _producer_add_block!(one_body, exchange, workspace.center_interval,
        workspace.center_interval, center, producer, true)
    @inbounds for local_row in eachindex(workspace.center_interval)
        density[first(workspace.center_interval)+local_row-1] =
            center[local_row, local_row]
    end
    _producer_branch!(one_body, exchange, density, covariance, left_boundary,
        left_links, workspace.left_intervals, true, producer, workspace)
    _producer_branch!(one_body, exchange, density, covariance, right_boundary,
        right_links, workspace.right_intervals, false, producer, workspace)
    _producer_center_branch!(one_body, exchange, covariance, center_range,
        left_boundary, left_links, workspace.left_intervals, true, producer,
        workspace)
    _producer_center_branch!(one_body, exchange, covariance, center_range,
        right_boundary, right_links, workspace.right_intervals, false,
        producer, workspace)
    _producer_cross_branches!(one_body, exchange, covariance, left_boundary,
        right_boundary, left_links, right_links, workspace.left_intervals,
        workspace.right_intervals, producer, workspace)
    nothing
end

function _producer_hartree!(sum, alpha, beta, producer)
    @inbounds for i in eachindex(alpha)
        charge_i = alpha[i] + beta[i]
        for j = i:length(alpha)
            charge_j = alpha[j] + beta[j]
            value = if i == j
                0.5 * Float64(producer.interaction[i, i]) * charge_i^2
            else
                0.5 * (Float64(producer.interaction[i, j]) +
                    Float64(producer.interaction[j, i])) * charge_i*charge_j
            end
            _producer_add!(sum, value)
        end
    end
    nothing
end

function producer_energy(result::Union{HFDMRGResult,UHFDMRGResult},
        producer::ProducerHamiltonian,
        workspace::ProducerEnergyWorkspace=ProducerEnergyWorkspace(
            result, producer))
    lifecycle = _producer_lifecycle(result)
    _check_producer_workspace(workspace, lifecycle, producer)
    reads_before = workspace.projector_block_reads
    alpha_h1, beta_h1 = _ProducerSum(), _ProducerSum()
    alpha_exchange, beta_exchange = _ProducerSum(), _ProducerSum()
    hartree = _ProducerSum()
    _producer_spin!(alpha_h1, alpha_exchange, workspace.density_alpha,
        lifecycle, true, producer, workspace)
    if workspace.restricted
        copyto!(workspace.density_beta, workspace.density_alpha)
        beta_h1.value = alpha_h1.value
        beta_h1.correction = alpha_h1.correction
        beta_h1.absolute_sum = alpha_h1.absolute_sum
        beta_exchange.value = alpha_exchange.value
        beta_exchange.correction = alpha_exchange.correction
        beta_exchange.absolute_sum = alpha_exchange.absolute_sum
    else
        _producer_spin!(beta_h1, beta_exchange, workspace.density_beta,
            lifecycle, false, producer, workspace)
    end
    _producer_hartree!(hartree, workspace.density_alpha,
        workspace.density_beta, producer)
    ha, hb = _producer_value(alpha_h1), _producer_value(beta_h1)
    ja = _producer_value(hartree)
    ka, kb = _producer_value(alpha_exchange),
        _producer_value(beta_exchange)
    total = producer.constant + ha + hb + ja + ka + kb
    scale = abs(producer.constant) + alpha_h1.absolute_sum +
        beta_h1.absolute_sum + hartree.absolute_sum +
        alpha_exchange.absolute_sum + beta_exchange.absolute_sum
    envelope = 256eps(Float64) * max(scale, 1.0)
    pieces = 1 + length(workspace.left_intervals) +
        length(workspace.right_intervals)
    reads = pieces*(pieces+1)÷2 * (workspace.restricted ? 1 : 2)
    workspace.projector_block_reads += reads
    ProducerEnergy(ha, hb, ja, ka, kb, producer.constant, total, envelope,
        workspace.projector_block_reads-reads_before)
end
