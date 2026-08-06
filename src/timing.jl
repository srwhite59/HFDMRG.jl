"""
Private hierarchical wall-clock instrumentation for HFDMRG.

Collection and live printing are disabled by default. `@timeg` evaluates its
expression directly when disabled; enabled scopes are task-local and nested.
Repeated labels are merged when reports are rendered, preserving exact call
counts and separating inclusive `elapsed_seconds` from child-exclusive
`self_seconds`.
"""
module TimeG

export @timeg,
       timing_enabled,
       timing_live_enabled,
       set_timing!,
       set_timing_live!,
       set_timing_thresholds!,
       reset_timing_report!,
       current_timing_report,
       timing_report

struct TimingConfig
    enabled::Bool
    live_print::Bool
    expand_threshold_seconds::Float64
    drop_threshold_seconds::Float64
end

mutable struct TimingNode
    label::String
    elapsed_seconds::Float64
    self_seconds::Float64
    call_count::Int
    children::Vector{TimingNode}
end

mutable struct _TimingFrame
    label::String
    start_ns::UInt64
    children::Vector{TimingNode}
    child_elapsed_seconds::Float64
end

mutable struct _TimingState
    stack::Vector{_TimingFrame}
    roots::Vector{TimingNode}
end

struct TimingReport
    roots::Vector{TimingNode}
end

Base.show(io::IO, report::TimingReport) =
    print(io, "TimingReport(nroots=", length(report.roots), ")")

const _TIMING_CONFIG = Ref(TimingConfig(false, false, 10.0, 0.1))
const _STATE_KEY = :hfdmrg_timing_state

function _timing_state()
    storage = task_local_storage()
    if !haskey(storage, _STATE_KEY)
        state = _TimingState(_TimingFrame[], TimingNode[])
        task_local_storage(_STATE_KEY, state)
        return state
    end
    storage[_STATE_KEY]
end

timing_enabled() = _TIMING_CONFIG[].enabled
timing_live_enabled() = _TIMING_CONFIG[].live_print

function set_timing!(enabled::Bool)
    config = _TIMING_CONFIG[]
    _TIMING_CONFIG[] = TimingConfig(enabled, config.live_print,
        config.expand_threshold_seconds, config.drop_threshold_seconds)
    enabled
end

function set_timing_live!(enabled::Bool)
    config = _TIMING_CONFIG[]
    _TIMING_CONFIG[] = TimingConfig(config.enabled, enabled,
        config.expand_threshold_seconds, config.drop_threshold_seconds)
    enabled
end

function set_timing_thresholds!(;
        expand::Real = _TIMING_CONFIG[].expand_threshold_seconds,
        drop::Real = _TIMING_CONFIG[].drop_threshold_seconds)
    expand_value, drop_value = Float64(expand), Float64(drop)
    expand_value >= 0 || throw(ArgumentError("timing expand threshold must be nonnegative"))
    drop_value >= 0 || throw(ArgumentError("timing drop threshold must be nonnegative"))
    config = _TIMING_CONFIG[]
    _TIMING_CONFIG[] = TimingConfig(config.enabled, config.live_print,
        expand_value, drop_value)
end

function reset_timing_report!()
    state = _timing_state()
    empty!(state.stack); empty!(state.roots)
    nothing
end

function current_timing_report()
    state = _timing_state()
    isempty(state.stack) || throw(ArgumentError(
        "cannot snapshot timing report while timed scopes are active"))
    TimingReport(deepcopy(state.roots))
end

function _begin_scope(label)
    state = _timing_state()
    frame = _TimingFrame(String(label), time_ns(), TimingNode[], 0.0)
    push!(state.stack, frame)
    frame
end

function _finish_scope(frame)
    state = _timing_state()
    finished = pop!(state.stack)
    finished === frame || error("timing stack closed out of order")
    elapsed = (time_ns() - finished.start_ns) * 1.0e-9
    node = TimingNode(finished.label, elapsed,
        max(0.0, elapsed - finished.child_elapsed_seconds), 1, finished.children)
    if isempty(state.stack)
        push!(state.roots, node)
    else
        parent = state.stack[end]
        push!(parent.children, node)
        parent.child_elapsed_seconds += elapsed
    end
    config = _TIMING_CONFIG[]
    if config.live_print && elapsed >= config.drop_threshold_seconds
        println("  "^length(state.stack), finished.label, ": ",
            round(elapsed; sigdigits = 4), " seconds")
    end
    nothing
end

function _time_scope(f, label)
    !timing_enabled() && return f()
    frame = _begin_scope(label)
    try
        f()
    finally
        _finish_scope(frame)
    end
end

"""
    @timeg label expression

Record a named hierarchical wall-clock scope when timing is enabled. The
disabled branch evaluates `expression` inline, avoiding closure capture in the
sweep loop and preserving warmed allocation behavior.
"""
macro timeg(label, expression)
    enabled = GlobalRef(@__MODULE__, :timing_enabled)
    begin_scope = GlobalRef(@__MODULE__, :_begin_scope)
    finish_scope = GlobalRef(@__MODULE__, :_finish_scope)
    quote
        if $enabled()
            local frame = $begin_scope($(esc(label)))
            try
                $(esc(expression))
            finally
                $finish_scope(frame)
            end
        else
            $(esc(expression))
        end
    end
end

function _timing_merge_nodes(nodes::AbstractVector{TimingNode})
    grouped = Dict{String,TimingNode}(); order = String[]
    for node in nodes
        if !haskey(grouped, node.label)
            grouped[node.label] = TimingNode(node.label, 0.0, 0.0, 0, TimingNode[])
            push!(order, node.label)
        end
        merged = grouped[node.label]
        merged.elapsed_seconds += node.elapsed_seconds
        merged.self_seconds += node.self_seconds
        merged.call_count += node.call_count
        append!(merged.children, node.children)
    end
    [grouped[label] for label in order]
end

function _timing_report_lines!(lines, nodes, config, indent)
    for node in _timing_merge_nodes(nodes)
        node.elapsed_seconds >= config.drop_threshold_seconds || continue
        push!(lines, string("  "^indent, node.label, ": elapsed=",
            node.elapsed_seconds, " s, self=", node.self_seconds,
            " s, calls=", node.call_count))
        node.elapsed_seconds >= config.expand_threshold_seconds &&
            _timing_report_lines!(lines, node.children, config, indent + 1)
    end
    lines
end

function timing_report(io::IO, report::TimingReport = current_timing_report())
    println(io, "HFDMRG timing report")
    lines = _timing_report_lines!(String[], report.roots, _TIMING_CONFIG[], 0)
    foreach(line -> println(io, line), lines)
    report
end

timing_report(report::TimingReport = current_timing_report()) =
    sprint(io -> timing_report(io, report))

end
