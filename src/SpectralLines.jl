export transition,
    AbstractLine, GenericLine, BroadLine, NarrowLine, BroadBaseLine, MultiCompLine, LineComponent

const transitions_db = DataFrame()

function load_transitions()
    global transitions_db
    (nrow(transitions_db) > 0)  &&  (return nothing)

    d, c = csvread(joinpath(@__DIR__, "atomic_line_list.vac"), '|', commentchar='#', header_exists=true)
    out = DataFrame(collect(d), Symbol.(strip.(c)))
    for i in 1:ncol(out)
        if isa(out[1, i], String)
            out[:, i] .= strip.(out[:,i])
        end
    end
    delete!(out, findall((out.tid .== "")  .|
                         (out.Label .== "")))

    out[!, :tid] .= Symbol.(out.tid)

    append!(transitions_db, out)
    return nothing
end


function transition(tid::Symbol)
    load_transitions()
    i = findall(transitions_db.tid .== tid)
    if length(i) == 0
        error("No known transition with ID: $tid")
    end
    return transitions_db[i, :]
end





# Define spectral line types and corresponding singletons
abstract type SpecLineType end
struct Narrow    <: SpecLineType; end; const narrow    = Narrow()
struct Broad     <: SpecLineType; end; const broad     = Broad()
struct VeryBroad <: SpecLineType; end; const verybroad = VeryBroad()
struct Unknown   <: SpecLineType; end; const unknown   = Unknown()

# Emission line descriptor including transition identifier and line type decomposition.
abstract type EmLineDescription end

struct StdEmLine <: EmLineDescription
    tid::Symbol
    types::Vector{SpecLineType}
    function StdEmLine(tid::Symbol, T::Vararg{SpecLineType})
        @assert length(T) >= 1
        new(tid, [T...])
    end
end

struct CustomEmLine <: EmLineDescription
    λ::Float64
    types::Vector{SpecLineType}
    function CustomEmLine(λ::Float64, T::Vararg{SpecLineType})
        @assert length(T) >= 1
        new(λ, [T...])
    end
end

# Structure containing the actual GFit component for a single contribution to an emission line
struct EmLineComponent{SpecLineType}
    suffix::Symbol
    group::Symbol
    comp::AbstractComponent
end
