export transition,
    AbstractLine, BroadLine, NarrowLine, BroadBaseLine, CombinedLine,
    LineComponent

const transitions_db = DataFrame()

function validate_transitions_db()
    @assert length(unique(transitions_db.tid))   == nrow(transitions_db)
    @assert length(unique(transitions_db.Label)) == nrow(transitions_db)
end

function load_transitions(;force=false)
    global transitions_db
    if nrow(transitions_db) > 0
        if force
            empty!(transitions_db)
            select!(transitions_db, Not(1:ncol(transitions_db)))
        else
            validate_transitions_db()
            return nothing
        end
    end

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
    validate_transitions_db()
end


function transition(tid::Symbol)
    load_transitions()
    i = findall(transitions_db.tid .== tid)
    @assert length(i) == 1
    return transitions_db[i[1], :]
end


function custom_transition_default_tid(λ::Float64)
    tid = "T" * string(λ)
    tid = join(split(tid, "."), "p")
    return Symbol(tid)
end

function custom_transition(λ_vac_ang::Float64; tid::Union{Symbol, Nothing}=nothing)
    load_transitions()
    isnothing(tid)  &&  (tid = custom_transition_default_tid(λ_vac_ang))
    i = findall(transitions_db.tid .== tid)
    if length(i) == 0
        push!(transitions_db, [tid, string(tid), λ_vac_ang, fill("", ncol(transitions_db)-3)...])
    end
    return tid
end

transition(λ_vac_ang::Float64) = transition(custom_transition_default_tid(λ_vac_ang))



# Line type descriptors associated to an atomic transition
abstract type AbstractLine end

struct BroadLine       <: AbstractLine; tid::Symbol; end
struct NarrowLine      <: AbstractLine; tid::Symbol; end
struct BroadBaseLine   <: AbstractLine; tid::Symbol; end
struct CombinedLine    <: AbstractLine; tid::Symbol; types::Vector{Type}; end
struct UnkLine         <: AbstractLine; tid::Symbol; end

struct LineComponent
    orig::AbstractLine
    comp::AbstractComponent
    group::Symbol
end
