
const transitions_db = DataFrame()

function validate_transitions_db()
    @assert length(unique(transitions_db.Name))  == nrow(transitions_db)
    @assert length(unique(transitions_db.Label)) == nrow(transitions_db)
end

function load_transitions(;force=false)
    global transitions_db
    if nrow(transitions_db) > 0
        if force
            empty!(transitions_db)
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
    delete!(out, findall((out.Name .== "")  .|
                         (out.Label .== "")))
  
    out[!, :Name] .= Symbol.(out.Name)

    append!(transitions_db, out)
    validate_transitions_db()
end


function transition(name::Symbol)
    load_transitions()
    i = findall(transitions_db.Name .== name)
    @assert length(i) == 1
    return transitions_db[i[1], :]
end

function transition(λ_vacuum_ang::Float64)
    load_transitions()
    name = new_transition(λ_vacuum_ang)
    i = findall(transitions_db.Name .== name)
    @assert length(i) == 1
    return transitions_db[i[1], :]
end

function transition_default_name(λ::Float64)
    name = "T" * string(λ)
    name = join(split(name, "."), "p")
    return Symbol(name)
end

function new_transition(λ_vacuum_ang::Float64; name::Union{Symbol, Nothing}=nothing)
    load_transitions()
    isnothing(name)  &&  (name = transition_default_name(λ_vacuum_ang))
    i = findall(transitions_db.Name .== name)
    if length(i) == 0
        push!(transitions_db, [name, string(name), λ_vacuum_ang, fill("", ncol(transitions_db)-3)...])
    end
    return name
end

# Line type descriptors associated to an atomic transition
abstract type AbstractLineType end

struct BroadType       <: AbstractLineType; tid::Symbol; end
struct NarrowType      <: AbstractLineType; tid::Symbol; end
struct BroadBaseType   <: AbstractLineType; tid::Symbol; end
struct CombinedType    <: AbstractLineType; tid::Symbol; types::Vector{Type}; end


struct LineComponent
    orig::AbstractLineType
    comp::AbstractComponent
    groups::Symbol
end
