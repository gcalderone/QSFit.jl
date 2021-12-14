export transition, custom_transition,
    AbstractLine, GenericLine, BroadLine, NarrowLine, BroadBaseLine, AsymmTailLine, MultiCompLine, LineComponent

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


function custom_transition_default_tid(λ::Float64)
    tid = "T" * string(λ)
    tid = join(split(tid, "."), "p")
    return Symbol(tid)
end


function custom_transition(λ_vac_ang::Float64; tid::Union{Symbol, Nothing}=nothing)
    load_transitions()
    isnothing(tid)  &&  (tid = custom_transition_default_tid(λ_vac_ang))
    i = findall(transitions_db.tid .== tid)
    if length(i) > 0
        delete!(transitions_db, i)
    end
    t = [tid, string(tid), λ_vac_ang, fill("", ncol(transitions_db)-3)...]
    push!(transitions_db, t)
    return tid
end

transition(λ_vac_ang::Float64) = transition(custom_transition_default_tid(λ_vac_ang))



# Line type descriptors associated to a specific transition
abstract type AbstractLine end

macro define_line(name, suffix, group)
    return esc(:(
        struct $name <: AbstractLine;
        cname::Symbol;
        tid::Symbol;
        $name(tid; cname=tid) = new(Symbol(cname), tid);
        end;
        suffix(::$(name), multicomp::Bool) = multicomp  ?  $(QuoteNode(suffix))  :  Symbol("");
        group( ::$(name)) = $(QuoteNode(group));
    ))
end

@define_line GenericLine   _gen GenericLines
@define_line BroadLine     _br  BroadLines
@define_line NarrowLine    _na  NarrowLines
@define_line BroadBaseLine _bb  BroadBaseLines

struct AsymmTailLine <: AbstractLine
    cname::Symbol
    tid::Symbol
    side::Symbol

    function AsymmTailLine(tid::Symbol, side::Symbol; cname=tid)
        @assert side in (:red, :blue)
        new(Symbol(cname), tid, side)
    end
end
suffix(line::AsymmTailLine, multicomp::Bool) = line.side == :blue ?  :_bw  :  :_rw
group(::AsymmTailLine) = :AsymmTailLines

struct MultiCompLine <: AbstractLine
    cname::Symbol
    tid::Symbol
    types::Vector{Type}
    MultiCompLine(tid::Symbol, types::Vector{Type}; cname=tid) = new(Symbol(cname), tid, types)
end

struct LineComponent
    line::AbstractLine
    comp::AbstractComponent
    multicomp::Bool
end

function collect_LineComponent(source::QSO)
    out = OrderedDict{Symbol, LineComponent}()
    for line in known_spectral_lines(source)
        if isa(line, MultiCompLine)
            for t in line.types
                nn = Symbol(line.cname, suffix(t(line.tid), true))
                (nn in source.options[:skip_lines])  &&  continue
                @assert !haskey(out, nn)
                lc = LineComponent(source, t(line.tid), true)
                out[nn] = lc
            end
        else
            nn = Symbol(line.cname, suffix(line, false))
            (nn in source.options[:skip_lines])  &&  continue
            @assert !haskey(out, nn)
            out[nn] = LineComponent(source, line, false)
        end
    end
    return out
end


LineComponent(source::QSO, line::AbstractLine, multicomp::Bool) =
    error("No constructor method has been defined for LineComponent($(typeof(source)), $(typeof(line)), ::Bool)")
