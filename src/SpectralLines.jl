export StdEmLine, broad, narrow, verybroad, unknown

struct SpectralTransition
    id::Symbol
    label::String
    λ::Float64  # vacuum wavelength in Angstrom
    element::String
    ttype::String
    config::String
    term::String
    JJ::String
    levels::NTuple{2, Float64}
    notes::String
end

const known_transitions = OrderedDict{Symbol, Vector{SpectralTransition}}()
function transitions(tid::Symbol)
    global known_transitions

    if length(known_transitions) == 0
        d, c = csvread(joinpath(@__DIR__, "atomic_line_list.vac"), '|', commentchar='#', header_exists=true)
        for i in 1:length(d)
            isa(d[i][1], AbstractString)  ||  continue
            d[i] .= strip.(d[i])
        end

        db = Vector{SpectralTransition}()
        for i in 1:length(d[1])
            ll = split(d[9][i], '-')
            tt = SpectralTransition(Symbol(d[1][i]),
                                    d[2][i], d[3][i], d[4][i],
                                    d[5][i], d[6][i], d[7][i], d[8][i],
                                    (Meta.parse(ll[1]), Meta.parse(ll[2])),
                                    d[10][i])
            push!(db, tt)
        end
        db = db[sortperm(getfield.(db, :λ))]

        for tt in db
            if haskey(known_transitions, tt.id)
                push!(known_transitions[tt.id], tt)
            else
                known_transitions[tt.id] = [tt]
            end
        end
    end

    @assert haskey(known_transitions, tid) "Unknown transition identifier: $tid"
    return known_transitions[tid]
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
