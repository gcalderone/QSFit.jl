export StdEmLine, broad, narrow, verybroad, nuisance

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

    @assert haskey(known_transitions, tid) "Nuisance transition identifier: $tid"
    return known_transitions[tid]
end


# Emission line descriptor including transition identifier and line type decomposition.
abstract type EmLineDescription end

struct StdEmLine <: EmLineDescription
    tid::Symbol
    types::Vector{Val}
    function StdEmLine(tid::Symbol, T::Vararg{Symbol})
        @assert length(T) >= 1
        new(tid, [Val.(T)...])
    end
end

struct CustomEmLine <: EmLineDescription
    λ::Float64
    types::Vector{Val}
    function CustomEmLine(λ::Float64, T::Vararg{Symbol})
        @assert length(T) >= 1
        new(λ, [Val.(T)...])
    end
end

# Structure containing the actual GModelFit component for a single contribution to an emission line
struct EmLineComponent{Val}
    suffix::Symbol
    group::Symbol
    comp::AbstractComponent
    EmLineComponent{T}(suffix::Symbol, group::Symbol, comp::AbstractComponent) where T =
        new{T}(suffix, group, comp)
end
