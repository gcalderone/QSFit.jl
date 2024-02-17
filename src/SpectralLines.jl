export StdEmLine, broad, narrow, verybroad, nuisance

struct SpectralLine
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

const known_transitions = OrderedDict{Symbol, Vector{SpectralLine}}()
function transitions(tid::Symbol)
    global known_transitions

    if length(known_transitions) == 0
        d, c = csvread(joinpath(@__DIR__, "atomic_line_list.vac"), '|', commentchar='#', header_exists=true)
        for i in 1:length(d)
            isa(d[i][1], AbstractString)  ||  continue
            d[i] .= strip.(d[i])
        end

        db = Vector{SpectralLine}()
        for i in 1:length(d[1])
            ll = split(d[9][i], '-')
            tt = SpectralLine(Symbol(d[1][i]),
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


abstract type AbstractLine end
abstract type ForbiddenLine <: AbstractLine  end
abstract type PermittedLine <: AbstractLine  end
abstract type NarrowLine    <: PermittedLine end
abstract type BroadLine     <: PermittedLine end
abstract type VeryBroadLine <: BroadLine     end
abstract type NuisanceLine  <: AbstractLine  end

struct LineDescriptor{T}
    id::T
    linetypes::Vector{DataType}
    LineDescriptor(id::Symbol, ::Type{T}) where T <: AbstractLine = new{Symbol}(id, [T])
    LineDescriptor(λ::Real   , ::Type{T}) where T <: AbstractLine = new{Float64}(float(λ), [T])
    function LineDescriptor(id::Symbol, types::Vararg{DataType, N}) where N
        @assert N > 0
        @assert all(isa.(types, Type{<: AbstractLine}))
        new{Symbol}(id, [types...])
    end
    function LineDescriptor(λ::Real, types::Vararg{DataType, N}) where N
        @assert N > 0
        @assert all(isa.(types, Type{<: AbstractLine}))
        new{Float64}(float(λ), [types...])
    end
end

struct LineComponent
    id::Symbol
    linetype::DataType
    group::Symbol
    comp::AbstractComponent
end


line_suffix(::RRef{R}, ::Type{T}) where {R <: AbstractRecipe, T <: AbstractLine} =
    error("No line_suffix() method defined for recipe $R and type $T")

line_group(::RRef{R}, ::Type{T}) where {R <: AbstractRecipe, T <: AbstractLine} =
    error("No line_group() method defined for recipe $R and type $T")

line_descriptors(recipe::RRef{R}) where {R <: AbstractRecipe} =
    error("No line_descriptors() method defined for recipe $R and type $T")

line_component(recipe::RRef{R}, ::Type{T}, λ::Float64) where {R <: AbstractRecipe, T <: AbstractLine} =
    error("No line_component() method defined for recipe $R and type $T")


function line_components(recipe::RRef)
    lcs = OrderedDict{Symbol, LineComponent}()
    for line in line_descriptors(recipe)
        if isa(line, LineDescriptor{Symbol})
            tt = transitions(line.id)
            λ = getfield.(tt, :λ)
            if length(λ) > 1
                # TODO: take spectral resolution into account
                @warn "Considering average wavelength for the $(line.id) multiplet: " * string(mean(λ)) * "Å"
                λ = mean(λ)  # average lambda of multiplets
            else
                λ = λ[1]
            end
            id = line.id
        else
            @assert isa(line, LineDescriptor{Float64})
            λ = line.id
            id = Symbol(:l, λ, :A)
        end

        for linetype in line.linetypes
            cname = Symbol(id, line_suffix(recipe, linetype))
            comp  = line_component(recipe, linetype, λ)
            @assert !(cname in keys(lcs)) "Duplicate component name: $cname"
            lcs[cname] = LineComponent(id, linetype, line_group(recipe, linetype), comp)
        end
    end
    return lcs
end
