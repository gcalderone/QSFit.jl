export AbstractLine, ForbiddenLine, NarrowLine, BroadLine, VeryBroadLine, NuisanceLine, LineDescriptor

export line_group, line_cname, line_component

# ====================================================================
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


struct Multiplet
    list::Vector{SpectralTransition}
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
    out = known_transitions[tid]

    if length(out) == 1
        return out[1]
    end
    @warn "$tid is a multiplet"
    return Multiplet(out)
end


vacuum_wavelength(t::SpectralTransition) = t.λ
vacuum_wavelength(t::Multiplet) = mean(vacuum_wavelength.(t.list))
vacuum_wavelength(tid::Symbol) = vacuum_wavelength(transitions(tid))


# ====================================================================
abstract type AbstractLine end
abstract type ForbiddenLine <: AbstractLine  end
abstract type PermittedLine <: AbstractLine  end
abstract type NarrowLine    <: PermittedLine end
abstract type BroadLine     <: PermittedLine end
abstract type VeryBroadLine <: BroadLine     end
abstract type NuisanceLine  <: AbstractLine  end


struct LineDescriptor{T}
    id::T
    types::Vector{DataType}
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


replace_descriptor!(ld::Vector{<: LineDescriptor}, id::T, type::DataType) where T = replace_descriptor!(ld, id, [type])
function replace_descriptor!(ld::Vector{<: LineDescriptor}, id::T, types::Vector{DataType}) where T
    for i in 1:length(ld)
        if isa(ld[i], LineDescriptor{T})
            if ld[i].id == id
                ld[i] = LineDescriptor(id, types...)
                return nothing
            end
        end
    end
end


line_suffix(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: AbstractLine})  = ""
line_suffix(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: NarrowLine})    = "_na"
line_suffix(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: BroadLine})     = "_br"
line_suffix(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: VeryBroadLine}) = "_bb"

line_group( recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: ForbiddenLine}) = :NarrowLines
line_group( recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: NarrowLine})    = :NarrowLines
line_group( recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: BroadLine})     = :BroadLines
line_group( recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: VeryBroadLine}) = :VeryBroadLines
line_group( recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: NuisanceLine}) =  :NuisanceLines

line_profile(::Recipe{<: AbstractRecipeSpec}, ::Type{<: AbstractLine}, id::Val) = :gauss
line_profile(::Recipe{<: AbstractRecipeSpec}, ::Type{<: AbstractLine})          = :gauss

line_cname(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{T}, id::Symbol) where T <: AbstractLine = Symbol(   id    , line_suffix(recipe, T))
line_cname(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{T}, λ::Float64) where T <: AbstractLine = Symbol(:l, λ, :A, line_suffix(recipe, T))

function set_constraints!(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: ForbiddenLine}, comp::AbstractSpecLineComp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 2e3
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
end

function set_constraints!(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: NarrowLine}, comp::AbstractSpecLineComp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 1e3 # avoid confusion with the broad component
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
end

function set_constraints!(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: BroadLine}, comp::AbstractSpecLineComp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 900, 5e3, 1.5e4
    comp.voff.low, comp.voff.val, comp.voff.high = -3e3, 0, 3e3
end

function set_constraints!(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: VeryBroadLine}, comp::AbstractSpecLineComp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 1e4, 2e4, 3e4
    comp.voff.fixed = true
end

function set_constraints!(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{<: NuisanceLine}, comp::AbstractSpecLineComp)
    comp.norm.val = 0.
    comp.center.fixed = false;  comp.voff.fixed = true
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 600, 5e3, 1e4
end

line_component(::Val{:gauss}  , λ::Float64) = SpecLineGauss(λ)
line_component(::Val{:lorentz}, λ::Float64) = SpecLineLorentz(λ)
line_component(::Val{:voigt}  , λ::Float64) = SpecLineVoigt(λ)

function line_component(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{T}, id::Symbol) where T <: AbstractLine
    λ = vacuum_wavelength(id)
    profile = line_profile(recipe, T, Val(id))
    comp = line_component(Val(profile), λ)
    set_constraints!(recipe, T, comp)
    return comp
end

function line_component(recipe::Recipe{<: AbstractRecipeSpec}, ::Type{T}, λ::Float64) where T <: AbstractLine
    profile = line_profile(recipe, T)
    comp = line_component(Val(profile), λ)
    set_constraints!(recipe, T, comp)
    return comp
end
