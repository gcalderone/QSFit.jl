include("ATL.jl")

using .ATL

export ForbiddenLine, SemiForbiddenLine, NarrowLine, BroadLine, VeryBroadLine, NuisanceLine
export use_line!


abstract type LineTemplate end
abstract type ForbiddenLine     <: LineTemplate  end
abstract type SemiForbiddenLine <: LineTemplate  end
abstract type NarrowLine        <: LineTemplate  end
abstract type BroadLine         <: LineTemplate  end
abstract type VeryBroadLine     <: LineTemplate  end
abstract type NuisanceLine      <: LineTemplate  end


default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.Forbidden})     where N = (ForbiddenLine)
default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.SemiForbidden}) where N = (SemiForbiddenLine)
default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.Permitted})     where N = (NarrowLine, BroadLine)
default_line_templates(::CRecipe, ::ATL.UnidentifiedTransition)                      = (NuisanceLine)

line_suffix(::CRecipe, ::Type{<: LineTemplate})  = ""
line_suffix(::CRecipe, ::Type{<: NarrowLine})    = "_na"
line_suffix(::CRecipe, ::Type{<: BroadLine})     = "_br"
line_suffix(::CRecipe, ::Type{<: VeryBroadLine}) = "_bb"

line_group( ::CRecipe, ::Type{<: ForbiddenLine})     = :NarrowLines
line_group( ::CRecipe, ::Type{<: SemiForbiddenLine}) = :BroadLines
line_group( ::CRecipe, ::Type{<: NarrowLine})        = :NarrowLines
line_group( ::CRecipe, ::Type{<: BroadLine})         = :BroadLines
line_group( ::CRecipe, ::Type{<: VeryBroadLine})     = :VeryBroadLines
line_group( ::CRecipe, ::Type{<: NuisanceLine})      = :NuisanceLines

line_component(::CRecipe, center::Float64) = SpecLineGauss(center)

function line_component(recipe::CRecipe, center::Float64, template::Type{<: ForbiddenLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 2e3
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: SemiForbiddenLine})
    @track_recipe
    return line_component(recipe, center, BroadLine)
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: NarrowLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 1e3 # avoid confusion with the broad component
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: BroadLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 900, 5e3, 1.5e4
    comp.voff.low, comp.voff.val, comp.voff.high = -3e3, 0, 3e3
    return comp
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: VeryBroadLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 1e4, 2e4, 3e4
    comp.voff.fixed = true
    return comp
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: NuisanceLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.voff.fixed = true
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 600, 5e3, 1e4
    return comp
end


# ====================================================================
struct SpectralLine{T <: LineTemplate}
    tid::Symbol
    group::Symbol
    comp::AbstractSpecLineComp
    wavelength::Float64
end

function show(io::IO, line::SpectralLine)
    println(io, line.tid, " (", string(typeof(line.t)))
    show(io, line.comp)
end

# ====================================================================
struct LineSet <: AbstractDict{Symbol, SpectralLine}
    dict::OrderedDict{Symbol, SpectralLine}
    LineSet() = new(OrderedDict{Symbol, SpectralLine}())
    LineSet(cname::Symbol, line::SpectralLine) = new(OrderedDict{Symbol, SpectralLine}(cname => line))
end

haskey(s::LineSet, k::Symbol) = haskey(s.dict, k)
keys(s::LineSet) = keys(s.dict)
getindex(s::LineSet, k::Symbol) = s.dict[k]
function setindex!(s::LineSet, v::SpectralLine, k::Symbol)
    @assert !(k in keys(s.dict)) "Component already exists: $k"
    s.dict[k] = v
end
function iterate(s::LineSet, i=1)
    k = collect(keys(s))
    (i > length(k))  &&  return nothing
    return (k[i] => s[k[i]], i+1)
end
function merge!(s::LineSet, n::LineSet)
    for (cname, line) in n
        s[cname] = line
    end
end
delete!(s::LineSet, k::Symbol) = delete!(s.dict, k)

function sort_by_wavelength!(s::LineSet)
    dict = deepcopy(s.dict)
    ii = sortperm(getfield.(values(dict), :wavelength))
    empty!(s.dict)
    for k in collect(keys(dict))[ii]
        s.dict[k] = dict[k]
    end
    nothing
end


LineSet(recipe::CRecipe, tid::Symbol)             = LineSet(recipe, Val(tid))
LineSet(recipe::CRecipe, val::Val{tid}) where tid = LineSet(recipe, get_transition(tid), default_line_templates(recipe, get_transition(tid)))
LineSet(recipe::CRecipe, center::Float64)         = LineSet(recipe, ATL.UnidentifiedTransition(center), default_line_templates(recipe, ATL.UnidentifiedTransition(center)))

LineSet(recipe::CRecipe, tid::Symbol    , templates...)             = LineSet(recipe, Val(tid), templates...)
LineSet(recipe::CRecipe,    ::Val{tid}  , templates...) where {tid} = LineSet(recipe, get_transition(tid), templates)
LineSet(recipe::CRecipe, center::Float64, templates...)             = LineSet(recipe, ATL.UnidentifiedTransition(center), templates)

function LineSet(recipe::CRecipe, t::ATL.AbstractTransition, templates::Tuple)
    out = LineSet()
    for template in templates
        merge!(out, LineSet(recipe, t, template))
    end
    return out
end

# Singlets
function LineSet(recipe::CRecipe, t::ATL.Transition{1,T}, template::Type{<: LineTemplate}) where T
    @track_recipe
    tid = get_id(t)
    cname = Symbol(tid, line_suffix(recipe, template))
    group = line_group(recipe, template)
    wl = ATL.get_wavelengths(t)[1]
    comp = line_component(recipe, wl, template)
    return LineSet(cname, SpectralLine{template}(tid, group, comp, wl))
end

# Handle multiplets
function LineSet(recipe::CRecipe, t::ATL.Transition{N,T}, template::Type{<: LineTemplate}) where {N,T}
    @track_recipe
    tid = get_id(t)
    cname = Symbol(tid, line_suffix(recipe, template))
    group = line_group(recipe, template)
    wl = mean(ATL.get_wavelengths(t)) # using average wavelength for multiplet
    comp = line_component(recipe, wl, template)
    return LineSet(cname, SpectralLine{template}(tid, group, comp, wl))
end

# Unidentified lines
function LineSet(recipe::CRecipe, t::ATL.UnidentifiedTransition, template::Type{<: LineTemplate})
    @track_recipe
    tid = get_id(t)
    cname = Symbol(tid, line_suffix(recipe, template))
    group = line_group(recipe, template)
    wl = ATL.get_wavelengths(t)[1]
    comp = line_component(recipe, wl, template)
    return LineSet(cname, SpectralLine{template}(tid, group, comp, wl))
end
