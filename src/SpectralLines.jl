include("ATL.jl")

using .ATL

export ForbiddenLine, SemiForbiddenLine, NarrowLine, BroadLine, VeryBroadLine, NuisanceLine
export add_line!


abstract type LineTemplate end
abstract type ForbiddenLine     <: LineTemplate  end
abstract type SemiForbiddenLine <: LineTemplate  end
abstract type NarrowLine        <: LineTemplate  end
abstract type BroadLine         <: LineTemplate  end
abstract type VeryBroadLine     <: LineTemplate  end
abstract type NuisanceLine      <: LineTemplate  end


default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.Forbidden})     where N = [ForbiddenLine]
default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.SemiForbidden}) where N = [SemiForbiddenLine]
default_line_templates(::CRecipe, ::ATL.Transition{N, <: ATL.Permitted})     where N = [NarrowLine, BroadLine]
default_line_templates(::CRecipe, ::ATL.UnidentifiedTransition)                      = [NuisanceLine]

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
    t::ATL.AbstractTransition
    group::Symbol
    comp::AbstractSpecLineComp
end

function show(io::IO, line::SpectralLine)
    println(io, line.tid, " (", string(typeof(line.t)))
    show(io, line.comp)
end

add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, val::Val{tid}) where {tid}                                    = [add_line!(recipe, dict, val   , t) for t in default_line_templates(recipe, get_transition(tid))]
add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine},    ::Val{tid}  , template::Type{<: LineTemplate}) where {tid} =  add_line!(recipe, dict, get_transition(tid), template)

add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, tid::Symbol)                                                  = [add_line!(recipe, dict, tid   , t) for t in default_line_templates(recipe, get_transition(tid))]
add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, tid::Symbol    , template::Type{<: LineTemplate})             =  add_line!(recipe, dict, get_transition(tid), template)

add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, center::Float64)                                              = [add_line!(recipe, dict, center, t) for t in default_line_templates(recipe, ATL.UnidentifiedTransition(center))]
add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, center::Float64, template::Type{<: LineTemplate})             =  add_line!(recipe, dict, ATL.UnidentifiedTransition(center), template)

# Singlets
function add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, t::ATL.Transition{1,T}, template::Type{<: LineTemplate}) where T
    @track_recipe
    tid = get_id(t)
    cname = Symbol(tid, line_suffix(recipe, template))
    @assert !(cname in keys(dict)) "Duplicated component name: $cname"
    group = line_group(recipe, template)
    comp = line_component(recipe, ATL.get_wavelength(t), template)
    dict[cname] = SpectralLine{template}(tid, t, group, comp)
    nothing
end

# Handle multiplets
function add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, t::ATL.Transition{N,T}, template::Type{<: LineTemplate}) where {N,T}
    @track_recipe
    tid = get_id(t)
    cname = Symbol(tid, line_suffix(recipe, template))
    @assert !(cname in keys(dict)) "Duplicated component name: $cname"
    group = line_group(recipe, template)
    comp = line_component(recipe, ATL.get_wavelength(t), template)
    dict[cname] = SpectralLine{template}(tid, t, group, comp)
    nothing
end


# Unidentified lines
function add_line!(recipe::CRecipe, dict::OrderedDict{Symbol, SpectralLine}, t::ATL.UnidentifiedTransition, template::Type{<: LineTemplate}) where T
    @track_recipe
    tid = get_id(t)
    cname = tid
    @assert !(cname in keys(dict)) "Duplicated component name: $cname"
    group = line_group(recipe, template)
    comp = line_component(recipe, ATL.get_wavelength(t), template)
    dict[cname] = SpectralLine{template}(tid, t, group, comp)
    nothing
end
