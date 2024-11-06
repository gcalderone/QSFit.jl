include("ATL.jl")

using .ATL

export ForbiddenLine, SemiForbiddenLine, NarrowLine, BroadLine, VeryBroadLine, NuisanceLine
export get_lines_dict, set_lines_dict!, add_line!


abstract type LineTemplate end
abstract type ForbiddenLine     <: LineTemplate  end
abstract type SemiForbiddenLine <: LineTemplate  end
abstract type NarrowLine        <: LineTemplate  end
abstract type BroadLine         <: LineTemplate  end
abstract type VeryBroadLine     <: LineTemplate  end
abstract type NuisanceLine      <: LineTemplate  end

default_line_templates(::CRecipe, ::ATL.AbstractTransition{<: ATL.Forbidden})     = [ForbiddenLine]
default_line_templates(::CRecipe, ::ATL.AbstractTransition{<: ATL.SemiForbidden}) = [SemiForbiddenLine]
default_line_templates(::CRecipe, ::ATL.AbstractTransition{<: ATL.Permitted})     = [NarrowLine, BroadLine]
default_line_templates(::CRecipe, ::ATL.UnidentifiedTransition)                   = [NuisanceLine]

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
line_component(recipe::CRecipe, t::ATL.AbstractTransition, template::Type{<: LineTemplate}) = line_component(recipe, get_wavelength(t), template)

function line_component(recipe::CRecipe, center::Float64, template::Type{<: ForbiddenLine})
    @track_recipe
    comp = line_component(recipe, center)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 2e3
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end

function line_component(recipe::CRecipe, center::Float64, template::Type{<: SemiForbiddenLine})
    @track_recipe
    return line_component(recipe, BroadLine, t)
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

add_line!(recipe::CRecipe, tid::Symbol) = add_line!(recipe, get_transition(tid))
add_line!(recipe::CRecipe, tid::Symbol, templates::Vararg{DataType, N}) where N = add_line!(recipe, get_transition(tid), templates...)

add_line!(recipe::CRecipe, center::Float64) = add_line!(recipe, ATL.UnidentifiedTransition(center))
add_line!(recipe::CRecipe, center::Float64, templates::Vararg{DataType, N}) where N = add_line!(recipe, ATL.UnidentifiedTransition(center), templates...)

add_line!(recipe::CRecipe, t::ATL.AbstractTransition) = add_line!(recipe, t, default_line_templates(recipe, t)...)
function add_line!(recipe::CRecipe, t::ATL.AbstractTransition, templates::Vararg{DataType, N}) where N
    @track_recipe
    @assert N > 0
    dict = get_lines_dict(recipe)
    for template in templates
        @assert template <: LineTemplate
        cname = Symbol(get_id(t), line_suffix(recipe, template))
        @assert !(cname in keys(dict)) "Duplicated component name: $cname"
        group = line_group(recipe, template)
        comp = line_component(recipe, t, template)
        dict[cname] = SpectralLine{template}(get_id(t), t, group, comp)
    end
    return dict
end

function get_lines_dict(recipe::CRecipe)
    if !(:lines in propertynames(recipe))
        @track_recipe
        recipe.lines = OrderedDict{Symbol, SpectralLine}()
    end
    return recipe.lines
end

function set_lines_dict!(recipe::CRecipe)
    @track_recipe
    recipe.lines = OrderedDict{Symbol, SpectralLine}()
    for tid in ATL.get_transition_ids()
        add_line!(recipe, tid)
    end
    return get_lines_dict(recipe)
end
