include("ATL.jl")

using .ATL

export ForbiddenLine, SemiForbiddenLine, NarrowLine, BroadLine, VeryBroadLine, NuisanceLine


abstract type LineTemplate end
abstract type ForbiddenLine     <: LineTemplate  end
abstract type SemiForbiddenLine <: LineTemplate  end
abstract type NarrowLine        <: LineTemplate  end
abstract type BroadLine         <: LineTemplate  end
abstract type VeryBroadLine     <: LineTemplate  end
abstract type NuisanceLine      <: LineTemplate  end


default_line_templates(::CRecipe, ::ATL.Transition{<: ATL.Forbidden})     = [ForbiddenLine]
default_line_templates(::CRecipe, ::ATL.Transition{<: ATL.SemiForbidden}) = [SemiForbiddenLine]
default_line_templates(::CRecipe, ::ATL.Transition{<: ATL.Permitted})     = [NarrowLine, BroadLine]

line_suffix(::CRecipe, ::Type{<: LineTemplate})  = ""
line_suffix(::CRecipe, ::Type{<: NarrowLine})    = "_na"
line_suffix(::CRecipe, ::Type{<: BroadLine})     = "_br"
line_suffix(::CRecipe, ::Type{<: VeryBroadLine}) = "_bb"

line_group( ::CRecipe, ::Type{<: ForbiddenLine})     = :NarrowLines
line_group( ::CRecipe, ::Type{<: SemiForbiddenLine}) = :BroadLines
line_group( ::CRecipe, ::Type{<: NarrowLine})        = :NarrowLines
line_group( ::CRecipe, ::Type{<: BroadLine})         = :BroadLines
line_group( ::CRecipe, ::Type{<: VeryBroadLine})     = :VeryBroadLines

line_component(::CRecipe, center::Float64) = SpecLineGauss(center)
line_component(recipe::CRecipe, tid::Val{TID}) where TID = line_component(recipe, get_wavelength(get_transition(TID)))

function line_component(recipe::CRecipe, tid::Union{Val{TID}, Float64}, template::Type{<: ForbiddenLine}) where TID
    @track_recipe
    comp = line_component(recipe, tid)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 2e3
    comp.voff.low, comp.voff.val, comp.voff.high =-5e2,   0, 5e2
    return comp
end

function line_component(recipe::CRecipe, tid::Union{Val{TID}, Float64}, template::Type{<: SemiForbiddenLine}) where TID
    @track_recipe
    return line_component(recipe, tid, BroadLine)
end

function line_component(recipe::CRecipe, tid::Union{Val{TID}, Float64}, template::Type{<: NarrowLine}) where TID
    @track_recipe
    comp = line_component(recipe, tid)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 1e3 # avoid confusion with the broad component
    comp.voff.low, comp.voff.val, comp.voff.high =-5e2,  0, 5e2
    return comp
end

function line_component(recipe::CRecipe, tid::Union{Val{TID}, Float64}, template::Type{<: BroadLine}) where TID
    @track_recipe
    comp = line_component(recipe, tid)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 900, 5e3, 1.5e4
    comp.voff.low, comp.voff.val, comp.voff.high =-5e2,  0, 5e2
    return comp
end

function line_component(recipe::CRecipe, tid::Union{Val{TID}, Float64}, template::Type{<: VeryBroadLine}) where TID
    @track_recipe
    comp = line_component(recipe, tid)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 1e4, 2e4, 3e4
    comp.voff.fixed = true
    return comp
end

function line_component(recipe::CRecipe, wl::Float64, template::Type{<: NuisanceLine})
    @track_recipe
    comp = line_component(recipe, wl)
    comp.center.fixed = false
    comp.voff.fixed = true
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 600, 5e3, 1e4
    return comp
end


# ====================================================================
struct SpecLine{T <: LineTemplate}
    wavelength::Float64
    comp::GModelFit.AbstractComponent
    group::Symbol

    function SpecLine{T}(recipe::CRecipe, tid::Symbol) where T <: LineTemplate
        transition = get_transition(tid)
        wl = get_wavelength(transition)
        comp = line_component(recipe, Val(tid), T)
        return new{T}(wl, comp, line_group(recipe, T))
    end
    function SpecLine{T}(recipe::CRecipe, wl::Float64) where T <: LineTemplate
        comp = line_component(recipe, wl, T)
        return new{T}(wl, comp, line_group(recipe, T))
    end
end

struct SpecLineSet <: AbstractDict{Symbol, SpecLine}
    dict::OrderedDict{Symbol, SpecLine}
    SpecLineSet() = new(OrderedDict{Symbol, SpecLine}())
end

haskey(s::SpecLineSet, k::Symbol) = haskey(s.dict, k)
keys(s::SpecLineSet) = keys(s.dict)
values(s::SpecLineSet) = values(s.dict)
getindex(s::SpecLineSet, k::Symbol) = s.dict[k]
function setindex!(s::SpecLineSet, v::SpecLine, k::Symbol)
    @assert !(k in keys(s.dict)) "Component already exists: $k"
    s.dict[k] = v
end
function iterate(s::SpecLineSet, i=1)
    k = collect(keys(s))
    (i > length(k))  &&  return nothing
    return (k[i] => s[k[i]], i+1)
end
delete!(s::SpecLineSet, k::Symbol) = delete!(s.dict, k)

function sort_by_wavelength!(s::SpecLineSet)
    dict = deepcopy(s.dict)
    ii = sortperm(getfield.(values(dict), :wavelength))
    empty!(s.dict)
    for k in collect(keys(dict))[ii]
        s.dict[k] = dict[k]
    end
    nothing
end

add_line!(recipe::CRecipe, lines::SpecLineSet, tid::Symbol)  = add_line!(recipe, lines, tid, default_line_templates(recipe, get_transition(tid))...)

function add_line!(recipe::CRecipe, lines::SpecLineSet, tid::Symbol, templates...)
    for template in templates
        cname = Symbol(tid, line_suffix(recipe, template))
        lines[cname] = SpecLine{template}(recipe, tid)
    end
end

function add_line!(recipe::CRecipe, lines::SpecLineSet, wl::Float64, templates...)
    for template in templates
        cname = Symbol(@sprintf("l%.1f_", wl), line_suffix(recipe, template))
        lines[cname] = SpecLine{template}(recipe, wl)
    end
end
