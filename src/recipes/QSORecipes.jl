module QSORecipes

using Printf, DataStructures, Statistics
using Dierckx
using ..QSFit, GModelFit

import QSFit: init_recipe!, preprocess_spec!, line_profile, line_suffix, set_constraints!, analyze, reduce

export Type1

abstract type Type1 <: AbstractRecipeSpec end

line_profile(recipe::Recipe{<: Type1}, ::Type{<: AbstractLine}, id::Val) = recipe.line_profiles
line_profile(recipe::Recipe{<: Type1}, ::Type{<: AbstractLine})          = recipe.line_profiles

# Special cases for emission lines
abstract type MgIIBroadLine <: BroadLine end
function set_constraints!(recipe::Recipe{<: Type1}, ::Type{MgIIBroadLine}, comp::GModelFit.AbstractComponent)
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
end


abstract type OIIIBlueWing  <: NarrowLine end
line_suffix(recipe::Recipe{<: Type1}, ::Type{OIIIBlueWing}) = :_bw
function set_constraints!(recipe::Recipe{<: Type1}, ::Type{OIIIBlueWing}, comp::GModelFit.AbstractComponent)
    comp.voff.low, comp.voff.val, comp.voff.high = 0, 0, 2e3
end

function init_recipe!(recipe::Recipe{T}) where T <: Type1
    @invoke init_recipe!(recipe::Recipe{<: AbstractRecipeSpec})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.min_spectral_coverage = Dict(:default => 0.6,
                                        :Ironuv  => 0.3,
                                        :Ironopt => 0.3)

    recipe.host_template = Dict(:library=>"swire", :template=>"Ell5")
    recipe.host_template_ref_wavelength = 5500. # A
    recipe.use_host_template = true
    recipe.host_template_range = [4000., 7000.]

    recipe.use_balmer = true
    recipe.use_ironuv = true;      recipe.Ironuv_fwhm    = 3000.
    recipe.use_ironopt = true;     recipe.Ironoptbr_fwhm = 3000.;  recipe.Ironoptna_fwhm =  500.

    recipe.line_profiles = :gauss

    recipe.n_nuisance = 10
    recipe.nuisance_avoid = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]  # Angstrom
    recipe.nuisance_maxoffset_from_guess = 1e3  # km/s

    recipe.lines = [
        LineDescriptor(:Lyb        , NarrowLine, BroadLine),
        # LineDescriptor(:OV_1213   , ForbiddenLine)  # 1213.8A, Ferland+92, Shields+95,
        LineDescriptor(:Lya        , NarrowLine, BroadLine),
        # LineDescriptor(:OV_1218   , ForbiddenLine)  # 1218.3A, Ferland+92, Shields+95,
        LineDescriptor(:NV_1241    , ForbiddenLine),
        LineDescriptor(:OI_1306    , BroadLine),
        LineDescriptor(:CII_1335   , BroadLine),
        LineDescriptor(:SiIV_1400  , BroadLine),
        LineDescriptor(:CIV_1549   , NarrowLine, BroadLine),
        LineDescriptor(:HeII_1640  , BroadLine),
        LineDescriptor(:OIII_1664  , BroadLine),
        LineDescriptor(:AlIII_1858 , BroadLine),
        LineDescriptor(:CIII_1909  , BroadLine),
        LineDescriptor(:CII_2326   , BroadLine),
        LineDescriptor(2420.0      , BroadLine),
        LineDescriptor(:MgII_2798  , NarrowLine, MgIIBroadLine),
        LineDescriptor(:NeV_3345   , ForbiddenLine),
        LineDescriptor(:NeV_3426   , ForbiddenLine),
        LineDescriptor(:OII_3727   , ForbiddenLine),
        LineDescriptor(:NeIII_3869 , ForbiddenLine),
        LineDescriptor(:Hd         , BroadLine),
        LineDescriptor(:Hg         , BroadLine),
        LineDescriptor(:OIII_4363  , ForbiddenLine),
        LineDescriptor(:HeII_4686  , BroadLine),
        LineDescriptor(:Hb         , NarrowLine, BroadLine),
        LineDescriptor(:OIII_4959  , ForbiddenLine),
        LineDescriptor(:OIII_5007  , ForbiddenLine, OIIIBlueWing),
        LineDescriptor(:HeI_5876   , BroadLine),
        LineDescriptor(:OI_6300    , ForbiddenLine),
        LineDescriptor(:OI_6364    , ForbiddenLine),
        LineDescriptor(:NII_6549   , ForbiddenLine),
        LineDescriptor(:Ha         , NarrowLine, BroadLine, VeryBroadLine),
        LineDescriptor(:NII_6583   , ForbiddenLine),
        LineDescriptor(:SII_6716   , ForbiddenLine),
        LineDescriptor(:SII_6731   , ForbiddenLine)
    ]
end


struct LineInstance{T}
    id::T
    type::DataType
    group::Symbol
    comp::GModelFit.AbstractComponent
end


function dict_line_instances(recipe::Recipe{<: Type1}, linedescrs::Vector{LineDescriptor})
    out = OrderedDict{Symbol, LineInstance}()
    for ld in linedescrs
        for T in ld.types
            cname = line_cname(recipe, T, ld.id)
            @assert !(cname in keys(out)) "Duplicate component name: $cname"
            out[cname] = LineInstance(ld.id, T, line_group(recipe, T), line_component(recipe, T, ld.id))
        end
    end
    return out
end


function fit!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    resid.mzer.config.ftol = 1.e-6
    bestfit, stats = GModelFit.minimize!(resid)
    show(stats)
    # println(resid.mzer.result)
    return bestfit, stats
end


function add_qso_continuum!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)

    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(resid.data)) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    resid.meval.model[:QSOcont] = comp
    push!(resid.meval.model[:Continuum].list, :QSOcont)
    GModelFit.update!(resid.meval)
end


function add_host_galaxy!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)
    if recipe.use_host_template  &&
        (recipe.host_template_range[1] .< maximum(λ))  &&
        (recipe.host_template_range[2] .> minimum(λ))
        resid.meval.model[:Galaxy] = QSFit.hostgalaxy(recipe.host_template[:template],
                                                      library=recipe.host_template[:library],
                                                      refwl=recipe.host_template_ref_wavelength)
        push!(resid.meval.model[:Continuum].list, :Galaxy)

        # Split total flux between continuum and host galaxy
        refwl = recipe.host_template_ref_wavelength
        vv = Dierckx.Spline1D(λ, values(resid.data), k=1, bc="extrapolate")(refwl)
        @assert !isnan(vv) "Predicted L_λ at $(refwl)A is NaN"
        if vv <= 0
            @warn "Predicted L_λ at $(refwl)A is negative, set host galaxy guess value at zero."
            resid.meval.model[:Galaxy].norm.val   = 0.
        else
            resid.meval.model[:Galaxy].norm.val   = 1/2 * vv
            resid.meval.model[:QSOcont].norm.val *= 1/2 * vv / Dierckx.Spline1D(λ, GModelFit.last_evaluation(resid.meval, :QSOcont), k=1, bc="extrapolate")(refwl)
        end
        GModelFit.update!(resid.meval)
    end
end


function add_balmer_cont!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    if recipe.use_balmer
        resid.meval.model[:Balmer] = QSFit.balmercont(0.1, 0.5)
        push!(resid.meval.model[:Continuum].list, :Balmer)
        c = resid.meval.model[:Balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        resid.meval.model[:Balmer].norm.patch = @fd (m, v) -> v * m[:QSOcont].norm
        GModelFit.update!(resid.meval)
    end
end


function renorm_cont!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    freeze!(resid.meval.model, :QSOcont)
    c = resid.meval.model[:QSOcont]
    d = resid.meval.cevals[:QSOcont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println("Cont. norm. (before): ", c.norm.val)
        scaling = 0.99
        while c.norm.val * scaling > c.norm.low
            residuals = (GModelFit.last_evaluation(resid.meval) - values(resid.data)) ./ uncerts(resid.data)
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= scaling
            GModelFit.update!(resid.meval)
        end
        println("Cont. norm. (after) : ", c.norm.val)
    else
        println("Skipping cont. renormalization")
    end
    GModelFit.update!(resid.meval)
end


function guess_norm_factor!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals, name::Symbol; quantile=0.95)
    @assert resid.meval.model[name].norm.val != 0
    m = GModelFit.last_evaluation(resid.meval, name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    r = values(resid.data) - GModelFit.last_evaluation(resid.meval)
    ratio = resid.meval.model[name].norm.val / sum(m[i1:i2])
    off = sum(r[i1:i2]) * ratio
    resid.meval.model[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if resid.meval.model[name].norm.val < 0  # ensure line has positive normalization
        resid.meval.model[name].norm.val = abs(off)
    end
end


function add_iron_uv!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)
    if recipe.use_ironuv
        fwhm = recipe.Ironuv_fwhm
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironuv, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            resid.meval.model[:Ironuv] = comp
            resid.meval.model[:Ironuv].norm.val = 1.
            push!(resid.meval.model[:Iron].list, :Ironuv)
            GModelFit.update!(resid.meval)
            guess_norm_factor!(recipe, resid, :Ironuv)
            GModelFit.update!(resid.meval)
        else
            println("Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)
    if recipe.use_ironopt
        fwhm = recipe.Ironoptbr_fwhm
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironopt, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            resid.meval.model[:Ironoptbr] = comp
            resid.meval.model[:Ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.Ironoptna_fwhm
            resid.meval.model[:Ironoptna] = QSFit.ironopt_narrow(fwhm)
            resid.meval.model[:Ironoptna].norm.val = 1 # TODO: guess a sensible value
            resid.meval.model[:Ironoptna].norm.fixed = false
            push!(resid.meval.model[:Iron].list, :Ironoptbr)
            push!(resid.meval.model[:Iron].list, :Ironoptna)
            GModelFit.update!(resid.meval)
            guess_norm_factor!(recipe, resid, :Ironoptbr)
            GModelFit.update!(resid.meval)
        else
            println("Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    groups = OrderedDict{Symbol, Vector{Symbol}}()

    # Create model components
    for (cname, line) in recipe.lcs
        resid.meval.model[cname] = line.comp
        haskey(groups, line.group)  ||  (groups[line.group] = Vector{Symbol}())
        push!(groups[line.group], cname)
    end

    # Create reducers for groups
    for (group, lnames) in groups
        @assert group in [:BroadLines, :NarrowLines, :VeryBroadLines] "Unexpected group for emission lines: $group"
        resid.meval.model[group] = SumReducer(lnames)
        push!(resid.meval.model[:main].list, group)
    end
    GModelFit.update!(resid.meval)

    # Guess normalizations
    for group in [:BroadLines, :NarrowLines, :VeryBroadLines]  # Note: order is important
        if group in keys(groups)
            for cname in groups[group]
                guess_norm_factor!(recipe, resid, cname)
            end
        end
    end
end


function add_patch_functs!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    model = resid.meval.model
    # Patch parameters
    if haskey(model, :OIII_4959)  &&  haskey(model, :OIII_5007)
        # model[:OIII_4959].norm.patch = @fd m -> m[:OIII_5007].norm / 3
        model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(model, :OIII_5007)  &&  haskey(model, :OIII_5007_bw)
        model[:OIII_5007_bw].voff.patch = @fd (m, v) -> v + m[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm.patch = @fd (m, v) -> v + m[:OIII_5007].fwhm
    end
    if haskey(model, :OI_6300)  &&  haskey(model, :OI_6364)
        # model[:OI_6300].norm.patch = @fd m -> m[:OI_6364].norm / 3
        model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(model, :NII_6549)  &&  haskey(model, :NII_6583)
        # model[:NII_6549].norm.patch = @fd m -> m[:NII_6583].norm / 3
        model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(model, :SII_6716)  &&  haskey(model, :SII_6731)
        # model[:SII_6716].norm.patch = @fd m -> m[:SII_6731].norm / 3
        model[:SII_6716].voff.patch = :SII_6731
    end

    if haskey(model, :Hb_na)  &&  haskey(model, :Ha_na)
        model[:Hb_na].voff.patch = :Ha_na
    end

    # The following are required to avoid degeneracy with iron
    # template
    if haskey(model, :Hg)  &&  haskey(model, :Hb_br)
        model[:Hg].voff.patch = :Hb_br
        model[:Hg].fwhm.patch = :Hb_br
    end
    if haskey(model, :Hg_br)  &&  haskey(model, :Hb_br)
        model[:Hg_br].voff.patch = :Hb_br
        model[:Hg_br].fwhm.patch = :Hb_br
    end
    if haskey(model, :Hg_na)  &&  haskey(model, :Hb_na)
        model[:Hg_na].voff.patch = :Hb_na
        model[:Hg_na].fwhm.patch = :Hb_na
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(model, :Hb_br)  &&
        haskey(model, :Hb_bb)
        model[:Hb_bb].norm.high = 1
        model[:Hb_bb].norm.val  = 0.5
        model[:Hb_bb].norm.patch = @fd (m, v) -> v * m[:Hb_br].norm / m[:Hb_br].fwhm * m[:Hb_bb].fwhm
    end
    if  haskey(model, :Ha_br)  &&
        haskey(model, :Ha_bb)
        model[:Ha_bb].norm.high = 1
        model[:Ha_bb].norm.val  = 0.5
        model[:Ha_bb].norm.patch = @fd (m, v) -> v * m[:Ha_br].norm / m[:Ha_br].fwhm * m[:Ha_bb].fwhm
    end
end


function add_nuisance_lines!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    (recipe.n_nuisance > 0)  ||  (return nothing)

    # Prepare nuisance line components
    for i in 1:recipe.n_nuisance
        resid.meval.model[Symbol(:nuisance, i)] = line_component(recipe, NuisanceLine, 3000.)
    end
    resid.meval.model[line_group(recipe, NuisanceLine)] = SumReducer([Symbol(:nuisance, i) for i in 1:recipe.n_nuisance])
    push!(resid.meval.model[:main].list, line_group(recipe, NuisanceLine))
    GModelFit.update!(resid.meval)
    for j in 1:recipe.n_nuisance
        freeze!(resid.meval.model, Symbol(:nuisance, j))
    end
    GModelFit.update!(resid.meval)

    # Set "nuisance" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λ = coords(resid.meval.domain)
    λnuisance = Vector{Float64}()
    while true
        (length(λnuisance) >= recipe.n_nuisance)  &&  break
        GModelFit.update!(resid.meval)
        Δ = (values(resid.data) - GModelFit.last_evaluation(resid.meval)) ./ uncerts(resid.data)

        # Avoid considering again the same region (within 1A) TODO: within resolution
        for l in λnuisance
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        for rr in recipe.nuisance_avoid
            Δ[findall(rr[1] .< λ .< rr[2])] .= 0.
        end

        # Do not add lines close to from the edges since these may
        # affect qso_cont fitting
        Δ[findall((λ .< minimum(λ)*1.02)  .|
                  (λ .> maximum(λ)*0.98))] .= 0.
        iadd = argmax(Δ)
        (Δ[iadd] <= 0)  &&  break  # No residual is greater than 0, skip further residuals....
        push!(λnuisance, λ[iadd])

        cname = Symbol(:nuisance, length(λnuisance))
        resid.meval.model[cname].norm.val = 1.
        resid.meval.model[cname].center.val  = λ[iadd]

        # Allow to shift by a quantity equal to ...
        @assert recipe.nuisance_maxoffset_from_guess > 0
        resid.meval.model[cname].center.low  = λ[iadd] * (1 - recipe.nuisance_maxoffset_from_guess / 3e5)
        resid.meval.model[cname].center.high = λ[iadd] * (1 + recipe.nuisance_maxoffset_from_guess / 3e5)

        # In any case, we must stay out of avoidance regions
        for rr in recipe.nuisance_avoid
            @assert !(rr[1] .< λ[iadd] .< rr[2])
            if rr[1] .>= λ[iadd]
                resid.meval.model[cname].center.high = min(resid.meval.model[cname].center.high, rr[1])
            end
            if rr[2] .<= λ[iadd]
                resid.meval.model[cname].center.low  = max(resid.meval.model[cname].center.low, rr[2])
            end
        end

        thaw!(resid.meval.model, cname)
        fit!(recipe, resid)
        freeze!(resid.meval.model, cname)
    end
    GModelFit.update!(resid.meval)
end


function neglect_weak_features!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    # Disable "nuisance" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:recipe.n_nuisance
        cname = Symbol(:nuisance, ii)
        isfreezed(resid.meval.model, cname)  &&  continue
        if resid.meval.model[cname].norm.val == 0.
            freeze!(resid.meval.model, cname)
            needs_fitting = true
            println("Disabling $cname (norm. = 0)")
        elseif resid.meval.model[cname].norm.unc / resid.meval.model[cname].norm.val > 3
            resid.meval.model[cname].norm.val = 0.
            freeze!(resid.meval.model, cname)
            needs_fitting = true
            println("Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


function reduce(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(GModelFit.last_evaluation(resid.meval))
    for cname in keys(recipe.lcs)
        haskey(resid.meval.model, cname) || continue
        cont .-= GModelFit.last_evaluation(resid.meval, cname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for cname in keys(recipe.lcs)
        haskey(resid.meval.model, cname) || continue
        EW[cname] = QSFit.int_tabulated(coords(resid.meval.domain),
                                        GModelFit.last_evaluation(resid.meval, cname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end


include("QSORecipes_single.jl")
# TODO include("QSORecipes_multi.jl")
end
