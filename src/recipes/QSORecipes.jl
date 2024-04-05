module QSORecipes

using Printf, DataStructures, Statistics
using Dierckx
using ..QSFit, GModelFit

import QSFit: init_recipe!, preprocess_spec!, line_profile, line_suffix, set_constraints!, analyze, reduce

abstract type QSOGeneric <: AbstractRecipeSpec end

line_profile(recipe::Recipe{<: QSOGeneric}, ::Type{<: AbstractLine}, id::Val) = recipe.line_profiles
line_profile(recipe::Recipe{<: QSOGeneric}, ::Type{<: AbstractLine})          = recipe.line_profiles


function init_recipe!(recipe::Recipe{T}) where T <: QSOGeneric
    @invoke init_recipe!(recipe::Recipe{<: AbstractRecipeSpec})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.min_spectral_coverage = Dict(:default => 0.6)

    recipe.host_template = Dict(:library=>"swire", :template=>"Ell5")
    recipe.host_template_ref_wavelength = 5500. # A
    recipe.use_host_template = true
    recipe.host_template_range = [4000., 7000.]

    recipe.line_profiles = :gauss

    recipe.n_nuisance = 10
    recipe.nuisance_avoid = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]  # Angstrom
    recipe.nuisance_maxoffset_from_guess = 1e3  # km/s

    recipe.lines = LineDescriptor[]
end


struct LineInstance{T}
    id::T
    type::DataType
    group::Symbol
    comp::GModelFit.AbstractComponent
end


function dict_line_instances(recipe::Recipe{<: QSOGeneric}, linedescrs::Vector{LineDescriptor})
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


function fit!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
    resid.mzer.config.ftol = 1.e-6
    bestfit, stats = GModelFit.minimize!(resid)
    show(stats)
    # println(resid.mzer.result)
    return bestfit, stats
end


function add_qso_continuum!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function add_host_galaxy!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function renorm_cont!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function guess_norm_factor!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals, name::Symbol; quantile=0.95)
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


function add_emission_lines!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function add_patch_functs!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function add_nuisance_lines!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function neglect_weak_features!(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


function preprocess_spec!(recipe::Recipe{T}, spec::QSFit.Spectrum) where T <: QSOGeneric
    @invoke preprocess_spec!(recipe::Recipe{<: AbstractRecipeSpec}, spec)

    spec.good[findall(spec.x .< recipe.wavelength_range[1])] .= false
    spec.good[findall(spec.x .> recipe.wavelength_range[2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is
    not sufficient the component should not be added to the model,
    and corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println("Good samples before line coverage filter: ", count(spec.good) , " / ", length(spec.good))

    # Collect line components (neglecting the ones with insufficent coverage)
    lcs = dict_line_instances(recipe, recipe.lines)
    for loop in 1:2
        # The second pass is required to neglect lines whose coverage
        # has been affected by the neglected spectral samples.
        if loop == 2
            println()
            println("Updated coverage:")
        end
        for (cname, line) in lcs
            threshold = get(recipe.min_spectral_coverage, cname, recipe.min_spectral_coverage[:default])
            (λmin, λmax, coverage) = QSFit.spectral_coverage(spec.x[findall(spec.good)],
                                                             spec.resolution, line.comp)
            @printf("Line %-15s coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax)
            if coverage < threshold
            @printf(", threshold is < %5.3f, neglecting...", threshold)
                spec.good[λmin .<= spec.x .< λmax] .= false
                delete!(lcs, cname)
            end
            println()
        end
    end
    println("Good samples after line coverage filter: ", count(spec.good) , " / ", length(spec.good))

    # Sort lines according to center wavelength, and save list in recipe
    kk = collect(keys(lcs))
    vv = collect(values(lcs))
    ii = sortperm(getfield.(getfield.(getfield.(vv, :comp), :center), :val))
    recipe.lcs = OrderedDict(Pair.(kk[ii], vv[ii]))
end


function reduce(recipe::Recipe{<: QSOGeneric}, resid::GModelFit.Residuals)
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


include("QSORecipes_Type1.jl")
include("QSORecipes_Type2.jl")

end
