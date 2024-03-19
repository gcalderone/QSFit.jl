function preprocess_spec!(recipe::Recipe{T}, spec::QSFit.Spectrum) where T <: Type1
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


function analyze(recipe::Recipe{<: Type1}, spec::Spectrum, resid::GModelFit.Residuals)
    recipe.spec = spec  # TODO: remove
    model = resid.meval.model
    model[:main] = SumReducer()
    select_maincomp!(model, :main)

    println("\nFit continuum components...")
    model[:Continuum] = SumReducer()
    push!(model[:main].list, :Continuum)
    add_qso_continuum!(recipe, resid)
    add_host_galaxy!(recipe, resid)
    add_balmer_cont!(recipe, resid)
    fit!(recipe, resid)
    renorm_cont!(recipe, resid)
    freeze!(model, :QSOcont)
    haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
    haskey(model, :Balmer)  &&  freeze!(model, :Balmer)
    GModelFit.update!(resid.meval)

    println("\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!(recipe, resid)
    add_iron_opt!(recipe, resid)

    if length(model[:Iron].list) > 0
        fit!(recipe, resid)
        haskey(model, :Ironuv   )  &&  freeze!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  freeze!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  freeze!(model, :Ironoptna)
    end
    GModelFit.update!(resid.meval)

    println("\nFit known emission lines...")
    add_emission_lines!(recipe, resid)
    add_patch_functs!(recipe, resid)

    fit!(recipe, resid)
    for cname in keys(recipe.lcs)
        freeze!(model, cname)
    end

    println("\nFit nuisance emission lines...")
    add_nuisance_lines!(recipe, resid)

    println("\nLast run with all parameters free...")
    thaw!(model, :QSOcont)
    haskey(model, :Galaxy   )  &&  thaw!(model, :Galaxy)
    haskey(model, :Balmer   )  &&  thaw!(model, :Balmer)
    haskey(model, :Ironuv   )  &&  thaw!(model, :Ironuv)
    haskey(model, :Ironoptbr)  &&  thaw!(model, :Ironoptbr)
    haskey(model, :Ironoptna)  &&  thaw!(model, :Ironoptna)
    for cname in keys(recipe.lcs)
        thaw!(model, cname)
    end
    for j in 1:recipe.n_nuisance
        cname = Symbol(:nuisance, j)
        if model[cname].norm.val > 0
            thaw!(model, cname)
        else
            freeze!(model, cname)
        end
    end
    bestfit, stats = fit!(recipe, resid)

    if neglect_weak_features!(recipe, resid)
        println("\nRe-run fit...")
        bestfit, stats = fit!(recipe, resid)
    end

    return bestfit, stats
end
