
function analyze(recipe::RRef{<: Type1Recipe}, state::QSFit.State)
    model = Model(:main => SumReducer())
    state.meval = GModelFit.ModelEval(model, domain(state.data))

    println("\nFit continuum components...")
    model[:Continuum] = SumReducer()
    push!(model[:main].list, :Continuum)
    add_qso_continuum!(recipe, state)
    add_host_galaxy!(recipe, state)
    add_balmer_cont!(recipe, state)
    fit!(recipe, state)
    renorm_cont!(recipe, state)
    freeze!(model, :QSOcont)
    haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
    haskey(model, :Balmer)  &&  freeze!(model, :Balmer)
    GModelFit.update!(state.meval)

    println("\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!(recipe, state)
    add_iron_opt!(recipe, state)

    if length(model[:Iron].list) > 0
        fit!(recipe, state)
        haskey(model, :Ironuv   )  &&  freeze!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  freeze!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  freeze!(model, :Ironoptna)
    end
    GModelFit.update!(state.meval)

    println("\nFit known emission lines...")
    add_emission_lines!(recipe,state)
    add_patch_functs!(recipe, state)

    fit!(recipe, state)
    for cname in keys(state.user[:lcs])
        freeze!(model, cname)
    end

    println("\nFit nuisance emission lines...")
    add_nuisance_lines!(recipe, state)

    println("\nLast run with all parameters free...")
    thaw!(model, :QSOcont)
    haskey(model, :Galaxy   )  &&  thaw!(model, :Galaxy)
    haskey(model, :Balmer   )  &&  thaw!(model, :Balmer)
    haskey(model, :Ironuv   )  &&  thaw!(model, :Ironuv)
    haskey(model, :Ironoptbr)  &&  thaw!(model, :Ironoptbr)
    haskey(model, :Ironoptna)  &&  thaw!(model, :Ironoptna)
    for cname in keys(state.user[:lcs])
        thaw!(model, cname)
    end
    for j in 1:recipe.options[:n_nuisance]
        cname = Symbol(:nuisance, j)
        if model[cname].norm.val > 0
            thaw!(model, cname)
        else
            freeze!(model, cname)
        end
    end
    bestfit, fitstats = fit!(recipe, state)

    if neglect_weak_features!(recipe, state)
        println("\nRe-run fit...")
        bestfit, fitstats = fit!(recipe, state)
    end

    return bestfit, fitstats
end
