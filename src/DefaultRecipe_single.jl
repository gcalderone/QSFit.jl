
function analyze(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    model = state.model

    model[:Continuum] = SumReducer()
    model[:main] = SumReducer()
    push!(model[:main].list, :Continuum)
    select_maincomp!(model, :main)

    println(state.logio, "\nFit continuum components...")
    add_qso_continuum!(recipe, state)
    add_host_galaxy!(recipe, state)
    add_balmer_cont!(recipe, state)
    fit!(recipe, state)
    renorm_cont!(recipe, state)
    freeze!(model, :qso_cont)
    haskey(model, :galaxy)  &&  freeze!(model, :galaxy)
    haskey(model, :balmer)  &&  freeze!(model, :balmer)
    GModelFit.update!(model)

    println(state.logio, "\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!(recipe, state)
    add_iron_opt!(recipe, state)

    if length(model[:Iron].list) > 0
        fit!(recipe, state)
        haskey(model, :ironuv   )  &&  freeze!(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze!(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze!(model, :ironoptna)
    end
    GModelFit.update!(model)

    println(state.logio, "\nFit known emission lines...")
    add_emission_lines!(recipe,state)
    add_patch_functs!(recipe, state)

    fit!(recipe, state)
    for lname in keys(state.pspec.lcs)
        freeze!(model, lname)
    end

    println(state.logio, "\nFit unknown emission lines...")
    add_unknown_lines!(recipe, state)

    println(state.logio, "\nLast run with all parameters free...")
    thaw!(model, :qso_cont)
    haskey(model, :galaxy   )  &&  thaw!(model, :galaxy)
    haskey(model, :balmer   )  &&  thaw!(model, :balmer)
    haskey(model, :ironuv   )  &&  thaw!(model, :ironuv)
    haskey(model, :ironoptbr)  &&  thaw!(model, :ironoptbr)
    haskey(model, :ironoptna)  &&  thaw!(model, :ironoptna)
    for lname in keys(state.pspec.lcs)
        thaw!(model, lname)
    end
    for j in 1:recipe.options[:n_unk]
        cname = Symbol(:unk, j)
        if model[cname].norm.val > 0
            thaw!(model, cname)
        else
            freeze!(model, cname)
        end
    end
    bestfit, fitstats = fit!(recipe, state)

    if neglect_weak_features!(recipe, state)
        println(state.logio, "\nRe-run fit...")
        bestfit, fitstats = fit!(recipe, state)
    end

    println(state.logio)
    show(state.logio, fitstats)

    return bestfit, fitstats
end
