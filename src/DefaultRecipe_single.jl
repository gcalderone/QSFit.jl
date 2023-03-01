function run(job::JobState{T}) where T <: DefaultRecipe
    elapsed = time()
    model = job.model

    model[:Continuum] = SumReducer()
    model[:main] = SumReducer()
    push!(model[:main].list, :Continuum)
    select_maincomp!(model, :main)

    println(job.logio, "\nFit continuum components...")
    add_qso_continuum!(job)
    add_host_galaxy!(job)
    add_balmer_cont!(job)
    fit!(job)
    renorm_cont!(job)
    freeze!(model, :qso_cont)
    haskey(model, :galaxy)  &&  freeze!(model, :galaxy)
    haskey(model, :balmer)  &&  freeze!(model, :balmer)
    GFit.update!(model)

    println(job.logio, "\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!( job)
    add_iron_opt!(job)

    if length(model[:Iron].list) > 0
        fit!(job)
        haskey(model, :ironuv   )  &&  freeze!(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze!(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze!(model, :ironoptna)
    end
    GFit.update!(model)

    println(job.logio, "\nFit known emission lines...")
    add_emission_lines!(job)
    add_patch_functs!(job)

    fit!(job)
    for lname in keys(job.pspec.lcs)
        freeze!(model, lname)
    end

    println(job.logio, "\nFit unknown emission lines...")
    add_unknown_lines!(job)

    println(job.logio, "\nLast run with all parameters free...")
    thaw!(model, :qso_cont)
    haskey(model, :galaxy   )  &&  thaw!(model, :galaxy)
    haskey(model, :balmer   )  &&  thaw!(model, :balmer)
    haskey(model, :ironuv   )  &&  thaw!(model, :ironuv)
    haskey(model, :ironoptbr)  &&  thaw!(model, :ironoptbr)
    haskey(model, :ironoptna)  &&  thaw!(model, :ironoptna)
    for lname in keys(job.pspec.lcs)
        thaw!(model, lname)
    end
    for j in 1:job.options[:n_unk]
        cname = Symbol(:unk, j)
        if model[cname].norm.val > 0
            thaw!(model, cname)
        else
            freeze!(model, cname)
        end
    end
    bestfit, fitstats = fit!(job)

    if neglect_weak_features!(job)
        println(job.logio, "\nRe-run fit...")
        bestfit, fitstats = fit!(job)
    end

    println(job.logio)
    show(job.logio, fitstats)

    # Estimate line EWs
    EWs = estimate_line_EWs(job)

    out = JobResults(job, bestfit, fitstats, time() - elapsed)
    out.reduced[:EW] = EWs

    println(job.logio, "\nElapsed time: $(out.elapsed) s")
    close_log(job)
    return out
end
