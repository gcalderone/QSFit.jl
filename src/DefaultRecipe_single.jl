function qsfit(source::Source, _job::Job{T}) where T <: DefaultRecipe
    job = JobState{T}(source, _job)
    elapsed = time()
    model = job.model

    model[:Continuum] = SumReducer([])
    model[:main] = SumReducer([])
    push!(model[:main].list, :Continuum)
    select_reducer!(model, :main)
    delete!(model.revals, :default_sum)

    println(job.logio, "\nFit continuum components...")
    add_qso_continuum!(job)
    add_host_galaxy!(job)
    add_balmer_cont!(job)
    fitres = fit!(job)
    renorm_cont!(job)
    freeze(model, :qso_cont)
    haskey(model, :galaxy)  &&  freeze(model, :galaxy)
    haskey(model, :balmer)  &&  freeze(model, :balmer)
    evaluate!(model)

    println(job.logio, "\nFit iron templates...")
    model[:Iron] = SumReducer([])
    push!(model[:main].list, :Iron)
    add_iron_uv!( job)
    add_iron_opt!(job)

    if length(model[:Iron].list) > 0
        fitres = fit!(job)
        haskey(model, :ironuv   )  &&  freeze(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze(model, :ironoptna)
    end
    evaluate!(model)

    println(job.logio, "\nFit known emission lines...")
    add_emission_lines!(job)
    guess_emission_lines!(job)
    add_patch_functs!(job)

    fitres = fit!(job)
    for lname in keys(job.pspec.lcs)
        freeze(model, lname)
    end

    println(job.logio, "\nFit unknown emission lines...")
    add_unknown_lines!(job)

    println(job.logio, "\nLast run with all parameters free...")
    thaw(model, :qso_cont)
    haskey(model, :galaxy   )  &&  thaw(model, :galaxy)
    haskey(model, :balmer   )  &&  thaw(model, :balmer)
    haskey(model, :ironuv   )  &&  thaw(model, :ironuv)
    haskey(model, :ironoptbr)  &&  thaw(model, :ironoptbr)
    haskey(model, :ironoptna)  &&  thaw(model, :ironoptna)
    for lname in keys(job.pspec.lcs)
        thaw(model, lname)
    end
    for j in 1:source.options[:n_unk]
        cname = Symbol(:unk, j)
        if model[cname].norm.val > 0
            thaw(model, cname)
        else
            freeze(model, cname)
        end
    end
    fitres = fit!(job)

    if neglect_weak_features!(job)
        println(job.logio, "\nRe-run fit...")
        fitres = fit!(job)
    end

    println(job.logio)
    show(job.logio, fitres)

    out = QSFitResults(job, fitres)
    elapsed = time() - elapsed
    println(job.logio, "\nElapsed time: $elapsed s")
    close_log(job)
    return out
end
