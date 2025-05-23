function qsfit_multi(source::QSO{TRecipe}; ref_id=1) where TRecipe <: DefaultRecipe
    elapsed = time()
    @assert length(source.specs) > 1
    @assert 1 <= ref_id <= length(source.specs)
    pspecs = [PreparedSpectrum(source, id=id) for id in 1:length(source.specs)]

    multi = Vector{Model}()
    println(logio(source), "\nFit continuum components...")
    for id in 1:length(pspecs)
        pspec = pspecs[id]
        model = Model(pspec.domain)
        model[:Continuum] = SumReducer()
        model[:main] = SumReducer()
        push!(model[:main].list, :Continuum)
        select_maincomp!(model, :main)

        QSFit.add_qso_continuum!(source, pspec, model)
        QSFit.add_host_galaxy!(source, pspec, model)
        QSFit.add_balmer_cont!(source, pspec, model)

        push!(multi, model)
        if id != ref_id
            @try_patch! multi[id][:galaxy].norm = multi[ref_id][:galaxy].norm
        end
    end
    fit!(source, multi, pspecs)

    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        QSFit.renorm_cont!(source, pspec, model)
        freeze!(model, :qso_cont)
        haskey(model, :galaxy)  &&  freeze!(model, :galaxy)
        haskey(model, :balmer)  &&  freeze!(model, :balmer)
        GModelFit.update!(model)
    end
    evaluate!(multi)

    println(logio(source), "\nFit iron templates...")
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        model[:Iron] = SumReducer()
        push!(model[:main].list, :Iron)
        QSFit.add_iron_uv!( source, pspec, model)
        QSFit.add_iron_opt!(source, pspec, model)

        if length(model[:Iron].list) > 0
            fit!(source, model, pspec)
            haskey(model, :ironuv   )  &&  freeze!(model, :ironuv)
            haskey(model, :ironoptbr)  &&  freeze!(model, :ironoptbr)
            haskey(model, :ironoptna)  &&  freeze!(model, :ironoptna)
        end
        GModelFit.update!(model)
    end
    evaluate!(multi)

    println(logio(source), "\nFit known emission lines...")
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        QSFit.add_emission_lines!(source, pspec, model)
        QSFit.add_patch_functs!(source, pspec, model)

        if id != ref_id
            @patch! multi[id][:OIII_5007].norm = multi[ref_id][:OIII_5007].norm
        end
    end
    fit!(source, multi, pspecs)
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        for lname in keys(pspec.lcs)
            freeze!(model, lname)
        end
    end
    evaluate!(multi)

    println(logio(source), "\nFit nuisance emission lines...")
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        QSFit.add_nuisance_lines!(source, pspec, model)
    end
    evaluate!(multi)

    println(logio(source), "\nLast run with all parameters free...")
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        thaw!(model, :qso_cont)
        haskey(model, :galaxy   )  &&  thaw!(model, :galaxy)
        haskey(model, :balmer   )  &&  thaw!(model, :balmer)
        haskey(model, :ironuv   )  &&  thaw!(model, :ironuv)
        haskey(model, :ironoptbr)  &&  thaw!(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  thaw!(model, :ironoptna)
        for lname in keys(pspec.lcs)
            thaw!(model, lname)
        end
        for j in 1:source.options[:n_nuisance]
            cname = Symbol(:nuisance, j)
            if model[cname].norm.val > 0
                thaw!(model, cname)
            else
                freeze!(model, cname)
            end
        end
    end
    bestfit, fsumm = fit!(source, multi, pspecs)

    rerun = false
    for id in 1:length(pspecs)
        model = multi[id]
        pspec = pspecs[id]
        rerun = rerun || QSFit.neglect_weak_features!(source, pspec, model, bestfit, fsumm)
    end
    if rerun
        println(logio(source), "\nRe-run fit...")
        bestfit, fsumm = fit!(source, multi, pspecs)
    end
    println(logio(source))
    show(logio(source), fsumm)

    out = QSFit.QSFitMultiResults(source, pspecs, multi, bestfit, fsumm)
    elapsed = time() - elapsed
    println(logio(source), "\nElapsed time: $elapsed s")
    close_logio(source)
    return out
end
