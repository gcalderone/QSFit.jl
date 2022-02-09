function qsfit(source::QSO{TRecipe}) where TRecipe <: DefaultRecipe
    elapsed = time()
    @assert length(source.specs) == 1
    pspec = PreparedSpectrum(source, id=1)

    model = Model(pspec.domain)
    model[:Continuum] = SumReducer([])
    model[:main] = SumReducer([])
    push!(model[:main].list, :Continuum)
    select_reducer!(model, :main)
    delete!(model.revals, :default_sum)

    println(logio(source), "\nFit continuum components...")
    QSFit.add_qso_continuum!(source, pspec, model)
    QSFit.add_host_galaxy!(source, pspec, model)
    QSFit.add_balmer_cont!(source, pspec, model)
    fitres = fit!(source, model, pspec)
    QSFit.renorm_cont!(source, pspec, model)
    freeze(model, :qso_cont)
    haskey(model, :galaxy)  &&  freeze(model, :galaxy)
    haskey(model, :balmer)  &&  freeze(model, :balmer)
    evaluate!(model)

    println(logio(source), "\nFit iron templates...")
    model[:Iron] = SumReducer([])
    push!(model[:main].list, :Iron)
    QSFit.add_iron_uv!( source, pspec, model)
    QSFit.add_iron_opt!(source, pspec, model)

    if length(model[:Iron].list) > 0
        fitres = fit!(source, model, pspec)
        haskey(model, :ironuv   )  &&  freeze(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze(model, :ironoptna)
    end
    evaluate!(model)

    println(logio(source), "\nFit known emission lines...")
    QSFit.add_emission_lines!(source, pspec, model)
    QSFit.guess_emission_lines!(source, pspec, model)
    QSFit.add_patch_functs!(source, pspec, model)

    fitres = fit!(source, model, pspec)
    for lname in keys(pspec.lcs)
        freeze(model, lname)
    end

    println(logio(source), "\nFit unknown emission lines...")
    QSFit.add_unknown_lines!(source, pspec, model)

    println(logio(source), "\nLast run with all parameters free...")
    thaw(model, :qso_cont)
    haskey(model, :galaxy   )  &&  thaw(model, :galaxy)
    haskey(model, :balmer   )  &&  thaw(model, :balmer)
    haskey(model, :ironuv   )  &&  thaw(model, :ironuv)
    haskey(model, :ironoptbr)  &&  thaw(model, :ironoptbr)
    haskey(model, :ironoptna)  &&  thaw(model, :ironoptna)
    for lname in keys(pspec.lcs)
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
    fitres = fit!(source, model, pspec)

    if QSFit.neglect_weak_features!(source, pspec, model, fitres)
        println(logio(source), "\nRe-run fit...")
        fitres = fit!(source, model, pspec)
    end

    println(logio(source))
    show(logio(source), fitres)

    out = QSFitResults(source, pspec, model, fitres)
    elapsed = time() - elapsed
    println(logio(source), "\nElapsed time: $elapsed s")
    close_logio(source)
    return out
end
