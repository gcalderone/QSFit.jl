export Type2

abstract type Type2 <: QSOGeneric end

function init_recipe!(recipe::Recipe{T}) where T <: Type2
    @invoke init_recipe!(recipe::Recipe{<: QSOGeneric})

    recipe.lines = [
        LineDescriptor(:Lya       , NarrowLine),
        LineDescriptor(:NV_1241   , NarrowLine),
        LineDescriptor(:CIV_1549  , NarrowLine),
        LineDescriptor(:CIII_1909 , NarrowLine),
        LineDescriptor(:MgII_2798 , NarrowLine),
        LineDescriptor(:NeV_3426  , NarrowLine),
        LineDescriptor(:OII_3727  , NarrowLine),
        LineDescriptor(:NeIII_3869, NarrowLine),
        LineDescriptor(:Hg        , NarrowLine),
        LineDescriptor(:Hb        , NarrowLine),
        LineDescriptor(:OIII_4959 , NarrowLine, BlueWing),
        LineDescriptor(:OIII_5007 , NarrowLine, BlueWing),
        LineDescriptor(:OI_6300   , NarrowLine),
        LineDescriptor(:OI_6364   , NarrowLine),
        LineDescriptor(:NII_6549  , NarrowLine),
        LineDescriptor(:Ha        , NarrowLine),
        LineDescriptor(:NII_6583  , NarrowLine),
        LineDescriptor(:SII_6716  , NarrowLine),
        LineDescriptor(:SII_6731  , NarrowLine)
    ]
end


function add_qso_continuum!(recipe::Recipe{<: Type2}, resid::GModelFit.Residuals)
    @invoke add_qso_continuum!(recipe::Recipe{<: QSOGeneric}, resid)
    resid.meval.model[:QSOcont].alpha.val = -1.8
    GModelFit.update!(resid.meval)
end


function set_constraints!(recipe::Recipe{<: Type2}, ::Type{NarrowLine}, comp::QSFit.AbstractSpecLineComp)
    @invoke set_constraints!(recipe::Recipe{<: QSOGeneric}, NarrowLine, comp)
    comp.fwhm.low = 10
end


function set_constraints!(recipe::Recipe{<: Type2}, ::Type{NuisanceLine}, comp::QSFit.AbstractSpecLineComp)
    @invoke set_constraints!(recipe::Recipe{<: QSOGeneric}, NuisanceLine, comp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 10, 500, 2e3
end


function add_patch_functs!(recipe::Recipe{<: Type2}, resid::GModelFit.Residuals)
    model = resid.meval.model
    # Patch parameters
    if haskey(model, :OIII_4959) && haskey(model, :OIII_5007)
        # model[:OIII_4959].norm.patch = @fd m -> m[:OIII_5007].norm / 3
        model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(model, :OIII_5007) && haskey(model, :OIII_5007_bw)
        model[:OIII_5007_bw].voff.patch = @fd (m, v) -> v + m[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm.patch = @fd (m, v) -> v + m[:OIII_5007].fwhm
        #model[:OIII_5007_bw].norm.patch = @fd (m, v) -> v + m[:OIII_5007].norm
    end
    if haskey(model, :OIII_4959_bw) && haskey(model, :OIII_5007_bw)
        model[:OIII_4959_bw].voff.patch = :OIII_5007_bw
        model[:OIII_4959_bw].fwhm.patch = :OIII_5007_bw
        # model[:OIII_4959_bw].norm.patch = @fd m -> m[:OIII_5007_bw].norm / 3
    end
    if haskey(model, :OI_6300) && haskey(model, :OI_6364)
        # model[:OI_6300].norm.patch = @fd m -> m[:OI_6364].norm / 3
        model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(model, :NII_6549) && haskey(model, :NII_6583)
        # model[:NII_6549].norm.patch = @fd m -> m[:NII_6583].norm / 3
        model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(model, :SII_6716) && haskey(model, :SII_6731)
        # model[:SII_6716].norm.patch = @fd m -> m[:SII_6731].norm / 3
        model[:SII_6716].voff.patch = :SII_6731
    end
end


function analyze(recipe::Recipe{<: Type2}, spec::Spectrum, resid::GModelFit.Residuals)
    recipe.spec = spec  # TODO: remove
    model = resid.meval.model
    model[:main] = SumReducer()
    select_maincomp!(model, :main)

    println("\nFit continuum components...")
    model[:Continuum] = SumReducer()
    push!(model[:main].list, :Continuum)
    add_qso_continuum!(recipe, resid)
    add_host_galaxy!(recipe, resid)
    fit!(recipe, resid)
    renorm_cont!(recipe, resid)
    freeze!(model, :QSOcont)
    haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
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
