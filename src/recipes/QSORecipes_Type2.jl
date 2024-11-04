export Type2

abstract type Type2 <: QSOGeneric end

function init_recipe!(recipe::Recipe{T}) where T <: Type2
    @print_current_function
    @invoke init_recipe!(recipe::Recipe{<: QSOGeneric})

    recipe.n_nuisance = 2
    push!(recipe.nuisance_avoid,
          4959 .+ [-1,1] .* 25,
          5008 .+ [-1,1] .* 25)
end

function set_lines_dict!(recipe::Recipe{T}) where T <: Type2
    @print_current_function
    (:lines in propertynames(recipe))  &&  (return get_lines_dict(recipe))
    add_line!(recipe, :Lya       , NarrowLine)
    add_line!(recipe, :NV_1241   , ForbiddenLine)
    add_line!(recipe, :CIV_1549  , NarrowLine)
    add_line!(recipe, :CIII_1909 , NarrowLine)
    add_line!(recipe, :MgII_2798 , NarrowLine)
    add_line!(recipe, :NeV_3426  , ForbiddenLine)
    add_line!(recipe, :OII_3727  , ForbiddenLine)
    add_line!(recipe, :NeIII_3869, ForbiddenLine)
    add_line!(recipe, :Hg        , NarrowLine)
    add_line!(recipe, :Hb        , NarrowLine)
    add_line!(recipe, :OIII_4959 , ForbiddenLine, BlueWing)
    add_line!(recipe, :OIII_5007 , ForbiddenLine, BlueWing)
    add_line!(recipe, :OI_6300   , ForbiddenLine)
    add_line!(recipe, :OI_6364   , ForbiddenLine)
    add_line!(recipe, :NII_6549  , ForbiddenLine)
    add_line!(recipe, :Ha        , NarrowLine)
    add_line!(recipe, :NII_6583  , ForbiddenLine)
    add_line!(recipe, :SII_6716  , ForbiddenLine)
    add_line!(recipe, :SII_6731  , ForbiddenLine)
    nothing
end


function add_qso_continuum!(recipe::Recipe{<: Type2}, resid::GModelFit.Residuals)
    @print_current_function
    @invoke add_qso_continuum!(recipe::Recipe{<: QSOGeneric}, resid)
    resid.meval.model[:QSOcont].alpha.val = -1.8
    GModelFit.update!(resid.meval)
end


function set_constraints!(recipe::Recipe{<: Type2}, ::Type{NarrowLine}, comp::QSFit.AbstractSpecLineComp)
    @print_current_function
    @invoke set_constraints!(recipe::Recipe{<: QSOGeneric}, NarrowLine, comp)
    comp.fwhm.low = 10
end


function set_constraints!(recipe::Recipe{<: Type2}, ::Type{NuisanceLine}, comp::QSFit.AbstractSpecLineComp)
    @print_current_function
    @invoke set_constraints!(recipe::Recipe{<: QSOGeneric}, NuisanceLine, comp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 10, 500, 2e3
end


function add_patch_functs!(recipe::Recipe{<: Type2}, resid::GModelFit.Residuals)
    @print_current_function
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
    @print_current_function
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

    if (:lines in propertynames(recipe))  &&  (length(recipe.lines) > 0)
        println("\nFit known emission lines...")
        add_emission_lines!(recipe, resid)
        add_patch_functs!(recipe, resid)

        fit!(recipe, resid)
        for cname in keys(recipe.lines)
            freeze!(model, cname)
        end
    end

    println("\nFit nuisance emission lines...")
    add_nuisance_lines!(recipe, resid)

    println("\nLast run with all parameters free...")
    thaw!(model, :QSOcont)
    haskey(model, :Galaxy   )  &&  thaw!(model, :Galaxy)
    for cname in keys(recipe.lines)
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
