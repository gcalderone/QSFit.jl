export Type2

abstract type Type2 <: QSOGeneric end

function init_recipe!(recipe::CRecipe{T}) where T <: Type2
    @track_recipe
    @invoke init_recipe!(recipe::CRecipe{<: QSOGeneric})

    recipe.n_nuisance = 2
    push!(recipe.nuisance_avoid,
          4959 .+ [-1,1] .* 25,
          5008 .+ [-1,1] .* 25)
end

function set_lines_dict!(recipe::CRecipe{T}) where T <: Type2
    @track_recipe
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


function add_qso_continuum!(recipe::CRecipe{<: Type2}, food::Food)
    @track_recipe
    @invoke add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, food.meval, data)
    food.model[:QSOcont].alpha.val = -1.8
    GModelFit.scan_model!(food.meval)
end


function set_constraints!(recipe::CRecipe{<: Type2}, ::Type{NarrowLine}, comp::QSFit.AbstractSpecLineComp)
    @track_recipe
    @invoke set_constraints!(recipe::CRecipe{<: QSOGeneric}, NarrowLine, comp)
    comp.fwhm.low = 10
end


function set_constraints!(recipe::CRecipe{<: Type2}, ::Type{NuisanceLine}, comp::QSFit.AbstractSpecLineComp)
    @track_recipe
    @invoke set_constraints!(recipe::CRecipe{<: QSOGeneric}, NuisanceLine, comp)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 10, 500, 2e3
end


function add_patch_functs!(recipe::CRecipe{<: Type2}, food::Food)
    @track_recipe
    # Patch parameters
    if haskey(food.model, :OIII_4959) && haskey(food.model, :OIII_5007)
        # food.model[:OIII_4959].norm.patch = @fd m -> m[:OIII_5007].norm / 3
        food.model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(food.model, :OIII_5007) && haskey(food.model, :OIII_5007_bw)
        food.model[:OIII_5007_bw].voff.patch = @fd (m, v) -> v + m[:OIII_5007].voff
        food.model[:OIII_5007_bw].fwhm.patch = @fd (m, v) -> v + m[:OIII_5007].fwhm
        #food.model[:OIII_5007_bw].norm.patch = @fd (m, v) -> v + m[:OIII_5007].norm
    end
    if haskey(food.model, :OIII_4959_bw) && haskey(food.model, :OIII_5007_bw)
        food.model[:OIII_4959_bw].voff.patch = :OIII_5007_bw
        food.model[:OIII_4959_bw].fwhm.patch = :OIII_5007_bw
        # food.model[:OIII_4959_bw].norm.patch = @fd m -> m[:OIII_5007_bw].norm / 3
    end
    if haskey(food.model, :OI_6300) && haskey(food.model, :OI_6364)
        # food.model[:OI_6300].norm.patch = @fd m -> m[:OI_6364].norm / 3
        food.model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(food.model, :NII_6549) && haskey(food.model, :NII_6583)
        # food.model[:NII_6549].norm.patch = @fd m -> m[:NII_6583].norm / 3
        food.model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(food.model, :SII_6716) && haskey(food.model, :SII_6731)
        # food.model[:SII_6716].norm.patch = @fd m -> m[:SII_6731].norm / 3
        food.model[:SII_6716].voff.patch = :SII_6731
    end
end


function analyze(recipe::CRecipe{<: Type2}, food::Food)
    @track_recipe
    food.model[:main] = SumReducer()
    food.meval = GModelFit.ModelEval(food.model, domain(data))

    println("\nFit continuum components...")
    food.model[:Continuum] = SumReducer()
    push!(food.model[:main].list, :Continuum)
    add_qso_continuum!(recipe, food)
    add_host_galaxy!(recipe, food)
    fit!(recipe, food)
    renorm_cont!(recipe, food)
    freeze!(food.model, :QSOcont)
    haskey(food.model, :Galaxy)  &&  freeze!(food.model, :Galaxy)
    GModelFit.scan_model!(food.meval)

    if (:lines in propertynames(recipe))  &&  (length(recipe.lines) > 0)
        println("\nFit known emission lines...")
        add_emission_lines!(recipe, food)
        add_patch_functs!(recipe, food)

        fit!(recipe, food)
        for cname in keys(recipe.lines)
            freeze!(food.model, cname)
        end
    end

    println("\nFit nuisance emission lines...")
    add_nuisance_lines!(recipe, food)

    println("\nLast run with all parameters free...")
    thaw!(food.model, :QSOcont)
    haskey(food.model, :Galaxy   )  &&  thaw!(food.model, :Galaxy)
    for cname in keys(recipe.lines)
        thaw!(food.model, cname)
    end
    for j in 1:recipe.n_nuisance
        cname = Symbol(:nuisance, j)
        if food.model[cname].norm.val > 0
            thaw!(food.model, cname)
        else
            freeze!(food.model, cname)
        end
    end
    bestfit, fsumm = fit!(recipe, food)

    if neglect_weak_features!(recipe, food)
        println("\nRe-run fit...")
        bestfit, fsumm = fit!(recipe, food)
    end

    return bestfit, fsumm
end
