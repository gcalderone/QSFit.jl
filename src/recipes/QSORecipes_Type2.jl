export Type2

abstract type Type2 <: QSOGeneric end

# Special cases for emission lines
function line_component(recipe::CRecipe{T}, tid::Val{TID}, ::Type{<: NarrowLine}) where {TID, T <: Type2}
    @track_recipe
    comp = @invoke line_component(recipe::CRecipe{<: supertype(T)}, tid::Val, NarrowLine)
    comp.fwhm.low = 10
    return comp
end

function line_component(recipe::CRecipe{T}, tid::Val{TID}, ::Type{<: NuisanceLine}) where {TID, T <: Type2}
    @track_recipe
    comp = @invoke line_component(recipe::CRecipe{<: supertype(T)}, tid::Val, NarrowLine)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 10, 500, 2e3
    return comp
end

function init_recipe!(recipe::CRecipe{T}) where T <: Type2
    @track_recipe
    @invoke init_recipe!(recipe::CRecipe{<: QSOGeneric})

    recipe.n_nuisance = 2
    push!(recipe.nuisance_avoid,
          4959 .+ [-1,1] .* 25,
          5008 .+ [-1,1] .* 25)
end

function lines_dict(recipe::CRecipe{T}) where T <: Type2
    @track_recipe
    out = SpecLineSet()
    add_line!(recipe, out, :Lya       , NarrowLine)
    add_line!(recipe, out, :NV_1241   , ForbiddenLine)
    add_line!(recipe, out, :CIV_1549  , NarrowLine)
    add_line!(recipe, out, :CIII_1909 , NarrowLine)
    add_line!(recipe, out, :MgII_2798 , NarrowLine)
    add_line!(recipe, out, :NeV_3426  , ForbiddenLine)
    add_line!(recipe, out, :OII_3727  , ForbiddenLine)
    add_line!(recipe, out, :NeIII_3869, ForbiddenLine)
    add_line!(recipe, out, :Hg        , NarrowLine)
    add_line!(recipe, out, :Hb        , NarrowLine)
    add_line!(recipe, out, :OIII_4959 , ForbiddenLine, BlueWing)
    add_line!(recipe, out, :OIII_5007 , ForbiddenLine, BlueWing)
    add_line!(recipe, out, :OI_6300   , ForbiddenLine)
    add_line!(recipe, out, :OI_6364   , ForbiddenLine)
    add_line!(recipe, out, :NII_6549  , ForbiddenLine)
    add_line!(recipe, out, :Ha        , NarrowLine)
    add_line!(recipe, out, :NII_6583  , ForbiddenLine)
    add_line!(recipe, out, :SII_6716  , ForbiddenLine)
    add_line!(recipe, out, :SII_6731  , ForbiddenLine)
    return out
end


function add_qso_continuum!(recipe::CRecipe{<: Type2}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    @invoke add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, fp, ith)
    getmodel(fp, ith)[:QSOcont].alpha.val = -1.8
end


function add_patch_functs!(recipe::CRecipe{<: Type2}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    model = getmodel(fp, ith)
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


function analyze(recipe::CRecipe{T}, data::Measures{1}) where T <: Type2
    bestfit, fsumm = analyze(recipe, [data])
    return bestfit[1], fsumm
end

function analyze(recipe::CRecipe{T}, data::Vector{Measures{1}}) where T <: Type2
    @track_recipe
    models = Vector{Model}()
    for i in 1:length(data)
        model = Model(:main => SumReducer())
        select_maincomp!(model, :main)
        push!(models, model)
    end
    fp = GModelFit.FitProblem(models, data)

    println("\nFit continuum components...")
    for model in models
        model[:Continuum] = SumReducer()
        push!(model[:main].list, :Continuum)
    end
    add_qso_continuum!(recipe, fp)
    add_host_galaxy!(recipe, fp)
    fit!(recipe, fp)
    for model in models
        freeze!(model, :QSOcont)
        haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
    end
    renorm_cont!(recipe, fp)

    println("\nFit known emission lines...")
    for i in 1:length(models)
        model = models[i]
        lines = recipe.specs[i].meta[:lines]
        for group in unique(getfield.(values(lines), :group))
            model[group] = SumReducer()
            push!(model[:main].list, group)
        end
        for (cname, line) in lines
            push!(model[line.group].list, cname)
        end
    end
    add_emission_lines!(recipe, fp)
    add_patch_functs!(recipe, fp)
    fit!(recipe, fp)
    for i in 1:length(models)
        model = models[i]
        lines = recipe.specs[i].meta[:lines]
        for (cname, line) in lines
            freeze!(model, cname)
        end
    end

    println("\nFit nuisance emission lines...")
    fit_nuisance_lines!(recipe, fp)

    println("\nLast run with all parameters free...")
    for i in 1:length(models)
        model = models[i]
        lines = recipe.specs[i].meta[:lines]
        thaw!(model, :QSOcont)
        haskey(model, :Galaxy   )  &&  thaw!(model, :Galaxy)
        for cname in keys(lines)
            thaw!(model, cname)
        end
        if :NuisanceLines in keys(model)
            for cname in model[:NuisanceLines].list
                if model[cname].norm.val > 0
                    thaw!(model, cname)
                else
                    freeze!(model, cname)
                end
            end
        end
    end
    bestfit, fsumm = fit!(recipe, fp)

    if any(neglect_weak_features!(recipe, fp))
        println("\nRe-run fit...")
        bestfit, fsumm = fit!(recipe, fp)
    end

    return bestfit, fsumm
end
