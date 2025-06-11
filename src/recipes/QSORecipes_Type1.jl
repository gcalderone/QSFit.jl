export Type1

abstract type Type1 <: QSOGeneric end

# Special cases for emission lines
function line_component(recipe::CRecipe{T}, tid::Val{:MgII_2798}, ::Type{<: BroadLine}) where T <: Type1
    @track_recipe
    comp = @invoke line_component(recipe::CRecipe{<: supertype(T)}, tid::Val, BroadLine)
    comp.fwhm.val = 3000
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end


function init_recipe!(recipe::CRecipe{T}) where T <: Type1
    @track_recipe
    @invoke init_recipe!(recipe::CRecipe{<: QSOGeneric})
    recipe.min_spectral_coverage[:Ironuv]  = 0.3
    recipe.min_spectral_coverage[:Ironopt] = 0.3

    recipe.use_balmer = true
    recipe.use_ironuv = true;      recipe.Ironuv_fwhm    = 3000.
    recipe.use_ironopt = true;     recipe.Ironoptbr_fwhm = 3000.;  recipe.Ironoptna_fwhm =  500.
end


function lines_dict(recipe::CRecipe{T}) where T <: Type1
    @track_recipe

    # TODO: Identify the following feature
    (:l2420p0 in ATL.get_transition_ids())  ||  ATL.register(ATL.Permitted, "l2420p0", (2420., NaN, NaN))

    out = SpecLineSet()
    add_line!(recipe, out, :Lyb)
    # add_line!(recipe, out, :OV_1213)  # 1213.8A, Ferland+92, Shields+95
    add_line!(recipe, out, :Lya)
    # add_line!(recipe, out, :OV_1218)  # 1218.3A, Ferland+92, Shields+95
    add_line!(recipe, out, :NV_1241     , NarrowLine)
    add_line!(recipe, out, :OI_1306     , BroadLine)
    add_line!(recipe, out, :CII_1335    , BroadLine)
    add_line!(recipe, out, :SiIV_1400   , BroadLine)
    add_line!(recipe, out, :CIV_1549    , NarrowLine, BroadLine)
    add_line!(recipe, out, :HeII_1640   , BroadLine)
    add_line!(recipe, out, :OIII_1664   , BroadLine)
    add_line!(recipe, out, :AlIII_1858  , BroadLine)
    add_line!(recipe, out, :CIII_1909   , BroadLine)
    add_line!(recipe, out, :CII_2326    , BroadLine)
    add_line!(recipe, out, :l2420p0     , BroadLine)
    add_line!(recipe, out, :MgII_2798   , NarrowLine, BroadLine)
    add_line!(recipe, out, :NeV_3345)
    add_line!(recipe, out, :NeV_3426)
    add_line!(recipe, out, :OII_3727)
    add_line!(recipe, out, :NeIII_3869)
    add_line!(recipe, out, :Hd          , BroadLine)
    add_line!(recipe, out, :Hg          , BroadLine)
    add_line!(recipe, out, :OIII_4363)
    add_line!(recipe, out, :HeII_4686   , BroadLine)
    add_line!(recipe, out, :Hb)
    add_line!(recipe, out, :OIII_4959)
    add_line!(recipe, out, :OIII_5007   , ForbiddenLine, BlueWing)
    add_line!(recipe, out, :HeI_5876    , BroadLine)
    add_line!(recipe, out, :OI_6300)
    add_line!(recipe, out, :OI_6364)
    add_line!(recipe, out, :NII_6549)
    add_line!(recipe, out, :Ha          , NarrowLine, BroadLine, VeryBroadLine)
    add_line!(recipe, out, :NII_6583)
    add_line!(recipe, out, :SII_6716)
    add_line!(recipe, out, :SII_6731)
    return out
end


function add_balmer_cont!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_balmer_cont!.(Ref(recipe), Ref(fp), 1:length(fp.multi))
end
function add_balmer_cont!(recipe::CRecipe{<: Type1}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    if recipe.use_balmer
        getmodel(fp, ith)[:Balmer] = QSFit.balmercont(0.1, 0.5)
        push!(getmodel(fp, ith)[:Continuum].list, :Balmer)
        c = getmodel(fp, ith)[:Balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        getmodel(fp, ith)[:Balmer].norm.patch = @fd (m, v) -> v * m[:QSOcont].norm
    end
end


function add_iron_uv!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_iron_uv!.(Ref(recipe), Ref(fp), 1:length(fp.multi))
end
function add_iron_uv!(recipe::CRecipe{<: Type1}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    λ = coords(getdomain(fp, ith))
    if recipe.use_ironuv
        fwhm = recipe.Ironuv_fwhm
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.specs[ith].resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironuv, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            getmodel(fp, ith)[:Ironuv] = comp
            getmodel(fp, ith)[:Ironuv].norm.val = 1.
            push!(getmodel(fp, ith)[:Iron].list, :Ironuv)
            guess_norm_factor!(recipe, fp, ith, :Ironuv)
        else
            println("Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_iron_opt!.(Ref(recipe), Ref(fp), 1:length(fp.multi))
end
function add_iron_opt!(recipe::CRecipe{<: Type1}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    λ = coords(getdomain(fp, ith))
    if recipe.use_ironopt
        fwhm = recipe.Ironoptbr_fwhm
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.specs[ith].resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironopt, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            getmodel(fp, ith)[:Ironoptbr] = comp
            getmodel(fp, ith)[:Ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.Ironoptna_fwhm
            getmodel(fp, ith)[:Ironoptna] = QSFit.ironopt_narrow(fwhm)
            getmodel(fp, ith)[:Ironoptna].norm.val = 1 # TODO: guess a sensible value
            getmodel(fp, ith)[:Ironoptna].norm.fixed = false
            push!(getmodel(fp, ith)[:Iron].list, :Ironoptbr)
            push!(getmodel(fp, ith)[:Iron].list, :Ironoptna)
            guess_norm_factor!(recipe, fp, ith, :Ironoptbr)
        else
            println("Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_patch_functs!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_patch_functs!.(Ref(recipe), Ref(fp), 1:length(fp.multi))
end
function add_patch_functs!(recipe::CRecipe{<: Type1}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    model = getmodel(fp, ith)
    # Patch parameters
    if haskey(model, :OIII_4959)  &&  haskey(model, :OIII_5007)
        # model[:OIII_4959].norm.patch = @fd m -> m[:OIII_5007].norm / 3
        model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(model, :OIII_5007)  &&  haskey(model, :OIII_5007_bw)
        model[:OIII_5007_bw].voff.patch = @fd (m, v) -> v + m[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm.patch = @fd (m, v) -> v + m[:OIII_5007].fwhm
    end
    if haskey(model, :OI_6300)  &&  haskey(model, :OI_6364)
        # model[:OI_6300].norm.patch = @fd m -> m[:OI_6364].norm / 3
        model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(model, :NII_6549)  &&  haskey(model, :NII_6583)
        # model[:NII_6549].norm.patch = @fd m -> m[:NII_6583].norm / 3
        model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(model, :SII_6716)  &&  haskey(model, :SII_6731)
        # model[:SII_6716].norm.patch = @fd m -> m[:SII_6731].norm / 3
        model[:SII_6716].voff.patch = :SII_6731
    end

    if haskey(model, :Hb_na)  &&  haskey(model, :Ha_na)
        model[:Hb_na].voff.patch = :Ha_na
    end

    # The following are required to avoid degeneracy with iron
    # template
    if haskey(model, :Hg)  &&  haskey(model, :Hb_br)
        model[:Hg].voff.patch = :Hb_br
        model[:Hg].fwhm.patch = :Hb_br
    end
    if haskey(model, :Hg_br)  &&  haskey(model, :Hb_br)
        model[:Hg_br].voff.patch = :Hb_br
        model[:Hg_br].fwhm.patch = :Hb_br
    end
    if haskey(model, :Hg_na)  &&  haskey(model, :Hb_na)
        model[:Hg_na].voff.patch = :Hb_na
        model[:Hg_na].fwhm.patch = :Hb_na
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(model, :Hb_br)  &&
        haskey(model, :Hb_bb)
        model[:Hb_bb].norm.high = 1
        model[:Hb_bb].norm.val  = 0.5
        model[:Hb_bb].norm.patch = @fd (m, v) -> v * m[:Hb_br].norm / m[:Hb_br].fwhm * m[:Hb_bb].fwhm
    end
    if  haskey(model, :Ha_br)  &&
        haskey(model, :Ha_bb)
        model[:Ha_bb].norm.high = 1
        model[:Ha_bb].norm.val  = 0.5
        model[:Ha_bb].norm.patch = @fd (m, v) -> v * m[:Ha_br].norm / m[:Ha_br].fwhm * m[:Ha_bb].fwhm
    end
end


function analyze(recipe::CRecipe{T}, data::Measures{1}) where T <: Type1
    bestfit, fsumm = analyze(recipe, [data])
    return bestfit[1], fsumm
end

function analyze(recipe::CRecipe{T}, data::Vector{Measures{1}}) where T <: Type1
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
    add_balmer_cont!(recipe, fp)
    fit!(recipe, fp)
    for model in models
        freeze!(model, :QSOcont)
        haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
        haskey(model, :Balmer)  &&  freeze!(model, :Balmer)
    end
    renorm_cont!(recipe, fp)

    println("\nFit iron templates...")
    for model in models
        model[:Iron] = SumReducer()
        push!(model[:main].list, :Iron)
    end
    add_iron_uv!(recipe, fp)
    add_iron_opt!(recipe, fp)
    fit!(recipe, fp)
    for model in models
        haskey(model, :Ironuv   )  &&  freeze!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  freeze!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  freeze!(model, :Ironoptna)
    end

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
        haskey(model, :Balmer   )  &&  thaw!(model, :Balmer)
        haskey(model, :Ironuv   )  &&  thaw!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  thaw!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  thaw!(model, :Ironoptna)
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
