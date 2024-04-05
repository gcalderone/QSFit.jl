export Type1

abstract type Type1 <: QSOGeneric end

# Special cases for emission lines
abstract type MgIIBroadLine <: BroadLine end
function set_constraints!(recipe::Recipe{<: Type1}, ::Type{MgIIBroadLine}, comp::GModelFit.AbstractComponent)
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
end


function init_recipe!(recipe::Recipe{T}) where T <: Type1
    @invoke init_recipe!(recipe::Recipe{<: QSOGeneric})
    recipe.min_spectral_coverage[:Ironuv]  = 0.3
    recipe.min_spectral_coverage[:Ironopt] = 0.3

    recipe.use_balmer = true
    recipe.use_ironuv = true;      recipe.Ironuv_fwhm    = 3000.
    recipe.use_ironopt = true;     recipe.Ironoptbr_fwhm = 3000.;  recipe.Ironoptna_fwhm =  500.

    recipe.lines = [
        LineDescriptor(:Lyb        , NarrowLine, BroadLine),
        # LineDescriptor(:OV_1213   , ForbiddenLine)  # 1213.8A, Ferland+92, Shields+95,
        LineDescriptor(:Lya        , NarrowLine, BroadLine),
        # LineDescriptor(:OV_1218   , ForbiddenLine)  # 1218.3A, Ferland+92, Shields+95,
        LineDescriptor(:NV_1241    , ForbiddenLine),
        LineDescriptor(:OI_1306    , BroadLine),
        LineDescriptor(:CII_1335   , BroadLine),
        LineDescriptor(:SiIV_1400  , BroadLine),
        LineDescriptor(:CIV_1549   , NarrowLine, BroadLine),
        LineDescriptor(:HeII_1640  , BroadLine),
        LineDescriptor(:OIII_1664  , BroadLine),
        LineDescriptor(:AlIII_1858 , BroadLine),
        LineDescriptor(:CIII_1909  , BroadLine),
        LineDescriptor(:CII_2326   , BroadLine),
        LineDescriptor(2420.0      , BroadLine),
        LineDescriptor(:MgII_2798  , NarrowLine, MgIIBroadLine),
        LineDescriptor(:NeV_3345   , ForbiddenLine),
        LineDescriptor(:NeV_3426   , ForbiddenLine),
        LineDescriptor(:OII_3727   , ForbiddenLine),
        LineDescriptor(:NeIII_3869 , ForbiddenLine),
        LineDescriptor(:Hd         , BroadLine),
        LineDescriptor(:Hg         , BroadLine),
        LineDescriptor(:OIII_4363  , ForbiddenLine),
        LineDescriptor(:HeII_4686  , BroadLine),
        LineDescriptor(:Hb         , NarrowLine, BroadLine),
        LineDescriptor(:OIII_4959  , ForbiddenLine),
        LineDescriptor(:OIII_5007  , ForbiddenLine, BlueWing),
        LineDescriptor(:HeI_5876   , BroadLine),
        LineDescriptor(:OI_6300    , ForbiddenLine),
        LineDescriptor(:OI_6364    , ForbiddenLine),
        LineDescriptor(:NII_6549   , ForbiddenLine),
        LineDescriptor(:Ha         , NarrowLine, BroadLine, VeryBroadLine),
        LineDescriptor(:NII_6583   , ForbiddenLine),
        LineDescriptor(:SII_6716   , ForbiddenLine),
        LineDescriptor(:SII_6731   , ForbiddenLine)
    ]    
end


function add_balmer_cont!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    if recipe.use_balmer
        resid.meval.model[:Balmer] = QSFit.balmercont(0.1, 0.5)
        push!(resid.meval.model[:Continuum].list, :Balmer)
        c = resid.meval.model[:Balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        resid.meval.model[:Balmer].norm.patch = @fd (m, v) -> v * m[:QSOcont].norm
        GModelFit.update!(resid.meval)
    end
end


function add_iron_uv!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)
    if recipe.use_ironuv
        fwhm = recipe.Ironuv_fwhm
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironuv, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            resid.meval.model[:Ironuv] = comp
            resid.meval.model[:Ironuv].norm.val = 1.
            push!(resid.meval.model[:Iron].list, :Ironuv)
            GModelFit.update!(resid.meval)
            guess_norm_factor!(recipe, resid, :Ironuv)
            GModelFit.update!(resid.meval)
        else
            println("Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(recipe::Recipe{<: Type1}, resid::GModelFit.Residuals)
    λ = coords(resid.meval.domain)
    if recipe.use_ironopt
        fwhm = recipe.Ironoptbr_fwhm
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironopt, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            resid.meval.model[:Ironoptbr] = comp
            resid.meval.model[:Ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.Ironoptna_fwhm
            resid.meval.model[:Ironoptna] = QSFit.ironopt_narrow(fwhm)
            resid.meval.model[:Ironoptna].norm.val = 1 # TODO: guess a sensible value
            resid.meval.model[:Ironoptna].norm.fixed = false
            push!(resid.meval.model[:Iron].list, :Ironoptbr)
            push!(resid.meval.model[:Iron].list, :Ironoptna)
            GModelFit.update!(resid.meval)
            guess_norm_factor!(recipe, resid, :Ironoptbr)
            GModelFit.update!(resid.meval)
        else
            println("Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function analyze(recipe::Recipe{<: Type1}, spec::Spectrum, resid::GModelFit.Residuals)
    recipe.spec = spec  # TODO: remove
    model = resid.meval.model
    model[:main] = SumReducer()
    select_maincomp!(model, :main)

    println("\nFit continuum components...")
    model[:Continuum] = SumReducer()
    push!(model[:main].list, :Continuum)
    add_qso_continuum!(recipe, resid)
    add_host_galaxy!(recipe, resid)
    add_balmer_cont!(recipe, resid)
    fit!(recipe, resid)
    renorm_cont!(recipe, resid)
    freeze!(model, :QSOcont)
    haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
    haskey(model, :Balmer)  &&  freeze!(model, :Balmer)
    GModelFit.update!(resid.meval)

    println("\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!(recipe, resid)
    add_iron_opt!(recipe, resid)

    if length(model[:Iron].list) > 0
        fit!(recipe, resid)
        haskey(model, :Ironuv   )  &&  freeze!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  freeze!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  freeze!(model, :Ironoptna)
    end
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
    haskey(model, :Balmer   )  &&  thaw!(model, :Balmer)
    haskey(model, :Ironuv   )  &&  thaw!(model, :Ironuv)
    haskey(model, :Ironoptbr)  &&  thaw!(model, :Ironoptbr)
    haskey(model, :Ironoptna)  &&  thaw!(model, :Ironoptna)
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
