export Type1

abstract type Type1 <: QSOGeneric end

# Special cases for emission lines
abstract type MgIIBroadLine <: BroadLine end
function line_component(recipe::CRecipe{<: Type1}, center::Float64, ::Type{<: MgIIBroadLine})
    @track_recipe
    comp = line_component(recipe, center, BroadLine)
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


function set_lines_dict!(recipe::CRecipe{T}) where T <: Type1
    @track_recipe
    (:lines in propertynames(recipe))  &&  (return get_lines_dict(recipe))
    add_line!(recipe, :Lyb)
    # add_line!(recipe, :OV_1213)  # 1213.8A, Ferland+92, Shields+95
    add_line!(recipe, :Lya)
    # add_line!(recipe, :OV_1218)  # 1218.3A, Ferland+92, Shields+95
    add_line!(recipe, :NV_1241    , NarrowLine)
    add_line!(recipe, :OI_1306    , BroadLine)
    add_line!(recipe, :CII_1335   , BroadLine)
    add_line!(recipe, :SiIV_1400  , BroadLine)
    add_line!(recipe, :CIV_1549   , NarrowLine, BroadLine)
    add_line!(recipe, :HeII_1640  , BroadLine)
    add_line!(recipe, :OIII_1664  , BroadLine)
    add_line!(recipe, :AlIII_1858 , BroadLine)
    add_line!(recipe, :CIII_1909  , BroadLine)
    add_line!(recipe, :CII_2326   , BroadLine)
    add_line!(recipe, QSFit.ATL.UnidentifiedTransition(2420.0), BroadLine)
    add_line!(recipe, :MgII_2798  , NarrowLine, MgIIBroadLine)
    add_line!(recipe, :NeV_3345)
    add_line!(recipe, :NeV_3426)
    add_line!(recipe, :OII_3727)
    add_line!(recipe, :NeIII_3869)
    add_line!(recipe, :Hd         , BroadLine)
    add_line!(recipe, :Hg         , BroadLine)
    add_line!(recipe, :OIII_4363)
    add_line!(recipe, :HeII_4686  , BroadLine)
    add_line!(recipe, :Hb)
    add_line!(recipe, :OIII_4959)
    add_line!(recipe, :OIII_5007  , ForbiddenLine, BlueWing)
    add_line!(recipe, :HeI_5876   , BroadLine)
    add_line!(recipe, :OI_6300)
    add_line!(recipe, :OI_6364)
    add_line!(recipe, :NII_6549)
    add_line!(recipe, :Ha         , NarrowLine, BroadLine, VeryBroadLine)
    add_line!(recipe, :NII_6583)
    add_line!(recipe, :SII_6716)
    add_line!(recipe, :SII_6731)
    return nothing
end


function add_balmer_cont!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_balmer_cont!.(Ref(recipe), dm)
    scan_and_evaluate!(dm)
end
function add_balmer_cont!(recipe::CRecipe{<: Type1}, meval::DataModel)
    @track_recipe
    if recipe.use_balmer
        getmodel(meval)[:Balmer] = QSFit.balmercont(0.1, 0.5)
        push!(getmodel(meval)[:Continuum].list, :Balmer)
        c = getmodel(meval)[:Balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        getmodel(meval)[:Balmer].norm.patch = @fd (m, v) -> v * m[:QSOcont].norm
    end
end


function add_iron_uv!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_iron_uv!.(Ref(recipe), dm)
end
function add_iron_uv!(recipe::CRecipe{<: Type1}, meval::DataModel)
    @track_recipe
    λ = coords(getdomain(meval))
    if recipe.use_ironuv
        fwhm = recipe.Ironuv_fwhm
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironuv, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            getmodel(meval)[:Ironuv] = comp
            getmodel(meval)[:Ironuv].norm.val = 1.
            push!(getmodel(meval)[:Iron].list, :Ironuv)
            scan_and_evaluate!(meval)
            guess_norm_factor!(recipe, meval, :Ironuv)
            scan_and_evaluate!(meval)
        else
            println("Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_iron_opt!.(Ref(recipe), dm)
    scan_and_evaluate!(dm)
end
function add_iron_opt!(recipe::CRecipe{<: Type1}, meval::DataModel)
    @track_recipe
    λ = coords(getdomain(meval))
    if recipe.use_ironopt
        fwhm = recipe.Ironoptbr_fwhm
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = QSFit.spectral_coverage(λ, recipe.spec.resolution, comp)
        threshold = get(recipe.min_spectral_coverage, :Ironopt, recipe.min_spectral_coverage[:default])
        if coverage >= threshold
            getmodel(meval)[:Ironoptbr] = comp
            getmodel(meval)[:Ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.Ironoptna_fwhm
            getmodel(meval)[:Ironoptna] = QSFit.ironopt_narrow(fwhm)
            getmodel(meval)[:Ironoptna].norm.val = 1 # TODO: guess a sensible value
            getmodel(meval)[:Ironoptna].norm.fixed = false
            push!(getmodel(meval)[:Iron].list, :Ironoptbr)
            push!(getmodel(meval)[:Iron].list, :Ironoptna)
            scan_and_evaluate!(meval)
            guess_norm_factor!(recipe, meval, :Ironoptbr)
        else
            println("Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_patch_functs!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_patch_functs!.(Ref(recipe), dm)
end
function add_patch_functs!(recipe::CRecipe{<: Type1}, meval::DataModel)
    @track_recipe
    model = getmodel(meval)
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


function analyze(recipe::CRecipe{T}, data::Vector{Measures{1}}) where T <: Type1
    @track_recipe
    model = Model(:main => SumReducer())
    select_maincomp!(model, :main)
    model[:Continuum] = SumReducer()
    push!(model[:main].list, :Continuum)
    meval = GModelFit.MEval(model, domain(data[1]))
    dm = [DataModel(meval, data, 1)]

    println("\nFit continuum components...")
    add_qso_continuum!(recipe, dm)
    add_host_galaxy!(recipe, dm)
    add_balmer_cont!(recipe, dm)
    fit!(recipe, dm)
    renorm_cont!(recipe, dm)
    freeze!(model, :QSOcont)
    haskey(model, :Galaxy)  &&  freeze!(model, :Galaxy)
    haskey(model, :Balmer)  &&  freeze!(model, :Balmer)
    scan_and_evaluate!(dm)

    println("\nFit iron templates...")
    model[:Iron] = SumReducer()
    push!(model[:main].list, :Iron)
    add_iron_uv!(recipe, dm)
    add_iron_opt!(recipe, dm)

    if length(model[:Iron].list) > 0
        fit!(recipe, dm)
        haskey(model, :Ironuv   )  &&  freeze!(model, :Ironuv)
        haskey(model, :Ironoptbr)  &&  freeze!(model, :Ironoptbr)
        haskey(model, :Ironoptna)  &&  freeze!(model, :Ironoptna)
    end
    scan_and_evaluate!(dm)

    println("\nFit known emission lines...")
    if (:lines in propertynames(recipe))  &&  (length(recipe.lines) > 0)
        add_emission_lines!(recipe, dm)
        add_patch_functs!(recipe, dm)

        fit!(recipe, dm)
        for cname in keys(recipe.lines)
            freeze!(model, cname)
        end
    end

    println("\nFit nuisance emission lines...")
    add_nuisance_lines!(recipe, dm)

    println("\nLast run with all parameters free...")
    thaw!(model, :QSOcont)
    haskey(model, :Galaxy   )  &&  thaw!(model, :Galaxy)
    haskey(model, :Balmer   )  &&  thaw!(model, :Balmer)
    haskey(model, :Ironuv   )  &&  thaw!(model, :Ironuv)
    haskey(model, :Ironoptbr)  &&  thaw!(model, :Ironoptbr)
    haskey(model, :Ironoptna)  &&  thaw!(model, :Ironoptna)
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
    scan_and_evaluate!(dm)
    bestfit, fsumm = fit!(recipe, dm)

    if any(neglect_weak_features!(recipe, dm))
        println("\nRe-run fit...")
        bestfit, fsumm = fit!(recipe, dm)
    end

    return bestfit, fsumm
end
