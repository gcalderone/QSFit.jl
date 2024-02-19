export DefaultRecipe

abstract type DefaultRecipe <: AbstractRecipe end


get_cosmology(recipe::RRef{T}, state::State) where T <: DefaultRecipe =
    cosmology(h=0.70, OmegaM=0.3)   #S11


line_suffix(recipe::RRef{T}, ::Type{<: AbstractLine})  where T <: DefaultRecipe  = ""
line_suffix(recipe::RRef{T}, ::Type{<: NarrowLine})    where T <: DefaultRecipe  = "_na"
line_suffix(recipe::RRef{T}, ::Type{<: BroadLine})     where T <: DefaultRecipe  = "_br"
line_suffix(recipe::RRef{T}, ::Type{<: VeryBroadLine}) where T <: DefaultRecipe  = "_bb"

line_group( recipe::RRef{T}, ::Type{<: ForbiddenLine}) where T <: DefaultRecipe  = :NarrowLines
line_group( recipe::RRef{T}, ::Type{<: NarrowLine})    where T <: DefaultRecipe  = :NarrowLines
line_group( recipe::RRef{T}, ::Type{<: BroadLine})     where T <: DefaultRecipe  = :BroadLines
line_group( recipe::RRef{T}, ::Type{<: VeryBroadLine}) where T <: DefaultRecipe  = :VeryBroadLines

line_descriptors(recipe::RRef{T}) where T <: DefaultRecipe = recipe.options[:lines]

line_component(::RRef{T}, ::Val{:gauss}  , λ::Float64) where T <: DefaultRecipe = SpecLineGauss(λ)
line_component(::RRef{T}, ::Val{:lorentz}, λ::Float64) where T <: DefaultRecipe = SpecLineLorentz(λ)
line_component(::RRef{T}, ::Val{:voigt}  , λ::Float64) where T <: DefaultRecipe = SpecLineVoigt(λ)

line_component(recipe::RRef{T}, λ::Float64) where T <: DefaultRecipe =
    line_component(recipe, Val(recipe.options[:line_profiles]), λ)

function line_component(recipe::RRef{T}, ::Type{ForbiddenLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, λ)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 2e3
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end

function line_component(recipe::RRef{T}, ::Type{NarrowLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, λ)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 100, 5e2, 1e3 # avoid confusion with the broad component
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end

function line_component(recipe::RRef{T}, ::Type{BroadLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, λ)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 900, 5e3, 1.5e4
    comp.voff.low, comp.voff.val, comp.voff.high = -3e3, 0, 3e3
    return comp
end

function line_component(recipe::RRef{T}, ::Type{VeryBroadLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, λ)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 1e4, 2e4, 3e4
    comp.voff.fixed = true
    return comp
end

function line_component(recipe::RRef{T}, ::Type{NuisanceLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, λ)
    comp.norm.val = 0.
    comp.center.fixed = false;  comp.voff.fixed = true
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 600, 5e3, 1e4
    return comp
end

# Special cases
abstract type MgIIBroadLine <: BroadLine end
function line_component(recipe::RRef{T}, ::Type{MgIIBroadLine}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, supertype(MgIIBroadLine), λ)
    comp.voff.low, comp.voff.val, comp.voff.high = -1e3, 0, 1e3
    return comp
end


abstract type OIIIBlueWing  <: NarrowLine end
function line_component(recipe::RRef{T}, ::Type{OIIIBlueWing}, λ::Float64) where T <: DefaultRecipe
    comp = line_component(recipe, supertype(OIIIBlueWing), λ)
    comp.voff.low, comp.voff.val, comp.voff.high = 0, 0, 2e3
    return comp
end
line_suffix(recipe::RRef{T}, ::Type{OIIIBlueWing}) where T <: DefaultRecipe = :_bw


function set_default_options!(recipe::RRef{T}) where {T <: DefaultRecipe}
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :Ironuv  => 0.3,
                                       :Ironopt => 0.3)

    out[:host_template] = Dict(:library=>"swire", :template=>"Ell5")
    out[:host_template_ref_wavelength] = 5500. # A
    out[:use_host_template] = true
    out[:host_template_range] = [4000., 7000.]

    out[:use_balmer] = true
    out[:use_ironuv] = true;      out[:Ironuv_fwhm]    = 3000.
    out[:use_ironopt] = true;     out[:Ironoptbr_fwhm] = 3000.;  out[:Ironoptna_fwhm] =  500.

    out[:line_profiles] = :gauss
    out[:line_broadening] = true
    out[:Iron_broadening] = true

    out[:n_nuisance] = 10
    out[:nuisance_avoid] = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]  # Angstrom
    out[:nuisance_maxoffset_from_guess] = 1e3  # km/s

    out[:lines] = [
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
        LineDescriptor(:OIII_5007  , ForbiddenLine, OIIIBlueWing),
        LineDescriptor(:HeI_5876   , BroadLine),
        LineDescriptor(:OI_6300    , ForbiddenLine),
        LineDescriptor(:OI_6364    , ForbiddenLine),
        LineDescriptor(:NII_6549   , ForbiddenLine),
        LineDescriptor(:Ha         , NarrowLine, BroadLine, VeryBroadLine),
        LineDescriptor(:NII_6583   , ForbiddenLine),
        LineDescriptor(:SII_6716   , ForbiddenLine),
        LineDescriptor(:SII_6731   , ForbiddenLine)
    ]
    empty!(recipe.options)
    for (k, v) in out
        recipe.options[k] = v
    end
    nothing
end


function PreparedSpectrum(recipe::RRef{T}, state::State, source::Source; id=1) where T <: DefaultRecipe
    data = deepcopy(source.specs[id])
    println(state.logio, "Spectrum: " * data.label)
    goodfraction = count(data.good) / length(data.good)
    println(state.logio, "  good fraction:: ", goodfraction)
    if goodfraction < 0.5
        error("Good fraction < 0.5")
    end
    println(state.logio, "  resolution: ", @sprintf("%.4g", data.resolution), " km / s (FWHM)")

    λ = data.λ ./ (1 + source.z)
    data.good[findall(λ .< recipe.options[:wavelength_range][1])] .= false
    data.good[findall(λ .> recipe.options[:wavelength_range][2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is
    not sufficient the component should not be added to the model,
    and corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println(state.logio, "Good samples before line coverage filter: ", count(data.good))

    # Collect line components (neglecting the ones with insufficent coverage)
    lcs = line_components(recipe)
    good = deepcopy(data.good)
    for (cname, line) in lcs
        # Line coverage test
        threshold = get(recipe.options[:min_spectral_coverage], cname, recipe.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, line.comp)

        print(state.logio, @sprintf("Line %-15s coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax))
        if coverage < threshold
            print(state.logio, @sprintf(", threshold is < %5.3f, neglecting...", threshold))
            good[λmin .<= λ .< λmax] .= false
            delete!(lcs, cname)
        end
        println(state.logio)
    end
    data.good .= good

    # Second pass to neglect lines whose coverage has been affected by
    # the neglected spectral samples.
    for (cname, line) in lcs
        threshold = get(recipe.options[:min_spectral_coverage], cname, recipe.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, line.comp)
        if coverage < threshold
            print(state.logio, @sprintf("Line %-15s updated coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax))
            println(state.logio, @sprintf(", threshold is < %5.3f, neglecting...", threshold))
            delete!(lcs, cname)
            data.good[λmin .<= λ .< λmax] .= false
        end
    end
    println(state.logio, "Good samples after line coverage filter: ", count(data.good))

    # Sort lines according to center wavelength
    kk = collect(keys(lcs))
    vv = collect(values(lcs))
    ii = sortperm(getfield.(getfield.(getfield.(vv, :comp), :center), :val))
    lcs = OrderedDict(Pair.(kk[ii], vv[ii]))

    # De-reddening
    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(state.logio, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ld = uconvert(u"cm", luminosity_dist(get_cosmology(recipe, state), source.z))
    if source.z == 0
        flux2lum = 1  # Input spectrum is already given in proper units, no need to multiply
    else
        flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())
        println(state.logio, "Using flux-to-lum. conversion factor: ", flux2lum)
    end

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(dom,
                   data.flux[ii] .* dered[ii] .* flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* flux2lum .* (1 + source.z))
    return PreparedSpectrum(data.resolution, dom, lum, lcs, flux2lum)
end


function fit!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    mzer = GModelFit.cmpfit()
    mzer.Δfitstat_threshold = 1.e-5

    bestfit, fitstats = fit(state.model, state.pspec.data, minimizer=mzer)
    # show(state.logio, state.model)
    show(state.logio, fitstats)
    # @gp :QSFit state.pspec.data model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
    return bestfit, fitstats
end


function add_qso_continuum!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))

    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(state.pspec.data)) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    state.model[:QSOcont] = comp
    push!(state.model[:Continuum].list, :QSOcont)
    GModelFit.update!(state.model)
end


function add_host_galaxy!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))
    if recipe.options[:use_host_template]  &&
        (recipe.options[:host_template_range][1] .< maximum(λ))  &&
        (recipe.options[:host_template_range][2] .> minimum(λ))
        state.model[:Galaxy] = QSFit.hostgalaxy(recipe.options[:host_template][:template],
                                                library=recipe.options[:host_template][:library],
                                                refwl=recipe.options[:host_template_ref_wavelength])
        push!(state.model[:Continuum].list, :Galaxy)

        # Split total flux between continuum and host galaxy
        refwl = recipe.options[:host_template_ref_wavelength]
        vv = Dierckx.Spline1D(λ, values(state.pspec.data), k=1, bc="extrapolate")(refwl)
        @assert !isnan(vv) "Predicted L_λ at $(refwl)A is NaN"
        @assert vv > 0 "Predicted L_λ at $(refwl)A is negative"
        # (vv <= 0)  &&  (vv = .1 * median(values(state.pspec.data)))
        state.model[:Galaxy].norm.val    = 1/2 * vv
        state.model[:QSOcont].norm.val *= 1/2 * vv / Dierckx.Spline1D(λ, state.model(:QSOcont), k=1, bc="extrapolate")(refwl)
        GModelFit.update!(state.model)
    end
end


function add_balmer_cont!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    if recipe.options[:use_balmer]
        state.model[:Balmer] = QSFit.balmercont(0.1, 0.5)
        push!(state.model[:Continuum].list, :Balmer)
        c = state.model[:Balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        state.model[:Balmer].norm.patch = @λ (m, v) -> v * m[:QSOcont].norm
        GModelFit.update!(state.model)
    end
end


function renorm_cont!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    freeze!(state.model, :QSOcont)
    c = state.model[:QSOcont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println(state.logio, "Cont. norm. (before): ", c.norm.val)
        while true
            residuals = (state.model() - values(state.pspec.data)) ./ uncerts(state.pspec.data)
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= 0.99
            GModelFit.update!(state.model)
        end
        println(state.logio, "Cont. norm. (after) : ", c.norm.val)
    else
        println(state.logio, "Skipping cont. renormalization")
    end
    GModelFit.update!(state.model)
    # @gp (domain(state.model), state.pspec.data) state.model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
end


function guess_norm_factor!(recipe::RRef{T}, state::State, name::Symbol; quantile=0.95) where T <: DefaultRecipe
    @assert state.model[name].norm.val != 0
    m = state.model(name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    resid = values(state.pspec.data) - state.model()
    ratio = state.model[name].norm.val / sum(m[i1:i2])
    off = sum(resid[i1:i2]) * ratio
    state.model[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if state.model[name].norm.val < 0  # ensure line has positive normalization
        state.model[name].norm.val = abs(off)
    end
end


function add_iron_uv!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))
    resolution = state.pspec.resolution
    if recipe.options[:use_ironuv]
        fwhm = recipe.options[:Ironuv_fwhm]
        if recipe.options[:Iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(recipe.options[:min_spectral_coverage], :Ironuv, recipe.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            state.model[:Ironuv] = comp
            state.model[:Ironuv].norm.val = 1.
            push!(state.model[:Iron].list, :Ironuv)
            GModelFit.update!(state.model)
            QSFit.guess_norm_factor!(recipe, state, :Ironuv)
            GModelFit.update!(state.model)
        else
            println(state.logio, "Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))
    resolution = state.pspec.resolution
    if recipe.options[:use_ironopt]
        fwhm = recipe.options[:Ironoptbr_fwhm]
        if recipe.options[:Iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(recipe.options[:min_spectral_coverage], :Ironopt, recipe.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            state.model[:Ironoptbr] = comp
            state.model[:Ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.options[:Ironoptna_fwhm]
            if recipe.options[:Iron_broadening]
                fwhm = sqrt(fwhm^2 + resolution^2)
            end
            state.model[:Ironoptna] = QSFit.ironopt_narrow(fwhm)
            state.model[:Ironoptna].norm.val = 1 # TODO: guess a sensible value
            state.model[:Ironoptna].norm.fixed = false
            push!(state.model[:Iron].list, :Ironoptbr)
            push!(state.model[:Iron].list, :Ironoptna)
            GModelFit.update!(state.model)
            QSFit.guess_norm_factor!(recipe, state, :Ironoptbr)
            GModelFit.update!(state.model)
        else
            println(state.logio, "Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()

    # Create model components
    for (cname, line) in state.pspec.lcs
        # All QSFit line profiles take spectral resolution into
        # account.  This is significantly faster than convolving the
        # whole model with an instrument response but has some
        # limitations:
        # - works only with built-in profiles (Gaussian, Lorentzian, Voigt);
        # - all components must be additive (i.e. no absorptions);
        # - further narrow components (besides known emission lines)
        #   will not be corrected for instrumental resolution.
        if recipe.options[:line_broadening]
            line.comp.resolution = state.pspec.resolution
        end

        state.model[cname] = line.comp
        haskey(groups, line.group)  ||  (groups[line.group] = Vector{Symbol}())
        push!(groups[line.group], cname)
    end

    # Create reducers for groups
    for (group, lnames) in groups
        @assert group in [:BroadLines, :NarrowLines, :VeryBroadLines] "Unexpected group for emission lines: $group"
        state.model[group] = SumReducer(lnames)
        push!(state.model[:main].list, group)
    end
    GModelFit.update!(state.model)

    # Guess normalizations
    for group in [:BroadLines, :NarrowLines, :VeryBroadLines]  # Note: order is important
        if group in keys(groups)
            for cname in groups[group]
                QSFit.guess_norm_factor!(recipe, state, cname)
            end
        end
    end
end

function add_patch_functs!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    # Patch parameters
    if haskey(state.model, :OIII_4959)  &&  haskey(state.model, :OIII_5007)
        # state.model[:OIII_4959].norm.patch = @λ m -> m[:OIII_5007].norm / 3
        state.model[:OIII_4959].voff.patch = :OIII_5007
    end
    if haskey(state.model, :OIII_5007)  &&  haskey(state.model, :OIII_5007_bw)
        state.model[:OIII_5007_bw].voff.patch = @λ (m, v) -> v + m[:OIII_5007].voff
        state.model[:OIII_5007_bw].fwhm.patch = @λ (m, v) -> v + m[:OIII_5007].fwhm
    end
    if haskey(state.model, :OI_6300)  &&  haskey(state.model, :OI_6364)
        # state.model[:OI_6300].norm.patch = @λ m -> m[:OI_6364].norm / 3
        state.model[:OI_6300].voff.patch = :OI_6364
    end
    if haskey(state.model, :NII_6549)  &&  haskey(state.model, :NII_6583)
        # state.model[:NII_6549].norm.patch = @λ m -> m[:NII_6583].norm / 3
        state.model[:NII_6549].voff.patch = :NII_6583
    end
    if haskey(state.model, :SII_6716)  &&  haskey(state.model, :SII_6731)
        # state.model[:SII_6716].norm.patch = @λ m -> m[:SII_6731].norm / 3
        state.model[:SII_6716].voff.patch = :SII_6731
    end

    if haskey(state.model, :Hb_na)  &&  haskey(state.model, :Ha_na)
        state.model[:Hb_na].voff.patch = :Ha_na
    end

    # The following are required to avoid degeneracy with iron
    # template
    if haskey(state.model, :Hg)  &&  haskey(state.model, :Hb_br)
        state.model[:Hg].voff.patch = :Hb_br
        state.model[:Hg].fwhm.patch = :Hb_br
    end
    if haskey(state.model, :Hg_br)  &&  haskey(state.model, :Hb_br)
        state.model[:Hg_br].voff.patch = :Hb_br
        state.model[:Hg_br].fwhm.patch = :Hb_br
    end
    if haskey(state.model, :Hg_na)  &&  haskey(state.model, :Hb_na)
        state.model[:Hg_na].voff.patch = :Hb_na
        state.model[:Hg_na].fwhm.patch = :Hb_na
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(state.model, :Hb_br)  &&
        haskey(state.model, :Hb_bb)
        state.model[:Hb_bb].norm.high = 1
        state.model[:Hb_bb].norm.val  = 0.5
        state.model[:Hb_bb].norm.patch = @λ (m, v) -> v * m[:Hb_br].norm / m[:Hb_br].fwhm * m[:Hb_bb].fwhm
    end
    if  haskey(state.model, :Ha_br)  &&
        haskey(state.model, :Ha_bb)
        state.model[:Ha_bb].norm.high = 1
        state.model[:Ha_bb].norm.val  = 0.5
        state.model[:Ha_bb].norm.patch = @λ (m, v) -> v * m[:Ha_br].norm / m[:Ha_br].fwhm * m[:Ha_bb].fwhm
    end
end


function add_nuisance_lines!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    (recipe.options[:n_nuisance] > 0)  ||  (return nothing)

    # Prepare nuisance line components
    for i in 1:recipe.options[:n_nuisance]
        state.model[Symbol(:nuisance, i)] = line_component(recipe, NuisanceLine, 3000.)
    end
    state.model[:NuisanceLines] = SumReducer([Symbol(:nuisance, i) for i in 1:recipe.options[:n_nuisance]])
    push!(state.model[:main].list, :NuisanceLines)
    GModelFit.update!(state.model)
    for j in 1:recipe.options[:n_nuisance]
        freeze!(state.model, Symbol(:nuisance, j))
    end
    GModelFit.update!(state.model)

    # Set "nuisance" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λ = coords(domain(state.model))
    λnuisance = Vector{Float64}()
    while true
        (length(λnuisance) >= recipe.options[:n_nuisance])  &&  break
        GModelFit.update!(state.model)
        Δ = (values(state.pspec.data) - state.model()) ./ uncerts(state.pspec.data)

        # Avoid considering again the same region (within 1A) TODO: within resolution
        for l in λnuisance
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        for rr in recipe.options[:nuisance_avoid]
            Δ[findall(rr[1] .< λ .< rr[2])] .= 0.
        end

        # Do not add lines close to from the edges since these may
        # affect qso_cont fitting
        Δ[findall((λ .< minimum(λ)*1.02)  .|
                  (λ .> maximum(λ)*0.98))] .= 0.
        iadd = argmax(Δ)
        (Δ[iadd] <= 0)  &&  break  # No residual is greater than 0, skip further residuals....
        push!(λnuisance, λ[iadd])

        cname = Symbol(:nuisance, length(λnuisance))
        state.model[cname].norm.val = 1.
        state.model[cname].center.val  = λ[iadd]

        # Allow to shift by a quantity equal to ...
        @assert recipe.options[:nuisance_maxoffset_from_guess] > 0
        state.model[cname].center.low  = λ[iadd] * (1 - recipe.options[:nuisance_maxoffset_from_guess] / 3e5)
        state.model[cname].center.high = λ[iadd] * (1 + recipe.options[:nuisance_maxoffset_from_guess] / 3e5)

        # In any case, we must stay out of avoidance regions
        for rr in recipe.options[:nuisance_avoid]
            @assert !(rr[1] .< λ[iadd] .< rr[2])
            if rr[1] .>= λ[iadd]
                state.model[cname].center.high = min(state.model[cname].center.high, rr[1])
            end
            if rr[2] .<= λ[iadd]
                state.model[cname].center.low  = max(state.model[cname].center.low, rr[2])
            end
        end

        thaw!(state.model, cname)
        fit!(recipe, state)
        freeze!(state.model, cname)
    end
    GModelFit.update!(state.model)
end


function neglect_weak_features!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    # Disable "nuisance" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:recipe.options[:n_nuisance]
        cname = Symbol(:nuisance, ii)
        isfreezed(state.model, cname)  &&  continue
        if state.model[cname].norm.val == 0.
            freeze!(state.model, cname)
            needs_fitting = true
            println(state.logio, "Disabling $cname (norm. = 0)")
        elseif state.model[cname].norm.unc / state.model[cname].norm.val > 3
            state.model[cname].norm.val = 0.
            freeze!(state.model, cname)
            needs_fitting = true
            println(state.logio, "Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


function reduce(recipe::RRef{T}, state::State) where T <: AbstractRecipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(state.model())
    for cname in keys(state.pspec.lcs)
        haskey(state.model, cname) || continue
        cont .-= state.model(cname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for cname in keys(state.pspec.lcs)
        haskey(state.model, cname) || continue
        EW[cname] = int_tabulated(coords(domain(state.model)),
                                  state.model(cname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end


include("DefaultRecipe_single.jl")
# TODO include("DefaultRecipe_multi.jl")
