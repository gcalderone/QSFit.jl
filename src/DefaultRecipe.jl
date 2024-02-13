export DefaultRecipe

abstract type DefaultRecipe <: AbstractRecipe end


function set_default_options!(recipe::RRef{T}) where {T <: DefaultRecipe}
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)

    out[:host_template] = Dict(:library=>"swire", :template=>"Ell5")
    out[:host_template_ref_wavelength] = 5500. # A
    out[:use_host_template] = true
    out[:host_template_range] = [4000., 7000.]

    out[:use_balmer] = true
    out[:use_ironuv] = true;      out[:ironuv_fwhm]    = 3000.
    out[:use_ironopt] = true;     out[:ironoptbr_fwhm] = 3000.;  out[:ironoptna_fwhm] =  500.

    out[:line_profiles] = :gauss
    out[:line_broadening] = true
    out[:iron_broadening] = true

    out[:n_nuisance] = 10
    out[:nuisance_avoid] = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]  # Angstrom
    out[:nuisance_maxoffset_from_guess] = 1e3  # km/s

    lines = OrderedDict{Symbol, EmLineDescription}()
    out[:lines] = lines

    lines[:Lyb         ] = StdEmLine(:Lyb       ,  :narrow, :broad)
    # lines[:OV_1213     ] = CustomEmLine(1213.8  ,  :narrow)  # Ferland+92, Shields+95
    lines[:Lya         ] = StdEmLine(:Lya       ,  :narrow, :broad)
    # lines[:OV_1218     ] = CustomEmLine(1218.3  ,  :narrow)  # Ferland+92, Shields+95
    lines[:NV_1241     ] = StdEmLine(:NV_1241   ,  :narrow )
    lines[:OI_1306     ] = StdEmLine(:OI_1306   ,  :broad  )
    lines[:CII_1335    ] = StdEmLine(:CII_1335  ,  :broad  )
    lines[:SiIV_1400   ] = StdEmLine(:SiIV_1400 ,  :broad  )
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  ,  :narrow, :broad)
    lines[:HeII_1640   ] = StdEmLine(:HeII_1640 ,  :broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  :broad  )
    lines[:AlIII_1858  ] = StdEmLine(:AlIII_1858,  :broad  )
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 ,  :broad  )
    lines[:CII_2326    ] = StdEmLine(:CII_2326  ,  :broad  )
    lines[:F2420       ] = CustomEmLine(2420.0  ,  :broad  )
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 ,  :narrow, :broad)
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  :narrow )
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  :narrow )
    lines[:OII_3727    ] = StdEmLine(:OII_3727  ,  :narrow )
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869,  :narrow )
    lines[:Hd          ] = StdEmLine(:Hd        ,  :broad  )
    lines[:Hg          ] = StdEmLine(:Hg        ,  :broad  )
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  :narrow )
    lines[:HeII_4686   ] = StdEmLine(:HeII_4686 ,  :broad  )
    lines[:Hb          ] = StdEmLine(:Hb        ,  :narrow, :broad)
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 ,  :narrow )
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  :narrow )
    lines[:HeI_5876    ] = StdEmLine(:HeI_5876  ,  :broad  )
    lines[:OI_6300     ] = StdEmLine(:OI_6300   ,  :narrow )
    lines[:OI_6364     ] = StdEmLine(:OI_6364   ,  :narrow )
    lines[:NII_6549    ] = StdEmLine(:NII_6549  ,  :narrow )
    lines[:Ha          ] = StdEmLine(:Ha        ,  :narrow, :broad, :verybroad)
    lines[:NII_6583    ] = StdEmLine(:NII_6583  ,  :narrow )
    lines[:SII_6716    ] = StdEmLine(:SII_6716  ,  :narrow )
    lines[:SII_6731    ] = StdEmLine(:SII_6731  ,  :narrow )

    empty!(recipe.options)
    for (k, v) in out
        recipe.options[k] = v
    end
    nothing
end


get_cosmology(recipe::RRef{T}, state::State) where T <: DefaultRecipe =
    cosmology(h=0.70, OmegaM=0.3)   #S11


function EmLineComponent(recipe::RRef{T}, state::State, λ::Float64) where T <: DefaultRecipe
    @assert recipe.options[:line_profiles] in [:gauss, :lorentz, :voigt]
    if recipe.options[:line_profiles] == :gauss
        comp = SpecLineGauss(λ)
    elseif recipe.options[:line_profiles] == :lorentz
        comp = SpecLineLorentz(λ)
    elseif recipe.options[:line_profiles] == :voigt
        comp = SpecLineVoigt(λ)
    else
        error("options[:line_profiles] must be one of: :gauss, :lorentz, :voigt.")
    end
    return comp
end


function EmLineComponent(recipe::RRef{T}, state::State, λ::Float64, ::Val{:narrow}) where T <: DefaultRecipe
    comp = EmLineComponent(recipe, state, λ)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3
    return EmLineComponent{Val{:narrow}}(:_na, :NarrowLines, comp)
end


function EmLineComponent(recipe::RRef{T}, state::State, λ::Float64, ::Val{:broad}) where T <: DefaultRecipe
    comp = EmLineComponent(recipe, state, λ)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3
    return EmLineComponent{Val{:broad}}(:_br, :BroadLines, comp)
end


function EmLineComponent(recipe::RRef{T}, state::State, λ::Float64, ::Val{:verybroad}) where T <: DefaultRecipe
    comp = EmLineComponent(recipe, state, λ)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return EmLineComponent{Val{:verybroad}}(:_bb, :VeryBroadLines, comp)
end


function EmLineComponent(recipe::RRef{T}, state::State, λ::Float64, ::Val{:nuisance}) where T <: DefaultRecipe
    comp = EmLineComponent(recipe, state, λ)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 600
    comp.fwhm.high = 1e4
    comp.voff.fixed = true
    return EmLineComponent{Val{:nuisance}}(Symbol(""), :Nuisance, comp)
end


function collect_linecomps(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    lcs = OrderedDict{Symbol, EmLineComponent}()
    for (name, linedesc) in recipe.options[:lines]
        if isa(linedesc, StdEmLine)
            tt = transitions(linedesc.tid)
            λ = getfield.(tt, :λ)
            if length(λ) > 1
                println(state.logio, "Considering average wavelength for the $(linedesc.tid) multiplet: " * string(mean(λ)) * "Å")
                λ = mean(λ)  # average lambda of multiplets
            else
                λ = λ[1]
            end
        else
            @assert isa(linedesc, CustomEmLine)
            λ = linedesc.λ
        end

        for ltype in linedesc.types
            lc = EmLineComponent(recipe, state, λ, ltype)

            if isa(linedesc, StdEmLine)
                if (ltype == Val(:narrow))  &&  (length(linedesc.types) > 1)
                    lc.comp.fwhm.high = 1e3
                end
                if (ltype == Val(:broad))  &&  (name == :MgII_2798)
                    lc.comp.voff.low  = -1e3
                    lc.comp.voff.high =  1e3
                end
                if (ltype == Val(:narrow))  &&  (name == :OIII_5007_bw)
                    lc.comp.fwhm.val  = 500
                    lc.comp.fwhm.high = 1e3
                    lc.comp.voff.low  = 0
                    lc.comp.voff.high = 2e3
                end
            end

            if length(linedesc.types) == 1
                cname = name
            else
                cname = Symbol(name, lc.suffix)
            end
            @assert !(cname in keys(lcs))
            lcs[cname] = lc
        end
    end
    return lcs
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

    # Collect line components (neglectiing the ones with insufficent coverage)
    lcs = collect_linecomps(recipe, state)
    good = deepcopy(data.good)
    for (cname, lc) in lcs
        # Line coverage test
        threshold = get(recipe.options[:min_spectral_coverage], cname, recipe.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, lc.comp)

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
    for (cname, lc) in lcs
        threshold = get(recipe.options[:min_spectral_coverage], cname, recipe.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, lc.comp)
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


# TODO function fit!(state::StateMulti{T}) where T <: DefaultRecipe
# TODO     mzer = minimizer(state)
# TODO     fitstats = fit!(state.models, getfield.(state.pspecs, :data), minimizer=mzer)
# TODO     show(state.logio, state.model)
# TODO     show(state.logio, fitstats)
# TODO     # for id in 1:length(state.model)
# TODO     #     @gp Symbol("QSFit$(id)") state.pspecs[id].data state.model[id]
# TODO     # end
# TODO     # printstyled(color=:blink, "Press ENTER to continue..."); readline()
# TODO     return fitstats
# TODO end


function add_qso_continuum!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))

    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(state.pspec.data)) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    state.model[:qso_cont] = comp
    push!(state.model[:Continuum].list, :qso_cont)
    GModelFit.update!(state.model)
end


function add_host_galaxy!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    λ = coords(domain(state.model))
    if recipe.options[:use_host_template]  &&
        (recipe.options[:host_template_range][1] .< maximum(λ))  &&
        (recipe.options[:host_template_range][2] .> minimum(λ))
        state.model[:galaxy] = QSFit.hostgalaxy(recipe.options[:host_template][:template],
                                                library=recipe.options[:host_template][:library],
                                                refwl=recipe.options[:host_template_ref_wavelength])
        push!(state.model[:Continuum].list, :galaxy)

        # Split total flux between continuum and host galaxy
        refwl = recipe.options[:host_template_ref_wavelength]
        vv = Dierckx.Spline1D(λ, values(state.pspec.data), k=1, bc="extrapolate")(refwl)
        @assert !isnan(vv) "Predicted L_λ at $(refwl)A is NaN"
        @assert vv > 0 "Predicted L_λ at $(refwl)A is negative"
        # (vv <= 0)  &&  (vv = .1 * median(values(state.pspec.data)))
        state.model[:galaxy].norm.val    = 1/2 * vv
        state.model[:qso_cont].norm.val *= 1/2 * vv / Dierckx.Spline1D(λ, state.model(:qso_cont), k=1, bc="extrapolate")(refwl)
        GModelFit.update!(state.model)
    end
end


function add_balmer_cont!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    if recipe.options[:use_balmer]
        state.model[:balmer] = QSFit.balmercont(0.1, 0.5)
        push!(state.model[:Continuum].list, :balmer)
        c = state.model[:balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        state.model[:balmer].norm.patch = @λ (m, v) -> v * m[:qso_cont].norm
        GModelFit.update!(state.model)
    end
end


function renorm_cont!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    freeze!(state.model, :qso_cont)
    c = state.model[:qso_cont]
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
        fwhm = recipe.options[:ironuv_fwhm]
        if recipe.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(recipe.options[:min_spectral_coverage], :ironuv, recipe.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            state.model[:ironuv] = comp
            state.model[:ironuv].norm.val = 1.
            push!(state.model[:Iron].list, :ironuv)
            GModelFit.update!(state.model)
            QSFit.guess_norm_factor!(recipe, state, :ironuv)
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
        fwhm = recipe.options[:ironoptbr_fwhm]
        if recipe.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(recipe.options[:min_spectral_coverage], :ironopt, recipe.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            state.model[:ironoptbr] = comp
            state.model[:ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = recipe.options[:ironoptna_fwhm]
            if recipe.options[:iron_broadening]
                fwhm = sqrt(fwhm^2 + resolution^2)
            end
            state.model[:ironoptna] = QSFit.ironopt_narrow(fwhm)
            state.model[:ironoptna].norm.val = 1 # TODO: guess a sensible value
            state.model[:ironoptna].norm.fixed = false
            push!(state.model[:Iron].list, :ironoptbr)
            push!(state.model[:Iron].list, :ironoptna)
            GModelFit.update!(state.model)
            QSFit.guess_norm_factor!(recipe, state, :ironoptbr)
            GModelFit.update!(state.model)
        else
            println(state.logio, "Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(recipe::RRef{T}, state::State) where T <: DefaultRecipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()

    # Create model components
    for (cname, lc) in state.pspec.lcs
        # All QSFit line profiles take spectral resolution into
        # account.  This is significantly faster than convolving the
        # whole model with an instrument response but has some
        # limitations:
        # - works only with built-in profiles (Gaussian, Lorentzian, Voigt);
        # - all components must be additive (i.e. no absorptions);
        # - further narrow components (besides known emission lines)
        #   will not be corrected for instrumental resolution.
        if recipe.options[:line_broadening]
            lc.comp.resolution = state.pspec.resolution
        end

        state.model[cname] = lc.comp
        haskey(groups, lc.group)  ||  (groups[lc.group] = Vector{Symbol}())
        push!(groups[lc.group], cname)
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
        state.model[Symbol(:nuisance, i)] = EmLineComponent(recipe, state, 3000., Val(:nuisance)).comp
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
    for (lname, lc) in state.pspec.lcs
        haskey(state.model, lname) || continue
        cont .-= state.model(lname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for (lname, lc) in state.pspec.lcs
        haskey(state.model, lname) || continue
        EW[lname] = int_tabulated(coords(domain(state.model)),
                                  state.model(lname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end


include("DefaultRecipe_single.jl")
# TODO include("DefaultRecipe_multi.jl")
