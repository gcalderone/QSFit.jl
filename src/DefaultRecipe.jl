export DefaultRecipe

abstract type DefaultRecipe <: AbstractRecipe end

# Use DefaultRecipe when no Job object is provided
run(source::Source) = run(source, QSFit.Job{DefaultRecipe}())


get_cosmology(         job::Job{T}     , args...; kws...) where T <: AbstractRecipe = get_cosmology(         T, job, args...; kws...)
EmLineComponent(       job::Job{T}     , args...; kws...) where T <: AbstractRecipe = EmLineComponent(       T, job, args...; kws...)
collect_linecomps(     job::Job{T}     , args...; kws...) where T <: AbstractRecipe = collect_linecomps(     T, job, args...; kws...)
StdSpectrum(           job::Job{T}     , args...; kws...) where T <: AbstractRecipe = StdSpectrum{T}(        T, job, args...; kws...)
minimizer(             job::JobState{T}, args...; kws...) where T <: AbstractRecipe = minimizer(             T, job, args...; kws...)
fit!(                  job::JobState{T}, args...; kws...) where T <: AbstractRecipe = fit!(                  T, job, args...; kws...)
add_qso_continuum!(    job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_qso_continuum!(    T, job, args...; kws...)
add_host_galaxy!(      job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_host_galaxy!(      T, job, args...; kws...)
add_balmer_cont!(      job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_balmer_cont!(      T, job, args...; kws...)
renorm_cont!(          job::JobState{T}, args...; kws...) where T <: AbstractRecipe = renorm_cont!(          T, job, args...; kws...)
guess_norm_factor!(    job::JobState{T}, args...; kws...) where T <: AbstractRecipe = guess_norm_factor!(    T, job, args...; kws...)
add_iron_uv!(          job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_iron_uv!(          T, job, args...; kws...)
add_iron_opt!(         job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_iron_opt!(         T, job, args...; kws...)
add_emission_lines!(   job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_emission_lines!(   T, job, args...; kws...)
guess_emission_lines!( job::JobState{T}, args...; kws...) where T <: AbstractRecipe = guess_emission_lines!( T, job, args...; kws...)
add_patch_functs!(     job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_patch_functs!(     T, job, args...; kws...)
add_unknown_lines!(    job::JobState{T}, args...; kws...) where T <: AbstractRecipe = add_unknown_lines!(    T, job, args...; kws...)
neglect_weak_features!(job::JobState{T}, args...; kws...) where T <: AbstractRecipe = neglect_weak_features!(T, job, args...; kws...)


function Options(::Type{DefaultRecipe})
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)

    out[:host_template] = Dict(:library=>"swire", :template=>"Ell5")
    out[:use_host_template] = true
    out[:host_template_range] = [4000., 7000.]

    out[:use_balmer] = true
    out[:use_ironuv] = true;      out[:ironuv_fwhm]    = 3000.
    out[:use_ironopt] = true;     out[:ironoptbr_fwhm] = 3000.;  out[:ironoptna_fwhm] =  500.

    out[:line_profiles] = :gauss
    out[:line_broadening] = true
    out[:iron_broadening] = true

    out[:n_unk] = 10
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]

    lines = OrderedDict{Symbol, EmLineDescription}()
    out[:lines] = lines

    lines[:Lyb         ] = StdEmLine(:Lyb       ,  narrow, broad)
    # lines[:OV_1213     ] = CustomEmLine(1213.8  ,  narrow)  # Ferland+92, Shields+95
    lines[:Lya         ] = StdEmLine(:Lya       ,  narrow, broad)
    # lines[:OV_1218     ] = CustomEmLine(1218.3  ,  narrow)  # Ferland+92, Shields+95
    lines[:NV_1241     ] = StdEmLine(:NV_1241   ,  narrow )
    lines[:OI_1306     ] = StdEmLine(:OI_1306   ,  broad  )
    lines[:CII_1335    ] = StdEmLine(:CII_1335  ,  broad  )
    lines[:SiIV_1400   ] = StdEmLine(:SiIV_1400 ,  broad  )
    lines[:CIV_1549    ] = StdEmLine(:CIV_1549  ,  narrow, broad)
    lines[:HeII_1640   ] = StdEmLine(:HeII_1640 ,  broad  )
    lines[:OIII_1664   ] = StdEmLine(:OIII_1664 ,  broad  )
    lines[:AlIII_1858  ] = StdEmLine(:AlIII_1858,  broad  )
    lines[:CIII_1909   ] = StdEmLine(:CIII_1909 ,  broad  )
    lines[:CII_2326    ] = StdEmLine(:CII_2326  ,  broad  )
    lines[:F2420       ] = CustomEmLine(2420.0  ,  broad  )
    lines[:MgII_2798   ] = StdEmLine(:MgII_2798 ,  narrow, broad)
    lines[:NeV_3345    ] = StdEmLine(:NeV_3345  ,  narrow )
    lines[:NeV_3426    ] = StdEmLine(:NeV_3426  ,  narrow )
    lines[:OII_3727    ] = StdEmLine(:OII_3727  ,  narrow )
    lines[:NeIII_3869  ] = StdEmLine(:NeIII_3869,  narrow )
    lines[:Hd          ] = StdEmLine(:Hd        ,  broad  )
    lines[:Hg          ] = StdEmLine(:Hg        ,  broad  )
    lines[:OIII_4363   ] = StdEmLine(:OIII_4363 ,  narrow )
    lines[:HeII_4686   ] = StdEmLine(:HeII_4686 ,  broad  )
    lines[:Hb          ] = StdEmLine(:Hb        ,  narrow, broad)
    lines[:OIII_4959   ] = StdEmLine(:OIII_4959 ,  narrow )
    lines[:OIII_5007   ] = StdEmLine(:OIII_5007 ,  narrow )
    lines[:OIII_5007_bw] = StdEmLine(:OIII_5007 ,  narrow )
    lines[:HeI_5876    ] = StdEmLine(:HeI_5876  ,  broad  )
    lines[:OI_6300     ] = StdEmLine(:OI_6300   ,  narrow )
    lines[:OI_6364     ] = StdEmLine(:OI_6364   ,  narrow )
    lines[:NII_6549    ] = StdEmLine(:NII_6549  ,  narrow )
    lines[:Ha          ] = StdEmLine(:Ha        ,  narrow, broad, verybroad)
    lines[:NII_6583    ] = StdEmLine(:NII_6583  ,  narrow )
    lines[:SII_6716    ] = StdEmLine(:SII_6716  ,  narrow )
    lines[:SII_6731    ] = StdEmLine(:SII_6731  ,  narrow )
    return out
end


get_cosmology(::Type{T}, job::Job) where T <: DefaultRecipe =
    cosmology(h=0.70, OmegaM=0.3)   #S11


function EmLineComponent(::Type{T}, job::Job, λ::Float64) where T <: DefaultRecipe
    @assert job.options[:line_profiles] in [:gauss, :lorentz, :voigt]
    if job.options[:line_profiles] == :gauss
        comp = SpecLineGauss(λ)
    elseif job.options[:line_profiles] == :lorentz
        comp = SpecLineLorentz(λ)
    elseif job.options[:line_profiles] == :voigt
        comp = SpecLineVoigt(λ)
    else
        error("options[:line_profiles] must be one of: :gauss, :lorentz, :voigt.")
    end
    return comp
end


function EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Narrow) where T <: DefaultRecipe
    comp = EmLineComponent(job, λ)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3
    return EmLineComponent{Narrow}(:_na, :NarrowLines, comp)
end


function EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Broad) where T <: DefaultRecipe
    comp = EmLineComponent(job, λ)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3
    return EmLineComponent{Broad}(:_br, :BroadLines, comp)
end


function EmLineComponent(::Type{T}, job::Job, λ::Float64, ::VeryBroad) where T <: DefaultRecipe
    comp = EmLineComponent(job, λ)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return EmLineComponent{VeryBroad}(:_bb, :VeryBroadLines, comp)
end


function EmLineComponent(::Type{T}, job::Job, λ::Float64, ::Unknown) where T <: DefaultRecipe
    comp = EmLineComponent(job, λ)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 600
    comp.fwhm.high = 1e4
    comp.voff.fixed = true
    return EmLineComponent{Unknown}(Symbol(""), :Unknown, comp)
end


function collect_linecomps(::Type{T}, job::Job) where T <: DefaultRecipe
    lcs = OrderedDict{Symbol, EmLineComponent}()
    for (name, linedesc) in job.options[:lines]
        if isa(linedesc, StdEmLine)
            tt = transitions(linedesc.tid)
            λ = getfield.(tt, :λ)
            if length(λ) > 1
                println(job.logio, "Considering average wavelength for the $(linedesc.tid) multiplet: " * string(mean(λ)) * "Å")
                λ = mean(λ)  # average lambda of multiplets
            else
                λ = λ[1]
            end
        else
            @assert isa(linedesc, CustomEmLine)
            λ = linedesc.λ
        end

        for ltype in linedesc.types
            lc = EmLineComponent(job, λ, ltype)

            if isa(linedesc, StdEmLine)
                if (ltype == Narrow)  &&  (length(linedesc.types) > 1)
                    lc.comp.fwhm.high = 1e3
                end
                if (ltype == Broad)  &&  (name == :MgII_2798)
                    lc.comp.voff.low  = -1e3
                    lc.comp.voff.high =  1e3
                end
                if (ltype == Narrow)  &&  (name == :OIII_5007_bw)
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


function StdSpectrum{T}(::Type{T}, job::Job, source::Source; id=1) where T <: DefaultRecipe
    data = deepcopy(source.specs[id])
    println(job.logio, "Spectrum: " * data.label)
    goodfraction = count(data.good) / length(data.good)
    println(job.logio, "  good fraction:: ", goodfraction)
    if goodfraction < 0.5
        error("Good fraction < 0.5")
    end
    println(job.logio, "  resolution: ", @sprintf("%.4g", data.resolution), " km / s (FWHM)")

    λ = data.λ ./ (1 + source.z)
    data.good[findall(λ .< job.options[:wavelength_range][1])] .= false
    data.good[findall(λ .> job.options[:wavelength_range][2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is
    not sufficient the component should not be added to the model,
    and corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println(job.logio, "Good samples before line coverage filter: ", count(data.good))

    # Collect line components (neglect components with insufficent coverage)
    lcs = collect_linecomps(job)
    good = deepcopy(data.good)
    for (cname, lc) in lcs
        # Line coverage test
        threshold = get(job.options[:min_spectral_coverage], cname, job.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, lc.comp)

        print(job.logio, @sprintf("Line %-15s coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax))
        if coverage < threshold
            print(job.logio, @sprintf(", threshold is < %5.3f, neglecting...", threshold))
            good[λmin .<= λ .< λmax] .= false
            delete!(lcs, cname)
        end
        println(job.logio)
    end
    data.good .= good

    # Second pass to neglect lines whose coverage has been affected by
    # the neglected spectral samples.
    for (cname, lc) in lcs
        threshold = get(job.options[:min_spectral_coverage], cname, job.options[:min_spectral_coverage][:default])
        (λmin, λmax, coverage) = spectral_coverage(λ[findall(data.good)], data.resolution, lc.comp)
        if coverage < threshold
            print(job.logio, @sprintf("Line %-15s updated coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax))
            println(job.logio, @sprintf(", threshold is < %5.3f, neglecting...", threshold))
            delete!(lcs, cname)
            data.good[λmin .<= λ .< λmax] .= false
        end
    end
    println(job.logio, "Good samples after line coverage filter: ", count(data.good))

    # Sort lines according to center wavelength
    kk = collect(keys(lcs))
    vv = collect(values(lcs))
    ii = sortperm(getfield.(getfield.(getfield.(vv, :comp), :center), :val))
    lcs = OrderedDict(Pair.(kk[ii], vv[ii]))

    # De-reddening
    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(job.logio, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ld = uconvert(u"cm", luminosity_dist(get_cosmology(job), source.z))
    flux2lum = 4pi * ld^2 * (scale_flux() * unit_flux()) / (scale_lum() * unit_lum())

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(dom,
                   data.flux[ii] .* dered[ii] .* flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* flux2lum .* (1 + source.z))

    return StdSpectrum{T}(data.resolution, dom, lum, lcs)
end


function minimizer(::Type{T}, job::JobState) where T <: DefaultRecipe
    mzer = GFit.cmpfit()
    mzer.Δfitstat_threshold = 1.e-5
    return mzer
end


function fit!(::Type{T}, job::JobState) where T <: DefaultRecipe
    mzer = minimizer(job)
    fitres = fit!(job.model, job.pspec.data, minimizer=mzer)
    # show(job.logio, job.model)
    show(job.logio, fitres)
    # @gp :QSFit job.pspec.data model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
    return fitres
end


# TODO function fit!(job::JobStateMulti{T}) where T <: DefaultRecipe
# TODO     mzer = minimizer(job)
# TODO     fitres = fit!(job.models, getfield.(job.pspecs, :data), minimizer=mzer)
# TODO     show(job.logio, job.model)
# TODO     show(job.logio, fitres)
# TODO     # for id in 1:length(job.model)
# TODO     #     @gp Symbol("QSFit$(id)") job.pspecs[id].data job.model[id]
# TODO     # end
# TODO     # printstyled(color=:blink, "Press ENTER to continue..."); readline()
# TODO     return fitres
# TODO end


function add_qso_continuum!(::Type{T}, job::JobState) where T <: DefaultRecipe
    λ = domain(job.model)[:]

    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = Spline1D(λ, job.pspec.data.val, k=1, bc="error")(comp.x0.val)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    job.model[:qso_cont] = comp
    push!(job.model[:Continuum].list, :qso_cont)
    evaluate!(job.model)
end


function add_host_galaxy!(::Type{T}, job::JobState) where T <: DefaultRecipe
    λ = domain(job.model)[:]
    if job.options[:use_host_template]  &&
        (job.options[:host_template_range][1] .< maximum(λ))  &&
        (job.options[:host_template_range][2] .> minimum(λ))
        job.model[:galaxy] = QSFit.hostgalaxy(job.options[:host_template])
        push!(job.model[:Continuum].list, :galaxy)

        # Split total flux between continuum and host galaxy
        vv = Spline1D(λ, job.pspec.data.val, k=1, bc="extrapolate")(5500.)
        job.model[:galaxy].norm.val    = 1/2 * vv
        job.model[:qso_cont].norm.val *= 1/2 * vv / Spline1D(λ, job.model(:qso_cont), k=1, bc="extrapolate")(5500.)
        evaluate!(job.model)
    end
end


function add_balmer_cont!(::Type{T}, job::JobState) where T <: DefaultRecipe
    if job.options[:use_balmer]
        job.model[:balmer] = QSFit.balmercont(0.1, 0.5)
        push!(job.model[:Continuum].list, :balmer)
        c = job.model[:balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        model = job.model # TODO: needed for @try_patch
        @try_patch! model[:balmer].norm *= model[:qso_cont].norm
        evaluate!(job.model)
    end
end


function renorm_cont!(::Type{T}, job::JobState) where T <: DefaultRecipe
    freeze(job.model, :qso_cont)
    c = job.model[:qso_cont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println(job.logio, "Cont. norm. (before): ", c.norm.val)
        while true
            residuals = (job.model() - job.pspec.data.val) ./ job.pspec.data.unc
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= 0.99
            evaluate!(job.model)
        end
        println(job.logio, "Cont. norm. (after) : ", c.norm.val)
    else
        println(job.logio, "Skipping cont. renormalization")
    end
    evaluate!(job.model)
    # @gp (domain(job.model), job.pspec.data) job.model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
end


function guess_norm_factor!(::Type{T}, job::JobState, name::Symbol; quantile=0.95) where T <: DefaultRecipe
    @assert job.model[name].norm.val != 0
    m = job.model(name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    resid = job.pspec.data.val - job.model()
    ratio = job.model[name].norm.val / sum(m[i1:i2])
    off = sum(resid[i1:i2]) * ratio
    job.model[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if job.model[name].norm.val < 0  # ensure line has positive normalization
        job.model[name].norm.val = abs(off)
    end
end


function add_iron_uv!(::Type{T}, job::JobState) where T <: DefaultRecipe
    λ = domain(job.model)[:]
    resolution = job.pspec.resolution
    if job.options[:use_ironuv]
        fwhm = job.options[:ironuv_fwhm]
        if job.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(job.options[:min_spectral_coverage], :ironuv, job.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            job.model[:ironuv] = comp
            job.model[:ironuv].norm.val = 1.
            push!(job.model[:Iron].list, :ironuv)
            evaluate!(job.model)
            QSFit.guess_norm_factor!(job, :ironuv)
            evaluate!(job.model)
        else
            println(job.logio, "Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(::Type{T}, job::JobState) where T <: DefaultRecipe
    λ = domain(job.model)[:]
    resolution = job.pspec.resolution
    if job.options[:use_ironopt]
        fwhm = job.options[:ironoptbr_fwhm]
        if job.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + resolution^2)
        end
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, resolution, comp)
        threshold = get(job.options[:min_spectral_coverage], :ironopt, job.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            job.model[:ironoptbr] = comp
            job.model[:ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = job.options[:ironoptna_fwhm]
            if job.options[:iron_broadening]
                fwhm = sqrt(fwhm^2 + resolution^2)
            end
            job.model[:ironoptna] = QSFit.ironopt_narrow(fwhm)
            job.model[:ironoptna].norm.val = 1 # TODO: guess a sensible value
            job.model[:ironoptna].norm.fixed = false
            push!(job.model[:Iron].list, :ironoptbr)
            push!(job.model[:Iron].list, :ironoptna)
            evaluate!(job.model)
            QSFit.guess_norm_factor!(job, :ironoptbr)
            evaluate!(job.model)
        else
            println(job.logio, "Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(::Type{T}, job::JobState) where T <: DefaultRecipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()
    for (cname, lc) in job.pspec.lcs
        # All QSFit line profiles take spectral resolution into
        # account.  This is significantly faster than convolving the
        # whole model with an instrument response but has some
        # limitations:
        # - works only with built-in profiles (Gaussian, Lorentzian, Voigt);
        # - all components must be additive (i.e. no absorptions);
        # - further narrow components (besides known emission lines)
        #   will not be corrected for instrumental resolution.
        if job.options[:line_broadening]
            lc.comp.resolution = job.pspec.resolution
        end

        job.model[cname] = lc.comp
        haskey(groups, lc.group)  ||  (groups[lc.group] = Vector{Symbol}())
        push!(groups[lc.group], cname)
    end
    for (group, lnames) in groups
        job.model[group] = SumReducer(lnames)
    end
    evaluate!(job.model)
end


function guess_emission_lines!(::Type{T}, job::JobState) where T <: DefaultRecipe
    groups_to_go = unique(getfield.(values(job.pspec.lcs), :group))
    for group in [:BroadLines, :NarrowLines, :VeryBroadLines]  # Note: order is important
        found = false
        for (cname, lc) in job.pspec.lcs
            (lc.group == group)  ||  continue
            QSFit.guess_norm_factor!(job, cname)
            found = true
        end
        if found
            push!(job.model[:main].list, group)
            deleteat!(groups_to_go, findfirst(groups_to_go .== group))
        end
        evaluate!(job.model)
    end

    # Ensure all groups have been considered
    @assert length(groups_to_go) == 0 "The following line groups have not been considered: $groups_to_go"
end


function add_patch_functs!(::Type{T}, job::JobState) where T <: DefaultRecipe
    model = job.model
    # Patch parameters
    @try_patch! begin
        # model[:OIII_4959].norm = model[:OIII_5007].norm / 3
        model[:OIII_4959].voff = model[:OIII_5007].voff
    end
    @try_patch! begin
        model[:OIII_5007_bw].voff += model[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm += model[:OIII_5007].fwhm
    end
    @try_patch! begin
        # model[:OI_6300].norm = model[:OI_6364].norm / 3
        model[:OI_6300].voff = model[:OI_6364].voff
    end
    @try_patch! begin
        # model[:NII_6549].norm = model[:NII_6583].norm / 3
        model[:NII_6549].voff = model[:NII_6583].voff
    end
    @try_patch! begin
        # model[:SII_6716].norm = model[:SII_6731].norm / 1.5
        model[:SII_6716].voff = model[:SII_6731].voff
    end

    @try_patch! model[:Hb_na].voff = model[:Ha_na].voff

    # The following are required to avoid degeneracy with iron
    # template
    @try_patch! begin
        model[:Hg].voff = model[:Hb_br].voff
        model[:Hg].fwhm = model[:Hb_br].fwhm
    end
    @try_patch! begin
        model[:Hg_br].voff = model[:Hb_br].voff
        model[:Hg_br].fwhm = model[:Hb_br].fwhm
    end
    @try_patch! begin
        model[:Hg_na].voff = model[:Hb_na].voff
        model[:Hg_na].fwhm = model[:Hb_na].fwhm
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(model, :Hb_br)  &&
        haskey(model, :Hb_bb)
        model[:Hb_bb].norm.high = 1
        model[:Hb_bb].norm.val  = 0.5
        @try_patch! model[:Hb_bb].norm *= model[:Hb_br].norm / model[:Hb_br].fwhm * model[:Hb_bb].fwhm
    end
    if  haskey(model, :Ha_br)  &&
        haskey(model, :Ha_bb)
        model[:Ha_bb].norm.high = 1
        model[:Ha_bb].norm.val  = 0.5
        @try_patch! model[:Ha_bb].norm *= model[:Ha_br].norm / model[:Ha_br].fwhm * model[:Ha_bb].fwhm
    end
end


function add_unknown_lines!(::Type{T}, job::JobState) where T <: DefaultRecipe
    (job.options[:n_unk] > 0)  ||  (return nothing)

    # Prepare unknown line components
    for i in 1:job.options[:n_unk]
        job.model[Symbol(:unk, i)] = EmLineComponent(job, 3000., unknown).comp
    end
    job.model[:UnkLines] = SumReducer([Symbol(:unk, i) for i in 1:job.options[:n_unk]])
    push!(job.model[:main].list, :UnkLines)
    evaluate!(job.model)
    for j in 1:job.options[:n_unk]
        freeze(job.model, Symbol(:unk, j))
    end
    evaluate!(job.model)

    # Set "unknown" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λ = domain(job.model)[:]
    λunk = Vector{Float64}()
    while true
        (length(λunk) >= job.options[:n_unk])  &&  break
        evaluate!(job.model)
        Δ = (job.pspec.data.val - job.model()) ./ job.pspec.data.unc

        # Avoid considering again the same region (within 1A) TODO: within resolution
        for l in λunk
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        for rr in job.options[:unk_avoid]
            Δ[findall(rr[1] .< λ .< rr[2])] .= 0.
        end

        # Do not add lines close to from the edges since these may
        # affect qso_cont fitting
        Δ[findall((λ .< minimum(λ)*1.02)  .|
                  (λ .> maximum(λ)*0.98))] .= 0.
        iadd = argmax(Δ)
        (Δ[iadd] <= 0)  &&  break  # No residual is greater than 0, skip further residuals....
        push!(λunk, λ[iadd])

        cname = Symbol(:unk, length(λunk))
        job.model[cname].norm.val = 1.
        job.model[cname].center.val  = λ[iadd]
        job.model[cname].center.low  = λ[iadd] - λ[iadd]/10. # allow to shift 10%
        job.model[cname].center.high = λ[iadd] + λ[iadd]/10.

        thaw(job.model, cname)
        fitres = fit!(job)
        freeze(job.model, cname)
    end
    evaluate!(job.model)
end


function neglect_weak_features!(::Type{T}, job::JobState) where T <: DefaultRecipe
    # Disable "unknown" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:job.options[:n_unk]
        cname = Symbol(:unk, ii)
        isfixed(job.model, cname)  &&  continue
        if job.model[cname].norm.val == 0.
            freeze(job.model, cname)
            needs_fitting = true
            println(job.logio, "Disabling $cname (norm. = 0)")
        elseif job.model[cname].norm.unc / job.model[cname].norm.val > 3
            job.model[cname].norm.val = 0.
            freeze(job.model, cname)
            needs_fitting = true
            println(job.logio, "Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


include("DefaultRecipe_single.jl")
# TODO include("DefaultRecipe_multi.jl")
