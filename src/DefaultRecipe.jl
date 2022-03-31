export DefaultRecipe, qsfit, qsfit_multi

abstract type DefaultRecipe <: AbstractRecipe end

function default_options(::Type{T}) where T <: DefaultRecipe
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)
    out[:skip_lines] = Symbol[]

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

    return out
end

function known_spectral_lines(source::QSO{T}) where T <: DefaultRecipe
    list = [
        MultiCompLine(                      :Lyb                     , [BroadLine, NarrowLine]),
        # NarrowLine( custom_transition(tid=:OV_1213     , 1213.8  )),  # Ferland+92, Shields+95
        MultiCompLine(                      :Lya                     , [BroadLine, NarrowLine]),
        # NarrowLine( custom_transition(tid=:OV_1218     , 1218.3  )),  # Ferland+92, Shields+95
        NarrowLine(                         :NV_1241                ),
        BroadLine(                          :OI_1306                ),
        BroadLine(                          :CII_1335               ),
        BroadLine(                          :SiIV_1400              ),
        MultiCompLine(                      :CIV_1549                , [BroadLine, NarrowLine]),
        BroadLine(                          :HeII_1640              ),
        BroadLine(                          :OIII_1664              ),
        BroadLine(                          :AlIII_1858             ),
        BroadLine(                          :CIII_1909              ),
        BroadLine(                          :CII_2326               ),
        BroadLine(    custom_transition(tid=:F2420       , 2420.0  )),
        MultiCompLine(                      :MgII_2798               , [BroadLine, NarrowLine]),
        NarrowLine(                         :NeV_3345               ),
        NarrowLine(                         :NeV_3426               ),
        NarrowLine(                         :OII_3727               ),
        NarrowLine(                         :NeIII_3869             ),
        BroadLine(                          :Hd                     ),
        BroadLine(                          :Hg                     ),
        NarrowLine(                         :OIII_4363              ),
        BroadLine(                          :HeII_4686              ),
        MultiCompLine(                      :Hb                      , [BroadLine, NarrowLine]),
        NarrowLine(                         :OIII_4959              ),
        NarrowLine(                         :OIII_5007              ),
        NarrowLine(   cname=:OIII_5007_bw,  :OIII_5007              ),
        BroadLine(                          :HeI_5876               ),
        NarrowLine(                         :OI_6300                ),
        NarrowLine(                         :OI_6364                ),
        NarrowLine(                         :NII_6549               ),
        MultiCompLine(                      :Ha                      , [BroadLine, NarrowLine, BroadBaseLine]),
        NarrowLine(                         :NII_6583               ),
        NarrowLine(                         :SII_6716               ),
        NarrowLine(                         :SII_6731               )]
    return list
end


function LineComponent(source::QSO{T}, line::GenericLine, multicomp::Bool) where T <: DefaultRecipe
    tt = transition(line.tid)
    λ = tt.LAMBDA_VAC_ANG
    if length(λ) > 1
        println(logio(source), "Considering average wavelength for the $(line.tid) multiplet: " * string(mean(λ)) * "Å")
        λ = mean(λ)  # average lambda o multiplets
    else
        λ = λ[1]
    end
    @assert source.options[:line_profiles] in [:gauss, :lorentz, :voigt]
    if source.options[:line_profiles] == :gauss
        comp = SpecLineGauss(λ)
    elseif source.options[:line_profiles] == :lorentz
        comp = SpecLineLorentz(λ)
    elseif source.options[:line_profiles] == :voigt
        comp = SpecLineVoigt(λ)
    else
        error("options[:line_profiles] must be one of: :gauss, :lorentz, :voigt.")
    end
    return LineComponent(line, comp, multicomp)
end

function LineComponent(source::QSO{T}, line::BroadLine, multicomp::Bool) where T <: DefaultRecipe
    comp = LineComponent(source, GenericLine(line.tid), multicomp).comp
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3

    if line.tid == :MgII_2798
        comp.voff.low  = -1e3
        comp.voff.high =  1e3
    end
    return LineComponent(line, comp, multicomp)
end

function LineComponent(source::QSO{T}, line::NarrowLine, multicomp::Bool) where T <: DefaultRecipe
    comp = LineComponent(source, GenericLine(line.tid), multicomp).comp
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = (multicomp  ?  1e3  :  2e3)
    comp.voff.low  = -1e3
    comp.voff.high =  1e3

    if line.cname == :OIII_5007_bw
        comp.fwhm.val  = 500
        comp.fwhm.high = 1e3
        comp.voff.low  = 0
        comp.voff.high = 2e3
    end
    return LineComponent(line, comp, multicomp)
end

function LineComponent(source::QSO{T}, line::BroadBaseLine, multicomp::Bool) where T <: DefaultRecipe
    comp = LineComponent(source, GenericLine(line.tid), multicomp).comp
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return LineComponent(line, comp, multicomp)
end

function PreparedSpectrum(source::QSO{T}; id=1) where T <: DefaultRecipe
    data = deepcopy(source.specs[id])
    println(logio(source), "Spectrum: " * data.label)
    println(logio(source), "  good fraction:: ", goodfraction(data))
    if goodfraction(data) < 0.5
        error("Good fraction < 0.5")
    end
    println(logio(source), "  resolution: ", @sprintf("%.4g", data.resolution), " km / s (FWHM)")

    λ = data.λ ./ (1 + source.z)
    data.good[findall(λ .< source.options[:wavelength_range][1])] .= false
    data.good[findall(λ .> source.options[:wavelength_range][2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is
    not sufficient the component should not be added to the model,
    and corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println(logio(source), "Good samples before line coverage filter: ", length(findall(data.good)))

    # Collect LineComponent objects
    lcs = collect_LineComponent(source)

    # Sort lnames according to the expected line FWHM.  This is
    # necessary to avoid cases where a single line has sufficient
    # coverage, but the check for a subsequent broad lines has not
    # (dropping the "good" flag also for the narrow one).
    lnames = collect(keys(lcs))[reverse(sortperm(getfield.(getfield.(getfield.(values(lcs), :comp), :fwhm), :val)))]
    for lname in lnames
        lc = lcs[lname]
        (λmin, λmax, coverage) = spectral_coverage(λ .* data.good, data.resolution, lc.comp)
        coverage = round(coverage * 1e3) / 1e3  # keep just 3 significant digits...
        threshold = get(source.options[:min_spectral_coverage], lname, source.options[:min_spectral_coverage][:default])
        print(logio(source), @sprintf("Line %-15s coverage: %4.2f (threshold: %4.2f)", lname, coverage, threshold))
        if coverage < threshold
            print(logio(source), @sprintf("  neglecting range: %10.5g < λ < %10.5g", λmin, λmax))
            ii = findall(λmin .<= λ .< λmax)
            data.good[ii] .= false
            delete!(lcs, lname)
        end
        println(logio(source))
    end
    println(logio(source), "Good samples after line coverage filtering: ", length(findall(data.good)))

    # De-reddening
    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(logio(source), "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(dom,
                   data.flux[ii] .* dered[ii] .* source.flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* source.flux2lum .* (1 + source.z))

    return PreparedSpectrum(id, data, dom, lum, lcs)
end

function minimizer(source::QSO{T}) where T <: DefaultRecipe
    mzer = GFit.cmpfit()
    mzer.Δfitstat_threshold = 1.e-5
    return mzer
end

function fit!(source::QSO{T}, model::Model, pspec::PreparedSpectrum) where T <: DefaultRecipe
    mzer = minimizer(source)
    fitres = fit!(model, pspec.data, minimizer=mzer)
    show(logio(source), model)
    show(logio(source), fitres)
    # @gp :QSFit pspec.data model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
    return fitres
end


function fit!(source::QSO{T}, multi::MultiModel, pspecs::Vector{PreparedSpectrum}) where T <: DefaultRecipe
    mzer = minimizer(source)
    fitres = fit!(multi, getfield.(pspecs, :data), minimizer=mzer)
    show(logio(source), multi)
    show(logio(source), fitres)
    # for id in 1:length(multi)
    #     @gp Symbol("QSFit$(id)") pspecs[id].data multi[id]
    # end
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
    return fitres
end


function add_qso_continuum!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]

    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = Spline1D(λ, pspec.data.val, k=1, bc="error")(comp.x0.val)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    model[:qso_cont] = comp
    push!(model[:Continuum].list, :qso_cont)
    evaluate!(model)
end


function add_host_galaxy!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]
    if source.options[:use_host_template]  &&
        (source.options[:host_template_range][1] .< maximum(λ))  &&
        (source.options[:host_template_range][2] .> minimum(λ))
        model[:galaxy] = QSFit.hostgalaxy(source.options[:host_template])
        push!(model[:Continuum].list, :galaxy)

        # Split total flux between continuum and host galaxy
        vv = Spline1D(λ, pspec.data.val, k=1, bc="extrapolate")(5500.)
        model[:galaxy].norm.val    = 1/2 * vv
        model[:qso_cont].norm.val *= 1/2 * vv / Spline1D(λ, model(:qso_cont), k=1, bc="extrapolate")(5500.)
        evaluate!(model)
    end
end


function add_balmer_cont!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    if source.options[:use_balmer]
        model[:balmer] = QSFit.balmercont(0.1, 0.5)
        push!(model[:Continuum].list, :balmer)
        c = model[:balmer]
        c.norm.val  = 0.1
        c.norm.fixed = false
        c.norm.high = 0.5
        c.ratio.val = 0.5
        c.ratio.fixed = false
        c.ratio.low  = 0.1
        c.ratio.high = 1
        @try_patch! model[:balmer].norm *= model[:qso_cont].norm
        evaluate!(model)
    end
end


function renorm_cont!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    freeze(model, :qso_cont)
    c = model[:qso_cont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println(logio(source), "Cont. norm. (before): ", c.norm.val)
        while true
            residuals = (model() - pspec.data.val) ./ pspec.data.unc
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= 0.99
            evaluate!(model)
        end
        println(logio(source), "Cont. norm. (after) : ", c.norm.val)
    else
        println(logio(source), "Skipping cont. renormalization")
    end
    evaluate!(model)
    # @gp (domain(model), pspec.data) model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
end


function guess_norm_factor!(pspec::PreparedSpectrum, model::Model, name::Symbol; quantile=0.95)
    @assert model[name].norm.val != 0
    m = model(name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    resid = pspec.data.val - model()
    ratio = model[name].norm.val / sum(m[i1:i2])
    off = sum(resid[i1:i2]) * ratio
    model[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if model[name].norm.val < 0  # ensure line has positive normalization
        model[name].norm.val = abs(off)
    end
end


function add_iron_uv!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]
    if source.options[:use_ironuv]
        fwhm = source.options[:ironuv_fwhm]
        if source.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + pspec.orig.resolution^2)
        end
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, pspec.orig.resolution, comp)
        threshold = get(source.options[:min_spectral_coverage], :ironuv, source.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            model[:ironuv] = comp
            model[:ironuv].norm.val = 1.
            push!(model[:Iron].list, :ironuv)
            evaluate!(model)
            QSFit.guess_norm_factor!(pspec, model, :ironuv)
            evaluate!(model)
        else
            println(logio(source), "Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]
    if source.options[:use_ironopt]
        fwhm = source.options[:ironoptbr_fwhm]
        if source.options[:iron_broadening]
            fwhm = sqrt(fwhm^2 + pspec.orig.resolution^2)
        end
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, pspec.orig.resolution, comp)
        threshold = get(source.options[:min_spectral_coverage], :ironopt, source.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            model[:ironoptbr] = comp
            model[:ironoptbr].norm.val = 1 # TODO: guess a sensible value
            fwhm = source.options[:ironoptna_fwhm]
            if source.options[:iron_broadening]
                fwhm = sqrt(fwhm^2 + pspec.orig.resolution^2)
            end
            model[:ironoptna] = QSFit.ironopt_narrow(fwhm)
            model[:ironoptna].norm.val = 1 # TODO: guess a sensible value
            model[:ironoptna].norm.fixed = false
            push!(model[:Iron].list, :ironoptbr)
            push!(model[:Iron].list, :ironoptna)
            evaluate!(model)
            QSFit.guess_norm_factor!(pspec, model, :ironoptbr)
            evaluate!(model)
        else
            println(logio(source), "Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()
    for (cname, lc) in pspec.lcs

        # All QSFit line progiles take spectral resolution into
        # account.  This is significantly faster than convolving the
        # whole model with an instrument response but has some
        # limitations:
        # - works only with built-in profiles (Gaussian, Lorentzian, Voigt);
        # - all components must be additive (i.e. no absorptions);
        # - further narrow components (besides known emission lines)
        #   will not be corrected for instrumental resolution.
        if source.options[:line_broadening]
            lc.comp.resolution = pspec.orig.resolution
        end

        model[cname] = lc.comp
        grp = group(lc.line)
        haskey(groups, grp)  ||  (groups[grp] = Vector{Symbol}())
        push!(groups[grp], cname)
    end
    for (group, lnames) in groups
        model[group] = SumReducer(lnames)
    end
    evaluate!(model)
end


function guess_emission_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    groups_to_go = unique([QSFit.group(v.line) for (k,v) in pspec.lcs])
    for grp in [:BroadLines, :NarrowLines, :BroadBaseLines, :AsymmTailLines]  # Note: order is important
        found = false
        for (cname, lc) in pspec.lcs
            (grp == group(lc.line))  ||  continue
            QSFit.guess_norm_factor!(pspec, model, cname)
            found = true
        end
        if found
            push!(model[:main].list, grp)
            deleteat!(groups_to_go, findfirst(groups_to_go .== grp))
        end
        evaluate!(model)
    end

    # Ensure all groups have been considered
    @assert length(groups_to_go) == 0 "The following line groups have not been considered: $groups_to_go"
end


function add_patch_functs!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
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


function default_unk_line(source::QSO{T}) where T <: DefaultRecipe
    comp = SpecLineGauss(5e3)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 600
    comp.fwhm.high = 1e4
    comp.voff.fixed = true
    return comp
end


function add_unknown_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    (source.options[:n_unk] > 0)  ||  (return nothing)
    λ = domain(model)[:]
    for i in 1:source.options[:n_unk]
        model[Symbol(:unk, i)] = default_unk_line(source)
    end
    model[:UnkLines] = SumReducer([Symbol(:unk, i) for i in 1:source.options[:n_unk]])
    push!(model[:main].list, :UnkLines)
    evaluate!(model)
    for j in 1:source.options[:n_unk]
        freeze(model, Symbol(:unk, j))
    end
    evaluate!(model)

    # Set "unknown" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λunk = Vector{Float64}()
    while true
        (length(λunk) >= source.options[:n_unk])  &&  break
        evaluate!(model)
        Δ = (pspec.data.val - model()) ./ pspec.data.unc

        # Avoid considering again the same region (within 1A) TODO: within resolution
        for l in λunk
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        for rr in source.options[:unk_avoid]
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
        model[cname].norm.val = 1.
        model[cname].center.val  = λ[iadd]
        model[cname].center.low  = λ[iadd] - λ[iadd]/10. # allow to shift 10%
        model[cname].center.high = λ[iadd] + λ[iadd]/10.

        thaw(model, cname)
        fitres = fit!(source, model, pspec)
        freeze(model, cname)
    end
    evaluate!(model)
end


function neglect_weak_features!(source::QSO{T}, pspec::PreparedSpectrum, model::Model, fitres::GFit.FitResult) where T <: DefaultRecipe
    # Disable "unknown" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:source.options[:n_unk]
        cname = Symbol(:unk, ii)
        isfixed(model, cname)  &&  continue
        if model[cname].norm.val == 0.
            freeze(model, cname)
            needs_fitting = true
            println(logio(source), "Disabling $cname (norm. = 0)")
        elseif model[cname].norm.unc / model[cname].norm.val > 3
            model[cname].norm.val = 0.
            freeze(model, cname)
            needs_fitting = true
            println(logio(source), "Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


include("DefaultRecipe_single.jl")
include("DefaultRecipe_multi.jl")
