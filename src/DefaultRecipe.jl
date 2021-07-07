abstract type DefaultRecipe <: AbstractRecipe end

function default_options(::Type{T}) where T <: DefaultRecipe
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)
    out[:skip_lines] = [:OIII_5007_bw] #, bb_Ha
    out[:host_template] = "Ell5"
    out[:use_host_template] = true
    out[:use_balmer] = true
    out[:use_ironuv] = true
    out[:use_ironopt] = true
    out[:n_unk] = 10
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50, 6565 .+ [-1,1] .* 150]
    out[:instr_broadening] = false
    out[:norm_integrated] = true
    return out
end

function known_spectral_lines(source::QSO{T}) where T <: DefaultRecipe
    list = OrderedDict{Symbol, AbstractSpectralLine}()
    list[:Lyb         ] = CombinedLine( 1026.0  , [BroadLine, NarrowLine])
    # list[:OV        ] = NarrowLine(   1213.8  )  # Ferland+92, Shields+95
    list[:Lya         ] = CombinedLine( 1215.24 , [BroadLine, NarrowLine])
    # list[:OV        ] = NarrowLine(   1218.3  )  # Ferland+92, Shields+95
    list[:NV_1241     ] = NarrowLine(   1240.81 )
    list[:OI_1306     ] = BroadLine(    1305.53 )
    list[:CII_1335    ] = BroadLine(    1335.31 )
    list[:SiIV_1400   ] = BroadLine(    1399.8  )
    list[:CIV_1549    ] = CombinedLine( 1549.48 , [BroadLine, NarrowLine])
    list[:HeII        ] = BroadLine(    1640.4  )
    list[:OIII        ] = BroadLine(    1665.85 )
    list[:AlIII       ] = BroadLine(    1857.4  )
    list[:CIII_1909   ] = BroadLine(    1908.734)
    list[:CII         ] = BroadLine(    2326.0  )
    list[:F2420       ] = BroadLine(    2420.0  )
    list[:MgII_2798   ] = CombinedLine( 2799.117, [BroadLine, NarrowLine])
    list[:NeVN        ] = NarrowLine(   3346.79 )
    list[:NeVI_3426   ] = NarrowLine(   3426.85 )
    list[:OII_3727    ] = NarrowLine(   3729.875)
    list[:NeIII_3869  ] = NarrowLine(   3869.81 )
    list[:Hd          ] = BroadLine(    4102.89 )
    list[:Hg          ] = BroadLine(    4341.68 )
    list[:OIII_4363   ] = NarrowLine(   4363.00 )  # TODO: Check wavelength is correct
    list[:HeII        ] = BroadLine(    4686.   )
    list[:Hb          ] = CombinedLine( 4862.68 , [BroadLine, NarrowLine])
    list[:OIII_4959   ] = NarrowLine(   4960.295)
    list[:OIII_5007   ] = NarrowLine(   5008.240)
    list[:OIII_5007_bw] = NarrowLine(   5008.240)
    list[:HeI_5876    ] = BroadLine(    5877.30 )
    list[:OI_6300     ] = NarrowLine(   6300.00 )  # TODO: Check wavelength is correct
    list[:OI_6364     ] = NarrowLine(   6364.00 )  # TODO: Check wavelength is correct
    list[:NII_6549    ] = NarrowLine(   6549.86 )
    list[:Ha          ] = CombinedLine( 6564.61 , [BroadLine, NarrowLine, BroadBaseLine])
    list[:NII_6583    ] = NarrowLine(   6585.27 )
    list[:SII_6716    ] = NarrowLine(   6718.29 )
    list[:SII_6731    ] = NarrowLine(   6732.67 )
    return list
end

function line_component(source::QSO{T}, name::Symbol, line::BroadLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3
    return LineComponent(nothing, line, comp, :BroadLines)
end

function line_component(source::QSO{T}, name::Symbol, line::NarrowLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3
    return LineComponent(nothing, line, comp, :NarrowLines)
end

function line_component(source::QSO{T}, name::Symbol, line::BroadBaseLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return LineComponent(nothing, line, comp, :BroadBaseLines)
end

function line_component(source::QSO{T}, name::Symbol, line::CombinedLine) where T <: DefaultRecipe
    out = OrderedDict{Symbol, LineComponent}()
    for tline in line.types
        lc = line_component(source, name, tline(line.λ))
        lc = LineComponent(line, lc.line, lc.comp, lc.reducer_name) # add parent type
        if tline == BroadLine
            out[Symbol(name, :_br)] = lc
        elseif tline == NarrowLine
            lc.comp.fwhm.high = 1e3
            out[Symbol(name, :_na)] = lc
        elseif tline == BroadBaseLine
            out[Symbol(name, :_bb)] = lc
        else
            error("Unsupported line type: $tline")
        end
    end
    return out
end

function line_component(source::QSO{T}, name::Symbol, line::UnkLine) where T <: DefaultRecipe
    comp = SpecLineGauss(line.λ)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 600
    comp.fwhm.high = 1e4
    comp.voff.fixed = true
    return LineComponent(nothing, typeof(line), comp, :UnknownLines)
end


function PreparedSpectrum(source::QSO{T}; id=1) where T <: DefaultRecipe
    data = source.specs[id]
    println(logio(source), "Spectrum: " * data.label)
    println(logio(source), "  good fraction:: ", goodfraction(data))
    if goodfraction(data) < 0.5
        error("Good fraction < 0.5")
    end
    println(logio(source), "  resolution: ", @sprintf("%.4g", data.resolution), " km / s")

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
    lcs = OrderedDict{Symbol, LineComponent}()
    for (lname, line) in known_spectral_lines(source)
        tmp = line_component(source, lname, line)
        isa(tmp, LineComponent)  &&  (tmp = OrderedDict(lname => tmp))
        @assert isa(tmp, OrderedDict{Symbol, LineComponent})
        for (lname2, lc) in tmp
            @assert !haskey(lcs, lname2)
            (λmin, λmax, coverage) = spectral_coverage(λ .* data.good, data.resolution, lc.comp)
            coverage = round(coverage * 1e3) / 1e3  # keep just 3 significant digits...
            threshold = get(source.options[:min_spectral_coverage], lname, source.options[:min_spectral_coverage][:default])
            print(logio(source), @sprintf("Line %-15s coverage: %4.2f (threshold: %4.2f)", lname2, coverage, threshold))
            if coverage < threshold
                print(logio(source), @sprintf("  neglecting range: %10.5g < λ < %10.5g", λmin, λmax))
                ii = findall(λmin .<= λ .< λmax)
                data.good[ii] .= false
            else
                lcs[lname2] = lc
            end
            println(logio(source))

        end
    end
    println(logio(source), "Good samples after line coverage filtering: ", length(findall(data.good)))

    # De-reddening
    dered = ccm_unred([1450, 3000, 5100.], source.mw_ebv)
    println(logio(source), "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, source.mw_ebv)

    ii = findall(data.good)
    dom = Domain(data.λ[ii] ./ (1 + source.z))
    lum = Measures(data.flux[ii] .* dered[ii] .* source.flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* source.flux2lum .* (1 + source.z))

    return PreparedSpectrum(id, data, dom, lum, lcs)
end


function minimizer(source::QSO{T}) where T <: DefaultRecipe
    mzer = GFit.cmpfit()
    mzer.Δfitstat_theshold = 1.e-5
    return mzer
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
        (minimum(λ) .< 5500 .< maximum(λ))
        model[:galaxy] = QSFit.hostgalaxy(source.options[:host_template])
        push!(model[:Continuum].list, :galaxy)

        # Split total flux between continuum and host galaxy
        vv = Spline1D(λ, pspec.data.val, k=1, bc="error")(5500.)
        model[:galaxy].norm.val  = 1/2 * vv
        model[:qso_cont].x0.val *= 1/2 * vv / Spline1D(λ, model(:qso_cont), k=1, bc="error")(5500.)
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
        @patch! model[:balmer].norm *= model[:qso_cont].norm
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
            (count(residuals .< 0) / length(residuals) > 0.9)  &&  break
            (c.norm.val < initialnorm / 5)  &&  break # give up
            c.norm.val *= 0.99
            evaluate!(model)
        end
        println(logio(source), "Cont. norm. (after) : ", c.norm.val)
    else
        println(logio(source), "Skipping cont. renormalization")
    end
    evaluate!(model)
end


function guess_norm_factor!(pspec::PreparedSpectrum, model::Model, name::Symbol; quantile=0.95)
    @assert model[name].norm.val != 0
    m = model(name)
    c = cumsum(m)
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    resid = pspec.data.val - model()
    ratio = model[name].norm.val / sum(m[i1:i2])
    model[name].norm.val += sum(resid[i1:i2]) * ratio
    if model[name].norm.val < 0
        @warn "$name component has negative normalization, set it to 0"
        model[name].norm.val = 0.
    end
end


function add_iron_uv!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]
    if source.options[:use_ironuv]
        fwhm = 3000.
        source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (pspec.orig.resolution * 2.355)^2))
        comp = QSFit.ironuv(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, pspec.orig.resolution, comp)
        threshold = get(source.options[:min_spectral_coverage], :ironuv, source.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            model[:ironuv] = comp
            model[:ironuv].norm.val = 1.
            push!(model[:Iron].list, :ironuv)
            evaluate!(model)
            guess_norm_factor!(pspec, model, :ironuv)
            evaluate!(model)
        else
            println(logio(source), "Ignoring ironuv component (threshold: $threshold)")
        end
    end
end


function add_iron_opt!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    λ = domain(model)[:]
    if source.options[:use_ironopt]
        fwhm = 3000.
        source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (pspec.orig.resolution * 2.355)^2))
        comp = QSFit.ironopt_broad(fwhm)
        (_1, _2, coverage) = spectral_coverage(λ, pspec.orig.resolution, comp)
        threshold = get(source.options[:min_spectral_coverage], :ironopt, source.options[:min_spectral_coverage][:default])
        if coverage >= threshold
            fwhm = 500.
            source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (pspec.orig.resolution * 2.355)^2))
            model[:ironoptbr] = comp
            model[:ironoptna] = QSFit.ironopt_narrow(fwhm)
            model[:ironoptbr].norm.val = 1 # TODO: guess a sensible value
            model[:ironoptna].norm.val = 0.0
            freeze(model, :ironoptna)  # will be freed during last run
            push!(model[:Iron].list, :ironoptbr)
            push!(model[:Iron].list, :ironoptna)
            evaluate!(model)
            guess_norm_factor!(pspec, model, :ironoptbr)
            evaluate!(model)
        else
            println(logio(source), "Ignoring ironopt component (threshold: $threshold)")
        end
    end
end


function add_emission_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()
    for (cname, lc) in pspec.lcs
        lc.comp.norm_integrated = source.options[:norm_integrated]

        # If instrumental broadening is not used and the line profile
        # is a Gaussian one take spectral resolution into account.
        # This is significantly faster than convolving with an
        # instrument response but has some limitations: - works only
        # with Gaussian profiles; - all components must be additive
        # (i.e. no absorptions) - further narrow components (besides
        # known emission lines) will not be corrected for instrumental
        # resolution
        if !source.options[:instr_broadening]
            if isa(lc.comp, QSFit.SpecLineGauss)
                lc.comp.spec_res_kms = pspec.orig.resolution
            else
                println(logio(source), "Line $cname is not a Gaussian profile: Can't take spectral resolution into account")
            end
        end

        model[cname] = lc.comp
        haskey(groups, lc.reducer_name)  ||  (groups[lc.reducer_name] = Vector{Symbol}())
        push!(groups[lc.reducer_name], cname)            
    end
    for (group, lnames) in groups
        model[group] = SumReducer(lnames)
    end
    
    if haskey(model, :MgII_2798)
        model[:MgII_2798].voff.low  = -1000
        model[:MgII_2798].voff.high =  1000
    end
    if haskey(model, :OIII_5007_bw)
        model[:OIII_5007_bw].fwhm.val  = 500
        model[:OIII_5007_bw].fwhm.low  = 1e2
        model[:OIII_5007_bw].fwhm.high = 1e3
        model[:OIII_5007_bw].voff.low  = 0
        model[:OIII_5007_bw].voff.high = 2e3
    end
    evaluate!(model)
end


function guess_emission_lines_values!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    for group in [:BroadLines, :NarrowLines, :BroadBaseLines]
        for (cname, lc) in pspec.lcs
            (lc.reducer_name == group)  ||  continue
            guess_norm_factor!(pspec, model, cname)
        end
        push!(model[:main].list, group)
        evaluate!(model)
    end
end


function add_patch_functs!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    # Patch parameters
    @patch! begin
         # model[:OIII_4959].norm = model[:OIII_5007].norm / 3
         model[:OIII_4959].voff = model[:OIII_5007].voff
     end
    @patch! begin
        model[:OIII_5007_bw].voff += model[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm += model[:OIII_5007].fwhm
    end
    @patch! begin
        # model[:OI_6300].norm = model[:OI_6364].norm / 3
        model[:OI_6300].voff = model[:OI_6364].voff
    end
    @patch! begin
        # model[:NII_6549].norm = model[:NII_6583].norm / 3
        model[:NII_6549].voff = model[:NII_6583].voff
    end
    @patch! begin
        # model[:SII_6716].norm = model[:SII_6731].norm / 1.5
        model[:SII_6716].voff = model[:SII_6731].voff
    end

    @patch! model[:na_Hb].voff = model[:na_Ha].voff

    # The following are required to avoid degeneracy with iron
    # template
    @patch! begin
        model[:Hg].voff = model[:br_Hb].voff
        model[:Hg].fwhm = model[:br_Hb].fwhm
    end
    @patch! begin
        model[:br_Hg].voff = model[:br_Hb].voff
        model[:br_Hg].fwhm = model[:br_Hb].fwhm
    end
    @patch! begin
        model[:na_Hg].voff = model[:na_Hb].voff
        model[:na_Hg].fwhm = model[:na_Hb].fwhm
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(model, :br_Hb)  &&
        haskey(model, :bb_Hb)
        model[:bb_Hb].norm.high = 1
        model[:bb_Hb].norm.val  = 0.5
        @patch! model[:bb_Hb].norm *= model[:br_Hb].norm / model[:br_Hb].fwhm * model[:bb_Hb].fwhm
    end
    if  haskey(model, :br_Ha)  &&
        haskey(model, :bb_Ha)
        model[:bb_Ha].norm.high = 1
        model[:bb_Ha].norm.val  = 0.5
        @patch! model[:bb_Ha].norm *= model[:br_Ha].norm / model[:br_Ha].fwhm * model[:bb_Ha].fwhm
    end
end


function add_unknown_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    (source.options[:n_unk] > 0)  &&  (return nothing)
    tmp = OrderedDict{Symbol, GFit.AbstractComponent}()
    for j in 1:source.options[:n_unk]
        tmp[Symbol(:unk, j)] = line_component(source, QSFit.UnkLine(5e3))
        tmp[Symbol(:unk, j)].norm_integrated = source.options[:norm_integrated]
    end
    for (cname, comp) in tmp
        model[cname] = comp
    end
    model[:UnkLines] = SumReducer(collect(keys(tmp)))
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
        bestfit = fit!(model, pspec.data, minimizer=mzer); show(logio(source), bestfit)
        freeze(model, cname)
    end
    evaluate!(model)
end


function neglect_weak_features!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    # Disable "unknown" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:source.options[:n_unk]
        cname = Symbol(:unk, ii)
        isfixed(model, cname)  &&  continue
        if bestfit[cname].norm.val == 0.
            freeze(model, cname)
            needs_fitting = true
            println(logio(source), "Disabling $cname (norm. = 0)")
        elseif bestfit[cname].norm.unc / bestfit[cname].norm.val > 3
            model[cname].norm.val = 0.
            freeze(model, cname)
            needs_fitting = true
            println(logio(source), "Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


function fit(source::QSO{TRecipe}) where TRecipe <: DefaultRecipe
    elapsed = time()
    @assert length(source.specs) == 1
    pspec = PreparedSpectrum(source, id=1)

    mzer = minimizer(source)
    model = Model(pspec.domain)
    model[:Continuum] = SumReducer([])
    model[:main] = SumReducer([])
    push!(model[:main].list, :Continuum)
    select_reducer!(model, :main)
    delete!(model.revals, :default_sum)

    # TODO if source.options[:instr_broadening]
    # TODO     GFit.set_instr_response!(model[1], (l, f) -> instrumental_broadening(l, f, source.spectra[id].resolution))
    # TODO end

    println(logio(source), "\nFit continuum components...")
    add_qso_continuum!(source, pspec, model)
    add_host_galaxy!(source, pspec, model)
    add_balmer_cont!(source, pspec, model)
    bestfit = fit!(model, pspec.data, minimizer=mzer);  show(logio(source), bestfit)
    renorm_cont!(source, pspec, model)
    freeze(model, :qso_cont)
    haskey(model, :galaxy)  &&  freeze(model, :galaxy)
    haskey(model, :balmer)  &&  freeze(model, :balmer)
    evaluate!(model)

    println(logio(source), "\nFit iron templates...")
    model[:Iron] = SumReducer([])
    push!(model[:main].list, :Iron)
    add_iron_uv!( source, pspec, model)
    add_iron_opt!(source, pspec, model)

    if length(model[:Iron].list) > 0
        bestfit = fit!(model, pspec.data, minimizer=mzer); show(logio(source), bestfit)
        haskey(model, :ironuv   )  &&  freeze(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze(model, :ironoptna)
    end
    evaluate!(model)

    println(logio(source), "\nFit known emission lines...")
    add_emission_lines!(source, pspec, model)
    guess_emission_lines_values!(source, pspec, model)
    add_patch_functs!(source, pspec, model)

    bestfit = fit!(model, pspec.data, minimizer=mzer); show(logio(source), bestfit)
    for lname in keys(pspec.lcs)
        freeze(model, lname)
    end

    println(logio(source), "\nFit unknown emission lines...")
    add_unknown_lines!(source, pspec, model)

    println(logio(source), "\nLast run with all parameters free...")
    thaw(model, :qso_cont)
    haskey(model, :galaxy   )  &&  thaw(model, :galaxy)
    haskey(model, :balmer   )  &&  thaw(model, :balmer)
    haskey(model, :ironuv   )  &&  thaw(model, :ironuv)
    haskey(model, :ironoptbr)  &&  thaw(model, :ironoptbr)
    haskey(model, :ironoptna)  &&  thaw(model, :ironoptna)
    for lname in keys(pspec.lcs)
        thaw(model, lname)
    end
    for j in 1:source.options[:n_unk]
        cname = Symbol(:unk, j)
        if model[cname].norm.val > 0
            thaw(model, cname)
        else
            freeze(model, cname)
        end
    end
    bestfit = fit!(model, pspec.data, minimizer=mzer)

    if neglect_weak_features!(source, pspec, model)
        println(logio(source), "\nRe-run fit...")
        bestfit = fit!(model, pspec.data, minimizer=mzer)
    end

    println(logio(source))
    show(logio(source), bestfit)

    out = QSFit.QSFitResults(source, pspec, model, bestfit)
    elapsed = time() - elapsed
    println(logio(source), "\nElapsed time: $elapsed s")
    QSFit.close_logio(source)
    return out
end
