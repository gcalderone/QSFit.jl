calibsum(calib, args...) = calib .* (.+(args...))

function multi_fit(source::QSO{TRecipe}; ref_id=1) where TRecipe <: DefaultRecipe
    Nspec = length(source.domain)

    elapsed = time()
    mzer = GFit.cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Initialize components and guess initial values
    println(source.log, "\nFit continuum components...")
    preds = Vector{Prediction}()
    for id in 1:Nspec
        λ = source.domain[id][:]
        pred = Prediction(source.domain[id], :Continuum => SumReducer([:qso_cont]),
                          :qso_cont => QSFit.qso_cont_component(TRecipe))

        if source.options[:instr_broadening]
            GFit.set_instr_response!(pred, (l, f) -> instrumental_broadening(l, f, source.spectra[id].resolution))
        end

        push!(preds, pred)
        c = pred[:qso_cont]
        c.x0.val = median(λ)
        c.norm.val = Spline1D(λ, source.data[id].val, k=1, bc="error")(c.x0.val)
    end
    model = Model(preds)

    for id in 1:Nspec
        λ = source.domain[id][:]

        # Host galaxy template
        if source.options[:use_host_template]   &&
            (minimum(λ) .< 5500 .< maximum(λ))
            add!(model[id], :Continuum => SumReducer([:qso_cont, :galaxy]),
                 :galaxy => QSFit.hostgalaxy(source.options[:host_template]))
            model[id][:galaxy].norm.val = Spline1D(λ, source.data[id].val, k=1, bc="error")(5500.)
            if id != ref_id
                model[id][:galaxy].norm.fixed = true
                @patch! model m -> m[id][:galaxy].norm = m[ref_id][:galaxy].norm
            end
        end

        # Balmer continuum and pseudo-continuum
        if source.options[:use_balmer]
            tmp = [:qso_cont, :balmer]
            (:galaxy in keys(model[id]))  &&  push!(tmp, :galaxy)
            add!(model[id], :Continuum => SumReducer(tmp),
                 :balmer => QSFit.balmercont(0.1, 0.5))
            c = model[id][:balmer]
            c.norm.val  = 0.1
            c.norm.fixed = false
            c.norm.high = 0.5
            c.ratio.val = 0.5
            c.ratio.fixed = false
            c.ratio.low  = 0.1
            c.ratio.high = 1
            @patch! model[id] m -> m[:balmer].norm *= m[:qso_cont].norm
        end
    end

    bestfit = fit!(model, source.data, minimizer=mzer);  show(source.log, bestfit)

    # QSO continuum renormalization
    for id in 1:Nspec
        freeze(model[id], :qso_cont)
        c = model[id][:qso_cont]
        initialnorm = c.norm.val
        if c.norm.val > 0
            println(source.log, "$id: Cont. norm. (before): ", c.norm.val)
            while true
                residuals = (model[id]() - source.data[id].val) ./ source.data[id].unc
                (count(residuals .< 0) / length(residuals) > 0.9)  &&  break
                (c.norm.val < initialnorm / 5)  &&  break # give up
                c.norm.val *= 0.99
                evaluate!(model)
            end
            println(source.log, "$id : Cont. norm. (after) : ", c.norm.val)
        else
            println(source.log, "$id: Skipping cont. renormalization")
        end
        freeze(model[id], :qso_cont)
        (:galaxy in keys(model[id]))  &&  freeze(model[id], :galaxy)
        (:balmer in keys(model[id]))  &&  freeze(model[id], :balmer)
    end
    evaluate!(model)

    # Fit iron templates
    println(source.log, "\nFit iron templates...")
    for id in 1:Nspec
        λ = source.domain[id][:]

        iron_components = Vector{Symbol}()
        if source.options[:use_ironuv]
            fwhm = 3000.
            source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (source.spectra[id].resolution * 2.355)^2))
            comp = QSFit.ironuv(fwhm)
            (_1, _2, coverage) = spectral_coverage(λ, source.spectra[id].resolution, comp)
            threshold = get(source.options[:min_spectral_coverage], :ironuv, source.options[:min_spectral_coverage][:default])
            if coverage >= threshold
                add!(model[id], :ironuv => comp)
                model[id][:ironuv].norm.val = 0.5
                push!(iron_components, :ironuv)
            else
                println(source.log, "Ignoring ironuv component on prediction $id (threshold: $threshold)")
            end
        end

        if source.options[:use_ironopt]
            fwhm = 3000.
            source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (source.spectra[id].resolution * 2.355)^2))
            comp = QSFit.ironopt_broad(fwhm)
            (_1, _2, coverage) = spectral_coverage(λ, source.spectra[id].resolution, comp)
            threshold = get(source.options[:min_spectral_coverage], :ironopt, source.options[:min_spectral_coverage][:default])
            if coverage >= threshold
                fwhm = 500.
                source.options[:instr_broadening]  ||  (fwhm = sqrt(fwhm^2 + (source.spectra[id].resolution * 2.355)^2))
                add!(model[id],
                     :ironoptbr => comp,
                     :ironoptna => QSFit.ironopt_narrow(fwhm))
                model[id][:ironoptbr].norm.val = 0.5
                model[id][:ironoptna].norm.val = 0.0
                freeze(model[id], :ironoptna)  # will be freed during last run
                push!(iron_components, :ironoptbr, :ironoptna)
            else
                println(source.log, "Ignoring ironopt component on prediction $id (threshold: $threshold)")
            end
        end
        if length(iron_components) > 0
            add!(model[id], :Iron => SumReducer(iron_components))
            add!(model[id], :main => SumReducer([:Continuum, :Iron]))
            evaluate!(model)
            bestfit = fit!(model, only_id=id, source.data, minimizer=mzer); show(source.log, bestfit)
        else
            add!(model[id], :Iron => @expr(m -> [0.]))
            add!(model[id], :main => SumReducer([:Continuum, :Iron]))
        end
        (:ironuv    in keys(model[id]))  &&  freeze(model[id], :ironuv)
        (:ironoptbr in keys(model[id]))  &&  freeze(model[id], :ironoptbr)
        (:ironoptna in keys(model[id]))  &&  freeze(model[id], :ironoptna)
    end
    evaluate!(model)

    # Add emission lines
    line_names = [collect(keys(source.line_names[id])) for id in 1:Nspec]
    line_groups = [unique(collect(values(source.line_names[id]))) for id in 1:Nspec]
    println(source.log, "\nFit known emission lines...")
    for id in 1:Nspec
        λ = source.domain[id][:]
        resid = source.data[id].val - model[id]()  # will be used to guess line normalization
        add!(model[id], source.line_comps[id])
        for (group, lnames) in QSFit.invert_dictionary(source.line_names[id])
            add!(model[id], group => SumReducer(lnames))
        end
        add!(model[id], :main => SumReducer([:Continuum, :Iron, line_groups[id]...]))

        if haskey(model[id], :MgII_2798)
            model[id][:MgII_2798].voff.low  = -1000
            model[id][:MgII_2798].voff.high =  1000
        end
        if haskey(model[id], :OIII_5007_bw)
            model[id][:OIII_5007_bw].fwhm.val  = 500
            model[id][:OIII_5007_bw].fwhm.low  = 1e2
            model[id][:OIII_5007_bw].fwhm.high = 1e3
            model[id][:OIII_5007_bw].voff.low  = 0
            model[id][:OIII_5007_bw].voff.high = 2e3
        end
        for cname in line_names[id]
            model[id][cname].norm_integrated = source.options[:norm_integrated]
        end

        # Guess values
        evaluate!(model)
        y = source.data[id].val - model[id]()
        for cname in line_names[id]
            c = model[id][cname]
            resid_at_line = Spline1D(λ, resid, k=1, bc="nearest")(c.center.val)
            c.norm.val *= abs(resid_at_line) / maximum(model[id](cname))

            # If instrumental broadening is not used and the line profile
            # is a Gaussian one take spectral resolution into account.
            # This is significantly faster than convolving with an
            # instrument response but has some limitations:
            # - works only with Gaussian profiles;
            # - all components must be additive (i.e. no absorptions)
            # - further narrow components (besides known emission lines)
            #   will not be corrected for instrumental resolution
            if !source.options[:instr_broadening]
                if isa(c, QSFit.SpecLineGauss)
                    c.spec_res_kms = source.spectra[id].resolution
                else
                    println(source.log, "Line $cname is not a Gaussian profile: Can't take spectral resolution into account")
                end
            end
        end

        # Patch parameters
        if  haskey(model[id], :OIII_4959)  &&
            haskey(model[id], :OIII_5007)
            model[id][:OIII_4959].norm.fixed = true
            model[id][:OIII_4959].voff.fixed = true
            @patch! model[id] m -> begin
                m[:OIII_4959].norm = m[:OIII_5007].norm / 3
                m[:OIII_4959].voff = m[:OIII_5007].voff
            end
        end
        if  haskey(model[id], :OIII_5007_bw)  &&
            haskey(model[id], :OIII_5007)
            @patch! model[id] m -> begin
                m[:OIII_5007_bw].voff += m[:OIII_5007].voff
                m[:OIII_5007_bw].fwhm += m[:OIII_5007].fwhm
            end
        end
        if  haskey(model[id], :OI_6300)  &&
            haskey(model[id], :OI_6364)
            # model[id][:OI_6300].norm.fixed = true
            model[id][:OI_6300].voff.fixed = true
            @patch! model[id] m -> begin
                # m[:OI_6300].norm = m[:OI_6364].norm / 3
                m[:OI_6300].voff = m[:OI_6364].voff
            end
        end
        if  haskey(model[id], :NII_6549)  &&
            haskey(model[id], :NII_6583)
            # model[id][:NII_6549].norm.fixed = true
            model[id][:NII_6549].voff.fixed = true
            @patch! model[id] m -> begin
                # m[:NII_6549].norm = m[:NII_6583].norm / 3
                m[:NII_6549].voff = m[:NII_6583].voff
            end
        end
        if  haskey(model[id], :SII_6716)  &&
            haskey(model[id], :SII_6731)
            # model[id][:SII_6716].norm.fixed = true
            model[id][:SII_6716].voff.fixed = true
            @patch! model[id] m -> begin
                # m[:SII_6716].norm = m[:SII_6731].norm / 1.5
                m[:SII_6716].voff = m[:SII_6731].voff
            end
        end

        if  haskey(model[id], :na_Ha)  &&
            haskey(model[id], :na_Hb)
            model[id][:na_Hb].voff.fixed = true
            @patch! model[id] m -> m[:na_Hb].voff = m[:na_Ha].voff
        end

        # The following are required to avoid degeneracy with iron
        # template
        if  haskey(model[id], :Hg)  &&
            haskey(model[id], :br_Hb)
            model[id][:Hg].voff.fixed = true
            model[id][:Hg].fwhm.fixed = true
            @patch! model[id] m -> begin
                m[:Hg].voff = m[:br_Hb].voff
                m[:Hg].fwhm = m[:br_Hb].fwhm
            end
        end
        if  haskey(model[id], :br_Hg)  &&
            haskey(model[id], :br_Hb)
            model[id][:br_Hg].voff.fixed = true
            model[id][:br_Hg].fwhm.fixed = true
            @patch! model[id] m -> begin
                m[:br_Hg].voff = m[:br_Hb].voff
                m[:br_Hg].fwhm = m[:br_Hb].fwhm
            end
        end
        if  haskey(model[id], :na_Hg)  &&
            haskey(model[id], :na_Hb)
            model[id][:na_Hg].voff.fixed = true
            model[id][:na_Hg].fwhm.fixed = true
            @patch! model[id] m -> begin
                m[:na_Hg].voff = m[:na_Hb].voff
                m[:na_Hg].fwhm = m[:na_Hb].fwhm
            end
        end

        # Ensure luminosity at peak of the broad base component is
        # smaller than the associated broad component:
        if  haskey(model[id], :br_Hb)  &&
            haskey(model[id], :bb_Hb)
            model[id][:bb_Hb].norm.high = 1
            model[id][:bb_Hb].norm.val  = 0.5
            @patch! model[id] m -> m[:bb_Hb].norm *= m[:br_Hb].norm / m[:br_Hb].fwhm * m[:bb_Hb].fwhm
        end
        if  haskey(model[id], :br_Ha)  &&
            haskey(model[id], :bb_Ha)
            model[id][:bb_Ha].norm.high = 1
            model[id][:bb_Ha].norm.val  = 0.5
            @patch! model[id] m -> m[:bb_Ha].norm *= m[:br_Ha].norm / m[:br_Ha].fwhm * m[:bb_Ha].fwhm
        end

        bestfit = fit!(model, only_id=id, source.data, minimizer=mzer); show(source.log, bestfit)

        for lname in line_names[id]
            freeze(model[id], lname)
        end
    end

    # Add unknown lines
    println(source.log, "\nFit unknown emission lines...")
    if source.options[:n_unk] > 0
        for id in 1:Nspec
            tmp = OrderedDict{Symbol, GFit.AbstractComponent}()
            for j in 1:source.options[:n_unk]
                tmp[Symbol(:unk, j)] = line_component(TRecipe, QSFit.UnkLine(5e3))
                tmp[Symbol(:unk, j)].norm_integrated = source.options[:norm_integrated]
            end
            add!(model[id], :UnkLines => SumReducer(collect(keys(tmp))), tmp)
            add!(model[id], :main => SumReducer([:Continuum, :Iron, line_groups[id]..., :UnkLines]))
            evaluate!(model)
            for j in 1:source.options[:n_unk]
                freeze(model[id], Symbol(:unk, j))
            end
        end
    else
        # Here we need a :UnkLines reducer, even when n_unk is 0
        for id in 1:Nspec
            add!(model[id], :UnkLines => @expr(m -> [0.]))
            add!(model[id], :main => SumReducer([:Continuum, :Iron, line_groups[id]...]))
        end
    end
    evaluate!(model)

    # Set "unknown" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    for id in 1:Nspec
        λ = source.domain[id][:]
        λunk = Vector{Float64}()
        while true
            (length(λunk) >= source.options[:n_unk])  &&  break
            evaluate!(model)
            Δ = (source.data[id].val - model[id]()) ./ source.data[id].unc

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
            model[id][cname].norm.val = 1.
            model[id][cname].center.val  = λ[iadd]
            model[id][cname].center.low  = λ[iadd] - λ[iadd]/10. # allow to shift 10%
            model[id][cname].center.high = λ[iadd] + λ[iadd]/10.

            thaw(model[id], cname)
            bestfit = fit!(model, only_id=id, source.data, minimizer=mzer); show(source.log, bestfit)
            freeze(model[id], cname)
        end
    end
    evaluate!(model)

    # ----------------------------------------------------------------
    # Last run with all parameters free
    println(source.log, "\nLast run with all parameters free...")
    for id in 1:Nspec
        thaw(model[id], :qso_cont)
        (:galaxy in keys(model[id]))        &&  thaw(model[id], :galaxy)
        (:balmer in keys(model[id]))        &&  thaw(model[id], :balmer)
        (:ironuv    in keys(model[id]))     &&  thaw(model[id], :ironuv)
        (:ironoptbr in keys(model[id]))     &&  thaw(model[id], :ironoptbr)
        (:ironoptna in keys(model[id]))     &&  thaw(model[id], :ironoptna)

        for lname in line_names[id]
            thaw(model[id], lname)
        end
        for j in 1:source.options[:n_unk]
            cname = Symbol(:unk, j)
            if model[id][cname].norm.val > 0
                thaw(model[id], cname)
            else
                freeze(model[id], cname)
            end
        end
    end
    bestfit = fit!(model, source.data, minimizer=mzer)

    # Disable "unknown" lines whose normalization uncertainty is larger
    # than 3 times the normalization
    needs_fitting = false
    for id in 1:Nspec
        for ii in 1:source.options[:n_unk]
            cname = Symbol(:unk, ii)
            isfixed(model[id], cname)  &&  continue
            if bestfit[id][cname].norm.val == 0.
                freeze(model[id], cname)
                needs_fitting = true
                println(source.log, "Disabling $cname (norm. = 0)")
            elseif bestfit[id][cname].norm.unc / bestfit[id][cname].norm.val > 3
                model[id][cname].norm.val = 0.
                freeze(model[id], cname)
                needs_fitting = true
                println(source.log, "Disabling $cname (unc. / norm. > 3)")
            end
        end
    end
    if needs_fitting
        println(source.log, "\nRe-run fit...")
        bestfit = fit!(model, source.data, minimizer=mzer)
    end

    println(source.log)
    show(source.log, bestfit)

    elapsed = time() - elapsed
    println(source.log, "\nElapsed time: $elapsed s")
    QSFit.close_log(source)

    return reduce(source, model, bestfit)
end
