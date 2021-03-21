calibsum(calib, args...) = calib .* (.+(args...))

function multi_fit(source::QSO{TRecipe}; ref_id=1) where TRecipe <: DefaultRecipe
    Nspec = length(source.domain)

    elapsed = time()
    mzer = GFit.cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # Arrays containing best fit values to be constrained across epochs
    galaxy_best = Vector{Float64}()
    galaxy_unc  = Vector{Float64}()
    OIII_best = Vector{Float64}()
    OIII_unc  = Vector{Float64}()

    # Initialize components and guess initial values
    println(source.log, "\nFit continuum components...")
    preds = Vector{Prediction}()
    for id in 1:Nspec
        λ = source.domain[id][:]
        pred = Prediction(source.domain[id], :Continuum => Reducer(sum, [:qso_cont]),
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
        if source.options[:use_host_template]
            add!(model[id], :Continuum => Reducer(sum, [:qso_cont, :galaxy]),
                 :galaxy => QSFit.hostgalaxy(source.options[:host_template]))
            model[id][:galaxy].norm.val = Spline1D(λ, source.data[id].val, k=1, bc="error")(5500.)
        end

        # Balmer continuum and pseudo-continuum
        if source.options[:use_balmer]
            tmp = [:qso_cont, :balmer]
            (:galaxy in keys(model[id]))  &&  push!(tmp, :galaxy)
            add!(model[id], :Continuum => Reducer(sum, tmp),
                 :balmer => QSFit.balmercont(0.1, 0.5))
            c = model[id][:balmer]
            c.norm.val  = 0.1
            c.norm.fixed = false
            c.norm.high = 0.5
            c.ratio.val = 0.5
            c.ratio.fixed = false
            c.ratio.low  = 0.3
            c.ratio.high = 1
            patch!(model) do m
                m[id][:balmer].norm *= m[id][:qso_cont].norm
            end
        end
    end

    for id in 1:Nspec
        bestfit = fit!(model, only_id=id, source.data, minimizer=mzer);  show(source.log, bestfit)
        if :galaxy in keys(model[id])
            push!(galaxy_best, bestfit[id][:galaxy].norm.val)
            push!(galaxy_unc , bestfit[id][:galaxy].norm.unc)
        end
    end

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
            add!(model[id], :Iron => Reducer(sum, iron_components))
            add!(model[id], :main => Reducer(sum, [:Continuum, :Iron]))
            evaluate!(model)
            bestfit = fit!(model, only_id=id, source.data, minimizer=mzer); show(source.log, bestfit)
        else
            add!(model[id], :Iron => Reducer(() -> [0.], Symbol[]))
            add!(model[id], :main => Reducer(sum, [:Continuum, :Iron]))
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
            add!(model[id], group => Reducer(sum, lnames))
        end
        add!(model[id], :main => Reducer(sum, [:Continuum, :Iron, line_groups[id]...]))

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
            model[id][:OIII_4959].voff.fixed = true
            patch!(model) do m
                m[id][:OIII_4959].voff = m[id][:OIII_5007].voff
            end
        end
        if  haskey(model[id], :NII_6549)  &&
            haskey(model[id], :NII_6583)
            model[id][:NII_6549].voff.fixed = true
            patch!(model) do m
                m[id][:NII_6549].voff = m[id][:NII_6583].voff
            end
        end
        if  haskey(model[id], :OIII_5007_bw)  &&
            haskey(model[id], :OIII_5007)
            patch!(model) do m
                m[id][:OIII_5007_bw].voff += m[id][:OIII_5007].voff
                m[id][:OIII_5007_bw].fwhm += m[id][:OIII_5007].fwhm
            end
        end

        if  haskey(model[id], :br_Hb)  &&
            haskey(model[id], :bb_Hb)
            # Ensure luminosity at peak of the broad base component is
            # smaller than the associated broad component:
            model[id][:bb_Hb].norm.high = 1
            model[id][:bb_Hb].norm.val  = 0.5
            patch!(model) do m
                m[id][:bb_Hb].norm *= m[id][:br_Hb].norm / m[id][:br_Hb].fwhm * m[id][:bb_Hb].fwhm
            end
        end

        if  haskey(model[id], :br_Ha)  &&
            haskey(model[id], :bb_Ha)
            # Ensure luminosity at peak of the broad base component is
            # smaller than the associated broad component:
            model[id][:bb_Ha].norm.high = 1
            model[id][:bb_Ha].norm.val  = 0.5
            patch!(model) do m
                m[id][:bb_Ha].norm *= m[id][:br_Ha].norm / m[id][:br_Ha].fwhm * m[id][:bb_Ha].fwhm
            end
        end

        #=
        model[id][:br_Hb].voff.fixed = 1
        model[id][:br_Hb].fwhm.fixed = 1
        patch!(model) do m
            m[id][:br_Hb].voff = m[id][:br_Ha].voff
            m[id][:br_Hb].fwhm = m[id][:br_Ha].fwhm
        end
        =#

        bestfit = fit!(model, only_id=id, source.data, minimizer=mzer); show(source.log, bestfit)
        push!(OIII_best, bestfit[id][:OIII_5007].norm.val)
        push!(OIII_unc , bestfit[id][:OIII_5007].norm.unc)

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
            add!(model[id], :UnkLines => Reducer(sum, collect(keys(tmp))), tmp)
            add!(model[id], :main => Reducer(sum, [:Continuum, :Iron, line_groups[id]..., :UnkLines]))
            evaluate!(model)
            for j in 1:source.options[:n_unk]
                freeze(model[id], Symbol(:unk, j))
            end
        end
    else
        # Here we need a :UnkLines reducer, even when n_unk is 0
        for id in 1:Nspec
            add!(model[id], :UnkLines => Reducer(() -> [0.], Symbol[]))
            add!(model[id], :main => Reducer(sum, [:Continuum, :Iron, line_groups[id]...]))
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
    # Constrain component normalization across epochs.  Note:
    # reference spectrum must have reliable estimation of all common
    # components
    if :galaxy in keys(model[ref_id])
        rval = [galaxy_best[ref_id], OIII_best[ref_id]]
        runc = [galaxy_unc[ ref_id], OIII_unc[ ref_id]]
    else
        rval = [OIII_best[ref_id]]
        runc = [OIII_unc[ ref_id]]
    end
    @assert all(isfinite.(rval))
    @assert all(rval .!= 0)
    @assert all(isfinite.(runc))
    @assert all(runc .!= 0)

    # Estimate calibration in all epochs w.r.t. reference epoch
    for id in 1:Nspec
        if id != ref_id
            if :galaxy in keys(model[id])
            val = [galaxy_best[id], OIII_best[id]]
            unc = [galaxy_unc[ id], OIII_unc[ id]]
            else
                val = [OIII_best[id]]
                unc = [OIII_unc[ id]]
            end
            j = findall((val .!= 0)  .&  (unc .!= 0))

            R = val[j] ./ rval[j]
            R_unc = (unc[j] ./ val[j] .+ runc[j] ./ rval[j]) .* R
            ratio = sum(R ./ R_unc) ./ sum(1 ./ R_unc)
            for cname in keys(model[id].cevals)
                (string(cname)[1:3] == "bb_")  &&  continue  # this is a patched param, no need to apply calib. factor
                (cname == :balmer)  &&  continue
                if :norm in propertynames(model[id][cname])
                    model[id][cname].norm.val /= ratio
                end
            end
            add!(model[id],
                 :main => Reducer(calibsum, [:calib, :Continuum, :Iron, line_groups[id]..., :UnkLines]),
                 :calib => ratio)
        end
    end
    evaluate!(model)

    for id in 1:Nspec
        if id != ref_id
            (:galaxy in keys(model[id]))  &&  (model[id][:galaxy].norm.fixed = true)
            model[id][:OIII_5007].norm.fixed = true
            model[id][:OIII_5007].fwhm.fixed = 1
            model[id][:OIII_5007].voff.fixed = 1
            patch!(model) do m
                (:galaxy in keys(model[id]))  &&  (m[id][   :galaxy].norm = m[ref_id][   :galaxy].norm)
                m[id][:OIII_5007].norm = m[ref_id][:OIII_5007].norm
                m[id][:OIII_5007].fwhm = m[ref_id][:OIII_5007].fwhm
                m[id][:OIII_5007].voff = m[ref_id][:OIII_5007].voff
            end
        end
    end
    evaluate!(model)

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
        if id != ref_id
            thaw(model[id], :calib)  # parameter is fixed in preds[ref_id]
        end
    end
    bestfit = fit!(model, source.data, minimizer=mzer); show(source.log, bestfit)

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
        bestfit = fit!(model, source.data, minimizer=mzer); show(source.log, bestfit)
    end

    println(source.log, "\nFinal model and bestfit:")
    show(source.log, model)
    println(source.log)
    show(source.log, bestfit)

    elapsed = time() - elapsed
    println(source.log, "\nElapsed time: $elapsed s")
    QSFit.close_log(source)

    QSFit.populate_metadata!(source, model)
    return (model, bestfit)
end
