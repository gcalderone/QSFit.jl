abstract type DefaultRecipe <: AbstractRecipe end

function default_options(::Type{T}) where T <: DefaultRecipe
    out = OrderedDict{Symbol, Any}()
    out[:wavelength_range] = [1215, 7.3e3]
    out[:min_spectral_coverage] = Dict(:default => 0.6,
                                       :ironuv  => 0.3,
                                       :ironopt => 0.3)
    out[:skip_lines] = [:OIII_5007_bw] #, Ha_bb
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
    list = [
        CombinedType( new_transition(name=:Lyb         , 1026.0   ), [BroadType, NarrowType]),
        NarrowType(   new_transition(name=:OV          , 1213.8  )),  # Ferland+92, Shields+95
        CombinedType( new_transition(name=:Lya         , 1215.24  ), [BroadType, NarrowType]),
        NarrowType(   new_transition(name=:OV          , 1218.3  )),  # Ferland+92, Shields+95
        NarrowType(   new_transition(name=:NV_1241     , 1240.81 )),
        BroadType(    new_transition(name=:OI_1306     , 1305.53 )),
        BroadType(    new_transition(name=:CII_1335    , 1335.31 )),
        BroadType(    new_transition(name=:SiIV_1400   , 1399.8  )),
        CombinedType( new_transition(name=:CIV_1549    , 1549.48  ), [BroadType, NarrowType]),
        BroadType(    new_transition(name=:HeII        , 1640.4  )),
        BroadType(    new_transition(name=:OIII        , 1665.85 )),
        BroadType(    new_transition(name=:AlIII       , 1857.4  )),
        BroadType(    new_transition(name=:CIII_1909   , 1908.734)),
        BroadType(    new_transition(name=:CII         , 2326.0  )),
        BroadType(    new_transition(name=:F2420       , 2420.0  )),
        CombinedType( new_transition(name=:MgII_2798   , 2799.117 ), [BroadType, NarrowType]),
        NarrowType(   new_transition(name=:NeVN        , 3346.79 )),
        NarrowType(   new_transition(name=:NeVI_3426   , 3426.85 )),
        NarrowType(   new_transition(name=:OII_3727    , 3729.875)),
        NarrowType(   new_transition(name=:NeIII_3869  , 3869.81 )),
        BroadType(    new_transition(name=:Hd          , 4102.89 )),
        BroadType(    new_transition(name=:Hg          , 4341.68 )),
        NarrowType(   new_transition(name=:OIII_4363   , 4363.00 )),  # TODO: Check wavelength is correct
        BroadType(    new_transition(name=:HeII        , 4686.   )),
        CombinedType( new_transition(name=:Hb          , 4862.68  ), [BroadType, NarrowType]),
        NarrowType(   new_transition(name=:OIII_4959   , 4960.295)),
        NarrowType(   new_transition(name=:OIII_5007   , 5008.240)),
        NarrowType(   new_transition(name=:OIII_5007_bw, 5008.240)),
        BroadType(    new_transition(name=:HeI_5876    , 5877.30 )),
        NarrowType(   new_transition(name=:OI_6300     , 6300.00 )),  # TODO: Check wavelength is correct
        NarrowType(   new_transition(name=:OI_6364     , 6364.00 )),  # TODO: Check wavelength is correct
        NarrowType(   new_transition(name=:NII_6549    , 6549.86 )),
        CombinedType( new_transition(name=:Ha          , 6564.61  ), [BroadType, NarrowType, BroadBaseType]),
        NarrowType(   new_transition(name=:NII_6583    , 6585.27 )),
        NarrowType(   new_transition(name=:SII_6716    , 6718.29 )),
        NarrowType(   new_transition(name=:SII_6731    , 6732.67 ))]
    return list
end


line_default_component(source::QSO{T}, ltype::AbstractLineType) where T <: DefaultRecipe =
    SpecLineGauss(transition(ltype.tid).LAMBDA_VAC_ANG)

function line_component(source::QSO{T}, ltype::BroadType) where T <: DefaultRecipe
    comp = line_default_component(source, ltype)
    comp.fwhm.val  = 5e3
    comp.fwhm.low  = 900
    comp.fwhm.high = 1.5e4
    comp.voff.low  = -3e3
    comp.voff.high =  3e3

    if ltype.tid == :MgII_2798
        comp.voff.low  = -1e3
        comp.voff.high =  1e3
    end
    return Dict(ltype.tid => LineComponent(ltype, comp, :BroadLines))
end

function line_component(source::QSO{T}, ltype::NarrowType) where T <: DefaultRecipe
    comp = line_default_component(source, ltype)
    comp.fwhm.val  = 5e2
    comp.fwhm.low  = 100
    comp.fwhm.high = 2e3
    comp.voff.low  = -1e3
    comp.voff.high =  1e3

    if ltype.tid == :OIII_5007_bw
        comp.fwhm.val  = 500
        comp.fwhm.high = 1e3
        comp.voff.low  = 0
        comp.voff.high = 2e3
    end
    return Dict(ltype.tid => LineComponent(ltype, comp, :NarrowLines))
end

function line_component(source::QSO{T}, ltype::BroadBaseType) where T <: DefaultRecipe
    comp = line_default_component(source, ltype)
    comp.fwhm.val  = 2e4
    comp.fwhm.low  = 1e4
    comp.fwhm.high = 3e4
    comp.voff.fixed = true
    return Dict(ltype.tid => LineComponent(ltype, comp, :BroadBaseLines))
end

function line_component(source::QSO{T}, ltype::CombinedType) where T <: DefaultRecipe
    out = OrderedDict{Symbol, LineComponent}()
    for t in ltype.types
        lc = collect(values(line_component(source, t(ltype.tid))))
        @assert length(lc) == 1
        lc = lc[1]
        lc = LineComponent(ltype, lc.comp, lc.group)
        if t == BroadType
            out[Symbol(ltype.tid, :_br)] = lc
        elseif t == NarrowType
            lc.comp.fwhm.high = 1e3
            out[Symbol(ltype.tid, :_na)] = lc
        elseif t == BroadBaseType
            out[Symbol(ltype.tid, :_bb)] = lc
        else
            error("Unsupported line type: $t")
        end
    end
    return out
end

# function line_component(source::QSO{T}, name::Symbol, line::UnkLine) where T <: DefaultRecipe
#     comp = SpecLineGauss(line.λ)
#     comp.norm.val = 0.
#     comp.center.fixed = false
#     comp.center.low = 0
#     comp.center.high = Inf
#     comp.fwhm.val  = 5e3
#     comp.fwhm.low  = 600
#     comp.fwhm.high = 1e4
#     comp.voff.fixed = true
#     return LineComponent(nothing, typeof(line), comp, :UnknownLines)
# end


function PreparedSpectrum(source::QSO{T}; id=1) where T <: DefaultRecipe
    data = deepcopy(source.specs[id])
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
    for line in known_spectral_lines(source)
        (line.tid in source.options[:skip_lines])  &&  continue
        for (lname, lc) in line_component(source, line)
            (lname in source.options[:skip_lines])  &&  continue
            @assert !haskey(lcs, lname)
            (λmin, λmax, coverage) = spectral_coverage(λ .* data.good, data.resolution, lc.comp)
            coverage = round(coverage * 1e3) / 1e3  # keep just 3 significant digits...
            threshold = get(source.options[:min_spectral_coverage], lname, source.options[:min_spectral_coverage][:default])
            print(logio(source), @sprintf("Line %-15s coverage: %4.2f (threshold: %4.2f)", lname, coverage, threshold))
            if coverage < threshold
                print(logio(source), @sprintf("  neglecting range: %10.5g < λ < %10.5g", λmin, λmax))
                ii = findall(λmin .<= λ .< λmax)
                data.good[ii] .= false
            else
                lcs[lname] = lc
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
    lum = Measures(Domain(data.λ[ii]),
                   data.flux[ii] .* dered[ii] .* source.flux2lum .* (1 + source.z),
                   data.err[ ii] .* dered[ii] .* source.flux2lum .* (1 + source.z))

    return PreparedSpectrum(id, data, dom, lum, lcs)
end


function fit!(source::QSO{T}, model::Model, pspec::PreparedSpectrum) where T <: DefaultRecipe
    mzer = GFit.cmpfit()
    mzer.Δfitstat_theshold = 1.e-5
    bestfit = fit!(model, pspec.data, minimizer=mzer)
    show(logio(source), bestfit)
    # @gp (domain(model), pspec.data) model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
    return bestfit
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
    # @gp (domain(model), pspec.data) model
    # printstyled(color=:blink, "Press ENTER to continue..."); readline()
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
        haskey(groups, lc.group)  ||  (groups[lc.group] = Vector{Symbol}())
        push!(groups[lc.group], cname)
    end
    for (group, lnames) in groups
        model[group] = SumReducer(lnames)
    end
    evaluate!(model)
end


function guess_emission_lines_values!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    for group in [:BroadLines, :NarrowLines, :BroadBaseLines]
        found = false
        for (cname, lc) in pspec.lcs
            (lc.group == group)  ||  continue
            guess_norm_factor!(pspec, model, cname)
            found = true
        end
        found  &&  push!(model[:main].list, group)
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

    @patch! model[:Hb_na].voff = model[:Ha_na].voff

    # The following are required to avoid degeneracy with iron
    # template
    @patch! begin
        model[:Hg].voff = model[:Hb_br].voff
        model[:Hg].fwhm = model[:Hb_br].fwhm
    end
    @patch! begin
        model[:Hg_br].voff = model[:Hb_br].voff
        model[:Hg_br].fwhm = model[:Hb_br].fwhm
    end
    @patch! begin
        model[:Hg_na].voff = model[:Hb_na].voff
        model[:Hg_na].fwhm = model[:Hb_na].fwhm
    end

    # Ensure luminosity at peak of the broad base component is
    # smaller than the associated broad component:
    if  haskey(model, :Hb_br)  &&
        haskey(model, :Hb_bb)
        model[:Hb_bb].norm.high = 1
        model[:Hb_bb].norm.val  = 0.5
        @patch! model[:Hb_bb].norm *= model[:Hb_br].norm / model[:Hb_br].fwhm * model[:Hb_bb].fwhm
    end
    if  haskey(model, :Ha_br)  &&
        haskey(model, :Ha_bb)
        model[:Ha_bb].norm.high = 1
        model[:Ha_bb].norm.val  = 0.5
        @patch! model[:Ha_bb].norm *= model[:Ha_br].norm / model[:Ha_br].fwhm * model[:Ha_bb].fwhm
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
        bestfit = fit!(source, model, pspec)
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
    bestfit = fit!(source, model, pspec)
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
        bestfit = fit!(source, model, pspec)
        haskey(model, :ironuv   )  &&  freeze(model, :ironuv)
        haskey(model, :ironoptbr)  &&  freeze(model, :ironoptbr)
        haskey(model, :ironoptna)  &&  freeze(model, :ironoptna)
    end
    evaluate!(model)

    println(logio(source), "\nFit known emission lines...")
    add_emission_lines!(source, pspec, model)
    guess_emission_lines_values!(source, pspec, model)
    add_patch_functs!(source, pspec, model)

    bestfit = fit!(source, model, pspec)
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
    bestfit = fit!(source, model, pspec)

    if neglect_weak_features!(source, pspec, model)
        println(logio(source), "\nRe-run fit...")
        bestfit = fit!(source, model, pspec)
    end

    println(logio(source))
    show(logio(source), bestfit)

    out = QSFit.QSFitResults(source, pspec, model, bestfit)
    elapsed = time() - elapsed
    println(logio(source), "\nElapsed time: $elapsed s")
    QSFit.close_logio(source)
    return out
end
