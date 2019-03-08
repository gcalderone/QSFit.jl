#module QSFIT

import DataFitting: AbstractDomain, Domain_1D, Domain_2D,
    Parameter, AbstractComponent, AbstractComponentData,
    cdata, evaluate!

export QSFit, adddata!, read_sdss_dr10, run, plot

using CMPFit, DataFitting, Gnuplot, ReusePatterns
using Serialization, Statistics, DataFrames, FFTW, DelimitedFiles, Interpolations, Printf
using Cosmology, Unitful, UnitfulAstro, FITSIO, Parameters

DataFitting.@enable_CMPFit
DataFitting.showsettings.fixedpars = false
const showstep = true

include("utils.jl")
include("ccm_unred.jl")
include("components/emline.jl")
include("components/powerlaw.jl")
include("components/hostgalaxy.jl")
include("components/ironopt.jl")
include("components/ironuv.jl")
include("components/balmercont.jl")
#include("components/test_components.jl")

qsfit_cosmology() = cosmology(h=0.70, OmegaM=0.3)   #S11

unit_λ() = UnitfulAstro.angstrom
unit_flux() = 1.e-17 * UnitfulAstro.erg / UnitfulAstro.s / UnitfulAstro.cm^2
unit_lum() =  1.e42  * UnitfulAstro.erg / UnitfulAstro.s

struct QSFitData
    label::String
    λ::Vector{Float64}
    flux::Vector{Float64}
    err::Vector{Float64}
    good::Vector{Bool}
    goodfraction::Float64
    median_flux::Float64
    median_err::Float64
    meta::Dict{Symbol, Any}
    function QSFitData(λ::Vector{T}, flux::Vector{T}, err::Vector{T},
                       good::Union{Nothing, Vector{Bool}}=nothing; label="") where T <: AbstractFloat
        if good == nothing
            good = fill(true, length(λ))
        end
        @assert length(λ) == length(flux) == length(err) == length(good)
        @assert issorted(λ)
        igood = findall(good)
        new(label, λ, flux, err, good, length(igood)/length(λ),
            median(flux[igood]), median(err[igood]), Dict{Symbol, Any}())
    end
end

QSFitData(wave::Vector{Quantity}, flux::Vector{Quantity}, err::Vector{Quantity}, good=nothing; label="") =
    QSFitData(getproperty.(uconvert.(Ref(unit_λ())   , wave), :val),
              getproperty.(uconvert.(Ref(unit_flux()), flux), :val),
              getproperty.(uconvert.(Ref(unit_flux()), err ), :val),
              good, label=label)

include("readers.jl")
include("lines.jl")
include("options.jl")
include("reduce.jl")
include("plot.jl")

struct QSFit
    name::String
    z::Float64
    ebv::Float64
    flux2lum::Float64
    lines::Vector{Line}
    log::IO
    options::QSFitOptions
    domain::Vector{DataFitting.Domain_1D}
    data::Vector{DataFitting.Measures_1D}

    function QSFit(name, z, ebv; log="", cosmo=qsfit_cosmology())
        @assert z > 0
        @assert ebv > 0
        ld = luminosity_dist(cosmo, float(z)) # * UnitfulAstro.Gpc
        ld = uconvert(u"cm", ld)
        flux2lum = 4pi * ld^2 * unit_flux() / unit_lum()
        @assert typeof(flux2lum) == Float64
        stream = stdout
        (log != "")  &&  (stream = open(log, "w"))
        return new(string(name), float(z), float(ebv), flux2lum,
                   qsfit_lines(), stream, QSFitOptions(),
                   Vector{DataFitting.Domain_1D}(), Vector{DataFitting.Measures_1D}())
    end
end


function adddata!(qsfit::QSFit, data::QSFitData)
    println(qsfit.log, "New data: " * data.label)
    println(qsfit.log, "  good fraction:: " * string(data.goodfraction))
    if data.goodfraction < 0.5
        @error "Good fraction < 0.5"
    end

    λ = data.λ ./ (1 + qsfit.z)

    ii = findall(λ .<= qsfit.options.min_wavelength)
    data.good[ii] .= false

    # Ignore lines on missing data
    let
        # Perform 3-5th test
        println(qsfit.log, "Good samples before 3/5th test: ", length(findall(data.good)))
        for line in qsfit.lines
            line.enabled  ||  continue
            if isa(line, EmissionLine)
                fwhm = (isa(line, NarrowLine)  ?  1e3  :  1.2e4) / 2.  # TODO: use current line FWHM
                λrange = line.λ .* (1 .+ [-1, 1] .* fwhm / 3.e5)
                δ = (λrange[2] - λrange[1]) / 5.
                count = 0
                for ii in 1:5
                    jj = findall((ii-1)*δ .<=  (λ .- λrange[1])  .< ii*δ)
                    (length(jj) > 0)  &&  (count += 1)
                end
                if count < 3
                    println(qsfit.log, "Disabling line: ", line.label)
                    line.enabled = false # Disable line and ignore data
                    ii = findall(λrange[1] .<= λ .< λrange[2])
                    data.good[ii] .= false
                else
                    println(qsfit.log, " Enabling line: ", line.label)
                end
            end
        end
        println(qsfit.log, "Good samples after  3/5th test: ", length(findall(data.good)))
    end

    dered = ccm_unred([1450, 3000, 5100.], qsfit.ebv)
    println(qsfit.log, "Dereddening factors @ 1450, 3000, 5100 AA: ", dered)
    dered = ccm_unred(data.λ, qsfit.ebv)

    igood = findall(data.good)
    dom = Domain(data.λ[igood] ./ (1 + qsfit.z))
    lum = Measures(data.flux[igood] .* dered[igood] .* qsfit.flux2lum .* (1 + qsfit.z),
                   data.err[ igood] .* dered[igood] .* qsfit.flux2lum .* (1 + qsfit.z))
    push!(qsfit.domain, dom)
    push!(qsfit.data, lum)
end


function run(qsfit::QSFit)
    # try
    elapsed = time()
    @assert length(qsfit.domain) == length(qsfit.data)

    mzer = cmpfit()
    mzer.config.ftol = mzer.config.gtol = mzer.config.xtol = 1.e-6

    # First dataset have a special meaning
    λ = qsfit.domain[1][1]

    # Prepare model
    model = Model()
    for ii in 1:length(qsfit.domain)
        add_dom!(model, qsfit.domain[ii])
    end

    # List of component names
    cnames = Dict{Symbol,   Vector{Symbol}}()
    cnames[:broadband]    = Vector{Symbol}()
    cnames[:narrow_lines] = Vector{Symbol}()
    cnames[:broad_lines]  = Vector{Symbol}()
    cnames[:unknown]      = Vector{Symbol}()

    # AGN continuum
    addcomp!(model, :continuum => powerlaw((minimum(λ) + maximum(λ)) / 2.))
    push!(cnames[:broadband], :continuum)

    # Host galaxy (disabled when z > qsfit.options.galaxy_enabled_z)
    if qsfit.options.use_galaxy  &&  (qsfit.z > qsfit.options.galaxy_enabled_z)
        println(qsfit.log, "Galaxy component is disabled")
        qsfit.options.use_galaxy = false
    end
    if qsfit.options.use_galaxy
        addcomp!(model, :galaxy => hostgalaxy(qsfit.options.galaxy_template))
        push!(cnames[:broadband], :galaxy)
    end

    # Balmer continuum
    if qsfit.options.use_balmer
        addcomp!(model, :balmer => balmercont(0.1, 0.5))
        push!(cnames[:broadband], :balmer)
    end

    # Iron blended lines
    if qsfit.options.use_ironuv  &&  (minimum(λ) > 2900)
        println(qsfit.log, "Iron UV is disabled")
        qsfit.options.use_ironuv = false
    end
    if qsfit.options.use_ironuv
        addcomp!(model, :ironuv => ironuv(3000))
        model.ironuv.norm.val = 0.
        model.ironuv.fixed = true
        push!(cnames[:broadband], :ironuv)
    end
    if qsfit.options.use_ironopt
        addcomp!(model, :ironoptbr => ironopt_broad(3000))
        addcomp!(model, :ironoptna => ironopt_narrow(500))
        model.ironoptbr.norm.val = 0.
        model.ironoptbr.fixed = true
        model.ironoptna.norm.val = 0.
        model.ironoptna.fixed = true
        push!(cnames[:broadband], :ironoptbr)
        push!(cnames[:broadband], :ironoptna)
    end
    addexpr!(model, :broadband, Meta.parse(join(cnames[:broadband], ".+")), cmp=false)

    # Emission lines
    for line in qsfit.lines
        if line.enabled
            cname = addEmLine(model, line)
            if     isa(line, NarrowLine)  ||  isa(line, CombinedNarrowLine); push!(cnames[:narrow_lines], cname)
            elseif isa(line,  BroadLine)  ||  isa(line, CombinedBroadLine ); push!(cnames[:broad_lines] , cname)
            else
                error("Unexpected line type: " * string(typeof(line)))
            end
            model[cname].norm.val = 0.
            model[cname].fixed = true
            if cname == :br_Ha_base
                model[cname].voff.val = 0
                model[cname].voff.fixed = true
            end
            if cname == :na_OIII_4959
                model[cname].voff.expr = "na_OIII_5007_voff"
                model[cname].voff.fixed = true
            end
            if cname == :na_OIII_5007_bw
                model[cname].fwhm.expr = "na_OIII_5007_fwhm + na_OIII_5007_bw_fwhm"
                model[cname].voff.expr = "na_OIII_5007_voff + na_OIII_5007_bw_voff"
            end
        end
    end
    addexpr!(model, :narrow_lines, Meta.parse(join(cnames[:narrow_lines], ".+")), cmp=false)
    addexpr!(model, :broad_lines , Meta.parse(join(cnames[:broad_lines] , ".+")), cmp=false)
    addexpr!(model, :known, :(broadband .+ narrow_lines .+ broad_lines), cmp=false)

    # Unknown lines
    for ii in 1:qsfit.options.unkLines
        line = UnknownLine(λ[1])
        line.label = "unk" * string(ii)
        cname = addEmLine(model, line)
        push!(cnames[:unknown], cname)
        model[cname].norm.val = 0
        model[cname].fixed = true
    end
    addexpr!(model, :unknown, Meta.parse(join(cnames[:unknown] , ".+")), cmp=false)
    addexpr!(model, :all, :(known .+ unknown))

    ##############################################

    # AGN continuum
    let
        L = qsfit.data[1].val
        model.continuum.norm.val = interpol(L, λ, model.continuum.x0.val)
        model.continuum.alpha.val = -1.5
        model.continuum.alpha.low  = -3
        model.continuum.alpha.high = 1 #--> -3, 1 in frequency
        if qsfit.z < qsfit.options.alpha1_fixed_z
            model.continuum.alpha.val = qsfit.options.alpha1_fixed_value
            model.continuum.alpha.fixed = true # avoid degeneracy with galaxy template
        end
    end

    # Host galaxy
    let
        if qsfit.options.use_galaxy
            L = qsfit.data[1].val
            model.galaxy.norm.val = interpol(L, λ, 5500)

            # If 5500 is outisde the available range compute value at
            # the edge (solve problems for e.g., spec-0411-51817-0198)
            if maximum(λ) < 5500
                model.galaxy.norm.val = interpol(L, λ, maximum(λ))
            end
            (model.galaxy.norm.val < 1e-4)  &&  (model.galaxy.norm.val = 1e-4)
        end
    end

    # Balmer continuum
    if qsfit.options.use_balmer
        if qsfit.z < qsfit.options.balmer_fixed_z
            model.balmer.norm.val   = 0.1
            model.balmer.norm.expr = "balmer_norm * Main.interpol1(continuum, domain[1], 3000.) / continuum_norm"
            model.balmer.norm.fixed = false
            model.balmer.norm.low   = 0
            model.balmer.norm.high  = 0.5
            model.balmer.ratio.val   = 0.5
            model.balmer.ratio.fixed = false
            model.balmer.ratio.low   = 0.3
            model.balmer.ratio.high  = 1
        else
            model.balmer.norm.val   = 0.1
            model.balmer.norm.expr = "balmer_norm * Main.interpol1(continuum, domain[1], 3000.)"
            model.balmer.norm.fixed = true
            model.balmer.ratio.val   = 0.3
            model.balmer.ratio.fixed = true
        end
    end
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())

    # Continuum renormalization
    let
        println(qsfit.log, "Cont. norm. (before): ", model.continuum.norm.val)
        check_fraction = -1.
        last_fraction = check_fraction
        yy = qsfit.data[1].val
        ee = qsfit.data[1].unc
        while true
            mm = model(1)

            residuals = (mm .- yy) ./ ee
            check_fraction = length(findall(residuals .< 0)) / length(yy)
            (last_fraction == check_fraction)  &&  break
            last_fraction = check_fraction
            (check_fraction > qsfit.options.cont_negative_fraction)  &&  break

            model.continuum.norm.val *= 0.99
            evaluate!(model)
        end
        println(qsfit.log, "Cont. norm. (after) : ", model.continuum.norm.val)
    end
    evaluate!(model)
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())
    model.continuum.fixed = true
    qsfit.options.use_galaxy  &&  (model.galaxy.fixed = true)
    qsfit.options.use_balmer  &&  (model.balmer.fixed = true)

    if qsfit.options.use_ironuv
        model.ironuv.norm.val = 1.
        model.ironuv.fixed = false
    end
    if qsfit.options.use_ironopt
        model.ironoptbr.norm.val = 0.1
        model.ironoptbr.fixed = false
        model.ironoptna.norm.val = 0.
        model.ironoptna.norm.fixed = true # will be freed during last run
    end
    evaluate!(model)
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())
    qsfit.options.use_ironuv    &&  (model.ironuv.fixed = true)
    qsfit.options.use_ironopt   &&  (model.ironoptbr.fixed = true)
    qsfit.options.use_ironopt   &&  (model.ironoptna.fixed = true)

    # Emission lines
    let
        x = model(:domain)
        y = qsfit.data[1].val - model(:broadband)
        for line in qsfit.lines
            if line.enabled
                yatline = interpol(y, x, line.λ)
                cname = Symbol(line.label)
                model[cname].norm.val = 1.
                model[cname].norm.val = abs(yatline) / maxvalue(model[cname])
                model[cname].fixed = false
            end
        end

        bestfit = fit!(model, qsfit.data, minimizer=mzer)
        showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())
        for line in qsfit.lines
            cname = Symbol(line.label)
            line.enabled  &&  (model[cname].fixed = true)
        end
    end

    # Unknown lines
    let
        # Set "unknown" line center wavelength where there is a
        # maximum in the fit residuals, and re-run a fit.
        iadd = Vector{Int}()
        run = true
        while run
            (length(iadd) >= qsfit.options.unkLines)  &&  break
            run  ||  break
            evaluate!(model)
            Δ = (qsfit.data[1].val - model()) ./ qsfit.data[1].unc
            @assert length(λ) == length(Δ)
            if length(iadd) > 0
                # Avoid considering again the same residuals and those
                # in the 2 neighbour samples
                Δ[vcat([i.+iadd for i in -2:2]...)] .= 0.
            end

            # Do not add lines within 6000 km/s from the edges since
            # these may influence continuum fitting (6000 km/s / speed
            # of light = 0.02).
            Δ[findall((λ .< minimum(λ)*1.02)  .|
                      (λ .> maximum(λ)*0.98))] .= 0.
            imax = argmax(Δ)
            if Δ[imax] <= 0
                printf(qsfit.log, "No residual is greater than 0, skip searching further residuals.")
                break
            end
            push!(iadd, imax)
            println(qsfit.log, "Adding \"unknown\" emission line at " * string(λ[imax]))
            cname = Symbol("unk", length(iadd))
            model[cname].norm.val = 1.
            model[cname].center.val  = λ[imax]
            model[cname].center.low  = λ[imax] - λ[imax]/10. # allow to move 10%
            model[cname].center.high = λ[imax] + λ[imax]/10.
            model[cname].fixed = false
            if false
                plot(qsfit, model)
                @gp :-        λ[imax]*[1,1] [0, 0.6] "w l lw 2"
                @gp :- :resid λ[imax]*[1,1] [-10, 10] "w l lw 2"
                @gp :- :resid [λ[1], λ[end]] Δ[imax]*[1, 1] "w l lw 2"
                sleep(0.3)
            end
            bestfit = fit!(model, qsfit.data, minimizer=mzer)
            model[cname].fixed = true
        end
    end
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())

    # Last run with all parameters free
    model.continuum.fixed = false
    qsfit.options.use_galaxy  &&  (model.galaxy.fixed = false)
    qsfit.options.use_balmer  &&  (model.balmer.fixed = false)
    qsfit.options.use_ironuv  &&  (model.ironuv.fixed = false)
    qsfit.options.use_ironopt &&  (model.ironoptbr.fixed = false)
    qsfit.options.use_ironopt &&  (model.ironoptna.fixed = false)
    qsfit.options.use_ironopt &&  (model.ironoptna.norm.fixed = false)

    for cname in cnames[:narrow_lines];  model[cname].fixed = false;  end
    for cname in cnames[:broad_lines] ;  model[cname].fixed = false;  end
    for cname in cnames[:unknown]     ;  model[cname].fixed = false;  end
    @info "last run"
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())

    elapsed = time() - elapsed
    @info elapsed
    showstep  &&   (show(model); show(bestfit); plot(qsfit, model); readline())
    return (model, bestfit)
    # catch err
    # finally
    #     if qsfit.log != stdout
    #         close(qsfit.log)
    #     end
    #     return nothing
    # end
end

#end  # module
