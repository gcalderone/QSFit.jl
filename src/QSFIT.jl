#module QSFIT

import DataFitting: AbstractDomain, Domain_1D, Domain_2D,
    Parameter, AbstractComponent, AbstractComponentData,
    cdata, evaluate!

using CMPFit, DataFitting, Gnuplot, ReusePatterns
using Serialization, Statistics, DataFrames, FFTW, DelimitedFiles, Interpolations, Printf
using Cosmology, Unitful, UnitfulAstro, FITSIO, Parameters

DataFitting.@enable_CMPFit

include("utils.jl")
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
    good::Vector{Int}
    goodfraction::Float64
    median_flux::Float64
    median_err::Float64
    meta::Dict{Symbol, Any}
    function QSFitData(λ::Vector{T}, flux::Vector{T}, err::Vector{T},
                       igood::Union{Nothing, Vector{Int}}=nothing; label="") where T <: AbstractFloat
        if igood == nothing
            igood = collect(1:length(λ))
        end
        @assert length(λ) == length(flux) == length(err)
        @assert (minimum(igood) >= 1)  &&  (maximum(igood) <= length(λ))
        @assert length(unique(igood)) == length(igood)
        ii = sortperm(λ)
        new(label, λ[ii], flux[ii], err[ii], sort(igood), length(igood)/length(λ),
            median(flux), median(err), Dict{Symbol, Any}())
    end
end

QSFitData(wave::Vector{Quantity}, flux::Vector{Quantity}, err::Vector{Quantity}, igood=nothing; label="") =
    QSFitData(getproperty.(uconvert.(Ref(unit_λ())   , wave), :val),
              getproperty.(uconvert.(Ref(unit_flux()), flux), :val),
              getproperty.(uconvert.(Ref(unit_flux()), err ), :val),
              igood, label=label)

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
        ld = luminosity_dist(cosmo, float(z))
        #ld = 1 * UnitfulAstro.Gpc
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

    ii = findall(data.λ .> qsfit.options.min_wavelength)
    if length(ii) != length(data.λ)
        mm = fill(false, length(data.λ))
        mm[data.good] .= true
        mm[ii] .= true
        empty!(data.good)
        append!(data.good, findall(mm))
    end

    # Ignore lines on missing data
    let
        # Perform 3-5th test
        mm = fill(false, length(data.λ))
        mm[data.good] .= true
        println(qsfit.log, "Good samples before 3/5th test: ", length(data.good))
        for line in qsfit.lines
            line.enabled  ||  continue
            if isa(line, EmissionLine)
                λrange = line.λ .* (1 .+ [-1, 1] .* line.fwhm / 3.e5)
                δ = (λrange[2] - λrange[1]) / 5.
                count = 0
                for ii in 1:5
                    jj = findall((ii-1)*δ .<=  (data.λ .- λrange[1])  .< ii*δ)
                    (length(jj) > 0)  &&  (count += 1)
                end
                if count < 3
                    println(qsfit.log, "Disabling line: ", line.label)
                    line.enabled = false # Disable line and ignore data
                    ii = findall(λrange[1] .<= data.λ .< λrange[2])
                    mm[ii] .= false
                end
            end
        end
        empty!(data.good)
        append!(data.good, findall(mm))
        println(qsfit.log, "Good samples after  3/5th test: ", length(data.good))
    end
    dom = Domain(data.λ ./ (1 + qsfit.z))
    lum = Measures(data.flux .* qsfit.flux2lum .* (1 + qsfit.z),
                   data.err  .* qsfit.flux2lum .* (1 + qsfit.z))
    push!(qsfit.domain, dom)
    push!(qsfit.data, lum)
end


function run(qsfit::QSFit)
    # try
    elapsed = time()
    @assert length(qsfit.domain) == length(qsfit.data)

    mzer = cmpfit()
    mzer.config.ftol = 1.e-6
    mzer.config.gtol = 1.e-6
    mzer.config.xtol = 1.e-6
    
    # First dataset have a special meaning
    λ = qsfit.domain[1][1]
    L = qsfit.data[1].val

    # Prepare model
    model = Model()

    # List of component names
    cnames = Dict{Symbol,   Vector{Symbol}}()
    cnames[:broadband]    = Vector{Symbol}()
    cnames[:narrow_lines] = Vector{Symbol}()
    cnames[:broad_lines]  = Vector{Symbol}()
    cnames[:unknown]      = Vector{Symbol}()

    # Add continuum components
    let
        l0 = median(λ)
        addcomp!(model, :continuum => powerlaw(l0))
        model.continuum.norm.val = interpol(L, λ, l0)

        model.continuum.alpha.val = -1.5
        model.continuum.alpha.low  = -3
        model.continuum.alpha.high = 1 #--> -3, 1 in frequency
        if qsfit.z < qsfit.options.alpha1_fixed_z
            model.continuum.alpha.val = qsfit.options.alpha1_fixed_value
            model.continuum.alpha.fixed = true # avoid degeneracy with galaxy template
        end
        push!(cnames[:broadband], :continuum)
    end

    # Host galaxy
    let
        # Galaxy component is disabled when z > qsfit.options.galaxy_enabled_z
        if qsfit.options.use_galaxy  &&  (qsfit.z > qsfit.options.galaxy_enabled_z)
            println(qsfit.log, "Galaxy component is disabled")
            qsfit.options.use_galaxy = false
        end

        if qsfit.options.use_galaxy
            addcomp!(model, :galaxy => hostgalaxy(qsfit.options.galaxy_template))
            model.galaxy.norm.val = interpol(L, λ, 5500)

            # If 5500 is outisde the available range compute value at
            # the edge (solve problems for e.g., spec-0411-51817-0198)
            if maximum(λ) < 5500
                model.galaxy.norm.val = interpol(L, λ, maximum(λ))
            end
            (model.galaxy.norm.val < 1e-4)  &&  (model.galaxy.norm.val = 1e-4)
            push!(cnames[:broadband], :galaxy)
        end
    end

    # Add Balmer continuum
    let
        if qsfit.options.use_balmer
            addcomp!(model, :balmer => balmercont(0.1, 0.5))
            if qsfit.z < qsfit.options.balmer_fixed_z
                model.balmer.norm.val   = 0.1
                model.balmer.norm.fixed = false
                model.balmer.norm.low   = 0
                model.balmer.norm.high  = 0.5
                model.balmer.ratio.val   = 0.3
                model.balmer.ratio.fixed = false
                model.balmer.ratio.low   = 0.3
                model.balmer.ratio.high  = 1
            else
                model.balmer.norm.val   = 0.1
                model.balmer.norm.fixed = true
                model.balmer.ratio.val   = 0.3
                model.balmer.ratio.fixed = true
            end
            push!(cnames[:broadband], :balmer)
        end
    end

    # Add domain(s) and expression for the first run
    for ii in 1:length(qsfit.domain)
        add_dom!(model, qsfit.domain[ii])
    end
    addexpr!(model, :broadband, Meta.parse(join(cnames[:broadband], ".+")))
    bestfit = fit!(model, qsfit.data, minimizer=mzer)

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
    model.continuum.fixed = true
    qsfit.options.use_galaxy  &&  (model.galaxy.fixed = true)
    qsfit.options.use_balmer  &&  (model.balmer.fixed = true)

    # Add Iron
    let
        if qsfit.options.use_ironuv  &&  (minimum(λ) > 2900)
            println(qsfit.log, "Iron UV is disabled")
            qsfit.options.use_ironuv = false
        end
        if qsfit.options.use_ironuv
            addcomp!(model, :ironuv => ironuv(3000))
            push!(cnames[:broadband], :ironuv)
        end
        if qsfit.options.use_ironopt
            addcomp!(model, :ironoptbr => ironopt_broad(3000))
            addcomp!(model, :ironoptna => ironopt_narrow(500))
            push!(cnames[:broadband], :ironoptbr)
            push!(cnames[:broadband], :ironoptna)
        end
        
        replaceexpr!(model, :broadband, Meta.parse(join(cnames[:broadband], ".+")))
        bestfit = fit!(model, qsfit.data, minimizer=mzer)
        qsfit.options.use_ironuv    &&  (model.ironuv.fixed = true)
        qsfit.options.use_ironopt   &&  (model.ironoptbr.fixed = true)
        qsfit.options.use_ironopt   &&  (model.ironoptna.fixed = true)
    end

    # Add emission lines
    let
        x = model(:domain)
        y = qsfit.data[1].val - model(:broadband)

        for line in qsfit.lines
            if line.enabled
                cname = addEmLine(model, line)
                if     isa(line, NarrowLine)  ||  isa(line, CombinedNarrowLine); push!(cnames[:narrow_lines], cname)
                elseif isa(line,  BroadLine)  ||  isa(line, CombinedBroadLine ); push!(cnames[:broad_lines] , cname)
                else
                    error("Unexpected line type: " * string(typeof(line)))
                end
                yatline = interpol(y, x, line.λ)
                model[cname].norm.val *= abs(yatline) / maxvalue(model[cname])
            end
        end
        addexpr!(model, :narrow_lines, Meta.parse(join(cnames[:narrow_lines], ".+")), cmp=false)
        addexpr!(model, :broad_lines , Meta.parse(join(cnames[:broad_lines] , ".+")), cmp=false)
        setflag!(model, :broadband, false)
        addexpr!(model, :known, :(broadband + narrow_lines + broad_lines))
        bestfit = fit!(model, qsfit.data, minimizer=mzer)

        for line in qsfit.lines
            cname = Symbol(line.label)
            line.enabled  &&  (model[cname].fixed = true)
        end
    end

    # Add unknown lines
    let
        # Set "unknown" line center wavelength where there is a
        # maximum in the fit residuals, and re-run a fit.

        # First add all unknown lines
        for ii in 1:qsfit.options.unkLines
            line = UnknownLine(λ[1])
            line.label = "unk" * string(ii)
            cname = addEmLine(model, line)
            push!(cnames[:unknown], cname)
            model[cname].norm.val = 0
            model[cname].fixed = true
        end
        addexpr!(model, :unknown, Meta.parse(join(cnames[:unknown] , ".+")), cmp=false)
        setflag!(model, :known, false)
        addexpr!(model, :all, :(known .+ unknown))

        # Now set the appropriate wavelength
        iadd = Vector{Int}()
        run = true
        while run
            (length(iadd) >= qsfit.options.unkLines)  &&  break
            run  ||  break
            evaluate!(model)
            Δ = (qsfit.data[1].val - model(:known)) ./ qsfit.data[1].unc
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
            model[cname].norm.val = 1
            model[cname].center.val  = λ[imax]
            model[cname].center.low  = λ[imax] - 100
            model[cname].center.high = λ[imax] + 100            
            model[cname].fixed = false
            bestfit = fit!(model, qsfit.data, minimizer=mzer)
            model[cname].fixed = true
        end
    end

    # Last run with all parameters free
    model.continuum.fixed = false
    qsfit.options.use_galaxy  &&  (model.galaxy.fixed = false)
    qsfit.options.use_balmer  &&  (model.balmer.fixed = false)
    qsfit.options.use_ironuv  &&  (model.ironuv.fixed = false)
    qsfit.options.use_ironopt &&  (model.ironoptbr.fixed = false)
    qsfit.options.use_ironopt &&  (model.ironoptna.fixed = false)
    for cname in cnames[:narrow_lines];  model[cname].fixed = false;  end
    for cname in cnames[:broad_lines] ;  model[cname].fixed = false;  end
    for cname in cnames[:unknown]     ;  model[cname].fixed = false;  end
    bestfit = fit!(model, qsfit.data, minimizer=mzer)
    elapsed = time() - elapsed
    @info elapsed
    return (model, bestfit)
    # catch err
    # finally
    #     if qsfit.log != stdout
    #         close(qsfit.log)
    #     end
    #     return nothing
    # end
end


#qsfit = QSFit("dd", 0.3806, 0.06846);
#adddata!(qsfit, read_sdss_dr10("spec-0752-52251-0323.fits"));

qsfit = QSFit("dd", 0.178064, 0.02);
adddata!(qsfit, read_sdss_dr10("spec-1959-53440-0066.fits"));
(model, res) = run(qsfit);  # 8.0 s
plot(qsfit, model)

# end  # module
