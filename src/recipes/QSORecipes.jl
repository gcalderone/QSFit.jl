module QSORecipes

using Printf, DataStructures, Statistics, Dates
using Dierckx
using ..QSFit, ..QSFit.ATL, GModelFit

import GModelFit: domain
import QSFit: init_recipe!, preprocess_spec!, SpecLineSet, add_line!, line_group, line_suffix, line_component, analyze, postanalysis, Results

abstract type QSOGeneric <: AbstractRecipe end

getmodel( fp::GModelFit.FitProblem, ith::Int) = fp.multi.v[ith].model
getdomain(fp::GModelFit.FitProblem, ith::Int) = fp.multi.v[ith].domain
getdata(  fp::GModelFit.FitProblem, ith::Int) = fp.data[ith]
function geteval(fp::GModelFit.FitProblem, ith::Int, cname=nothing)
    GModelFit.scan_model!(fp.multi)
    GModelFit.update_eval!(fp.multi)
    if isnothing(cname)
        return GModelFit.last_eval(fp.multi, ith)
    end
    return GModelFit.last_eval(fp.multi, ith, cname)
end
nmodels(fp::GModelFit.FitProblem) = length(fp.multi)

line_component(recipe::CRecipe{<: QSOGeneric}, center::Float64) = recipe.line_component(center)

abstract type BlueWing <: QSFit.LineTemplate end
line_suffix(recipe::CRecipe{<: QSOGeneric}, ::Type{BlueWing}) = :_bw
line_group( recipe::CRecipe{<: QSOGeneric}, ::Type{BlueWing}) = :NarrowLines

function line_component(recipe::CRecipe{<: QSOGeneric}, tid::Val{:OIII_5007}, ::Type{<: BlueWing})
    @track_recipe
    comp = line_component(recipe, tid, NarrowLine)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 0, 3e3, 5e3
    comp.voff.low, comp.voff.val, comp.voff.high = 0, 0, 2e3
    comp.voff.patch = @fd (m, v) -> v + m[:OIII_5007].voff
    comp.fwhm.patch = @fd (m, v) -> v + m[:OIII_5007].fwhm
    return comp
end

function line_component(recipe::CRecipe{<: QSOGeneric}, tid::Val{:OIII_4959}, ::Type{<: BlueWing})
    @track_recipe
    comp = line_component(recipe, tid, NarrowLine)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 0, 3e3, 5e3
    comp.voff.low, comp.voff.val, comp.voff.high = 0, 0, 2e3
    comp.voff.patch = @fd (m, v) -> v + m[:OIII_4959].voff
    comp.fwhm.patch = @fd (m, v) -> v + m[:OIII_4959].fwhm
    return comp
end


function init_recipe!(recipe::CRecipe{T}) where T <: QSOGeneric
    @track_recipe
    @invoke init_recipe!(recipe::CRecipe{<: AbstractRecipe})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.min_spectral_coverage = Dict(:default => 0.6)

    recipe.host_template = Dict(:library=>"swire", :template=>"Ell5", :ref_wavelength => 5500.) # wavelength is in Angstrom
    recipe.use_host_template = true
    recipe.host_template_range = [4000., 7000.]

    recipe.n_nuisance = 10
    recipe.nuisance_avoid = [4863 .+ [-1,1] .* 50,    # Angstrom
                             6565 .+ [-1,1] .* 150]
    recipe.nuisance_maxoffset_from_guess = 1e3  # km/s
    recipe.line_component = QSFit.SpecLineGauss
    recipe.solver = GModelFit.cmpfit()
    recipe.solver.config.ftol = 1.e-6
end


function preprocess_spec!(recipe::CRecipe{T}, spec::Spectrum) where T <: QSOGeneric
    @track_recipe
    @invoke preprocess_spec!(recipe::CRecipe{<: AbstractRecipe}, spec)
    (:specs in propertynames(recipe))   ||  (recipe.specs = Vector{Spectrum}())
    push!(recipe.specs, spec)  # save for later use

    spec.good[findall(spec.x .< recipe.wavelength_range[1])] .= false
    spec.good[findall(spec.x .> recipe.wavelength_range[2])] .= false

    #= Emission line are localized features whose parameter can be
    reliably estimated only if there are sufficient samples to
    constrain the corresponding parameters.  If data coverage is
    not sufficient the component should not be added to the model,
    and corresponding spectral samples should be ignored to avoid
    worsening the fit due to missing model components. =#
    println("Good samples before line coverage filter: ", count(spec.good) , " / ", length(spec.good))

    # Collect line components
    lines = lines_dict(recipe)
    for loop in 1:2
        # The second pass is required to neglect lines whose coverage
        # has been affected by the spectral samples discarded in the
        # first iteration.
        if loop == 2
            println()
            println("Updated coverage:")
        end
        toDelete = Symbol[]
        for (cname, line) in lines
            threshold = get(recipe.min_spectral_coverage, cname, recipe.min_spectral_coverage[:default])
            (λmin, λmax, coverage) = QSFit.spectral_coverage(spec.x[findall(spec.good)],
                                                             spec.resolution, line.comp)
            @printf("Line %-15s coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax)
            if coverage < threshold
            @printf(", threshold is < %5.3f, neglecting...", threshold)
                spec.good[λmin .<= spec.x .< λmax] .= false
                push!(toDelete, cname)
            end
            println()
        end
        delete!.(Ref(lines), toDelete)
        !(:OIII_4959 in keys(lines))  &&  (:OIII_4959_bw in keys(lines))  &&  delete!(lines, :OIII_4959_bw)
        !(:OIII_5007 in keys(lines))  &&  (:OIII_5007_bw in keys(lines))  &&  delete!(lines, :OIII_5007_bw)
    end
    println("Good samples after line coverage filter: ", count(spec.good) , " / ", length(spec.good))

    # Sort lines according to center wavelength, and save list in recipe
    QSFit.sort_by_wavelength!(lines)
    spec.meta[:lines] = lines
end


function fit!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    geteval(fp, 1)  # ensure model is updated
    (GModelFit.nfree(fp) == 0)  &&  (return nothing)
    bestfit, fsumm = GModelFit.fit!(fp, recipe.solver)
    show(fsumm)
    return bestfit, fsumm
end


function add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_qso_continuum!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
    #=
    for i in 2:nmodels(fp)
        m1 = getmodel(fp, 1)
        mi = getmodel(fp, i)
        if haskey(mi, :QSOcont)  &&  haskey(m1, :QSOcont)
            mi[:QSOcont].norm.mpatch  = @fd m -> m[1][:QSOcont].norm
            mi[:QSOcont].x0.mpatch    = @fd m -> m[1][:QSOcont].x0
            mi[:QSOcont].alpha.mpatch = @fd m -> m[1][:QSOcont].alpha
        end
    end
    =#
end
function add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    λ = coords(domain(getdata(fp, ith)))
    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(getdata(fp, ith))) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    getmodel(fp, ith)[:QSOcont] = comp
    push!(getmodel(fp, ith)[:Continuum].list, :QSOcont)
end


function add_host_galaxy!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_host_galaxy!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
    for i in 2:nmodels(fp)
        m1 = getmodel(fp, 1)
        mi = getmodel(fp, i)
        if haskey(mi, :Galaxy)  &&  haskey(m1, :Galaxy)
            mi[:Galaxy].norm.mpatch = @fd m -> m[1][:Galaxy].norm
        end
    end
end
function add_host_galaxy!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    λ = coords(getdomain(fp, ith))
    if recipe.use_host_template  &&
        (recipe.host_template_range[1] .< maximum(λ))  &&
        (recipe.host_template_range[2] .> minimum(λ))
        getmodel(fp, ith)[:Galaxy] = QSFit.hostgalaxy(recipe.host_template[:template],
                                                library=recipe.host_template[:library],
                                                refwl=recipe.host_template[:ref_wavelength])
        push!(getmodel(fp, ith)[:Continuum].list, :Galaxy)

        # Split total flux between continuum and host galaxy
        refwl = recipe.host_template[:ref_wavelength]
        vv = Dierckx.Spline1D(λ, values(getdata(fp, ith)), k=1, bc="extrapolate")(refwl)
        @assert !isnan(vv) "Predicted L_λ at $(refwl)A is NaN"
        if vv <= 0
            @warn "Predicted L_λ at $(refwl)A is negative, set host galaxy guess value at zero."
            getmodel(fp, ith)[:Galaxy].norm.val   = 0.
        else
            getmodel(fp, ith)[:Galaxy].norm.val   = 1/2 * vv
            getmodel(fp, ith)[:QSOcont].norm.val *= 1/2 * vv / Dierckx.Spline1D(λ, geteval(fp, ith, :QSOcont), k=1, bc="extrapolate")(refwl)
        end
    end
end


function renorm_cont!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    renorm_cont!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
end
function renorm_cont!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    c = getmodel(fp, ith)[:QSOcont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println("Cont. norm. (before): ", c.norm.val)
        scaling = 0.99
        while c.norm.val * scaling > c.norm.low
            residuals = (geteval(fp, ith) - values(getdata(fp, ith))) ./ uncerts(getdata(fp, ith))
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= scaling
        end
        println("Cont. norm. (after) : ", c.norm.val)
    else
        println("Skipping cont. renormalization")
    end
end


function guess_norm_factor!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    guess_norm_factor!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
end
function guess_norm_factor!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int, name::Symbol; quantile=0.95)
    @track_recipe
    @assert getmodel(fp, ith)[name].norm.val != 0
    m = geteval(fp, ith, name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    r = values(getdata(fp, ith)) - geteval(fp, ith)
    ratio = getmodel(fp, ith)[name].norm.val / sum(m[i1:i2])
    off = sum(r[i1:i2]) * ratio
    getmodel(fp, ith)[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if getmodel(fp, ith)[name].norm.val < 0  # ensure line has positive normalization
        getmodel(fp, ith)[name].norm.val = abs(off)
    end
end


function add_emission_lines!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    add_emission_lines!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
    for i in 2:nmodels(fp)
        m1 = getmodel(fp, 1)
        mi = getmodel(fp, i)
        if haskey(mi, :OIII_5007)  &&  haskey(m1, :OIII_5007)
            mi[:OIII_5007].norm.mpatch = @fd m -> m[1][:OIII_5007].norm
        end
        # if haskey(mi, :Ha_br)  &&  haskey(m1, :Hb_br)
        #     mi[:Ha_br].fwhm.mpatch = @fd m -> m[1][:Hb_br].fwhm
        # end
    end
end
function add_emission_lines!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    model = getmodel(fp, ith)
    lines = recipe.specs[ith].meta[:lines]

    # Create model components
    for (cname, line) in lines
        model[cname] = line.comp
        push!(model[line.group].list, cname)
    end

    # Patch functions
    add_patch_functs!(recipe, fp, ith)

    # Guess normalizations
    for group in [:BroadLines, :NarrowLines, :VeryBroadLines]  # Note: order is important
        for (cname, line) in lines
            (group == line.group)  ||  continue
            guess_norm_factor!(recipe, fp, ith, cname)
        end
    end
end


function fit_nuisance_lines!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    fit_nuisance_lines!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
end
function fit_nuisance_lines!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    (recipe.n_nuisance > 0)  ||  (return nothing)
    model = getmodel(fp, ith)

    # Prepare nuisance line components
    model[:NuisanceLines] = SumReducer([Symbol(:nuisance, i) for i in 1:recipe.n_nuisance])
    push!(model[:main].list, :NuisanceLines)
    for i in 1:recipe.n_nuisance
        model[Symbol(:nuisance, i)] = line_component(recipe, 3000., NuisanceLine)
        freeze!(model, Symbol(:nuisance, i))
    end

    # Set "nuisance" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λ = coords(getdomain(fp, ith))
    λnuisance = Vector{Float64}()
    while true
        (length(λnuisance) >= recipe.n_nuisance)  &&  break
        Δ = (values(getdata(fp, ith)) - geteval(fp, ith)) ./ uncerts(getdata(fp, ith))

        # Avoid considering again the same region (within 1A) TODO: within resolution
        for l in λnuisance
            Δ[findall(abs.(l .- λ) .< 1)] .= 0.
        end

        # Avoidance regions
        for rr in recipe.nuisance_avoid
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
        model[cname].norm.val = 1.
        model[cname].center.val  = λ[iadd]

        # Allow to shift by a quantity equal to ...
        @assert recipe.nuisance_maxoffset_from_guess > 0
        model[cname].center.low  = λ[iadd] * (1 - recipe.nuisance_maxoffset_from_guess / 3e5)
        model[cname].center.high = λ[iadd] * (1 + recipe.nuisance_maxoffset_from_guess / 3e5)

        # In any case, we must stay out of avoidance regions
        for rr in recipe.nuisance_avoid
            @assert !(rr[1] .< λ[iadd] .< rr[2])
            if rr[1] .>= λ[iadd]
                model[cname].center.high = min(model[cname].center.high, rr[1])
            end
            if rr[2] .<= λ[iadd]
                model[cname].center.low  = max(model[cname].center.low, rr[2])
            end
        end

        thaw!(model, cname)
        fit!(recipe, fp)
        freeze!(model, cname)
    end
end


function neglect_weak_features!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem)
    @track_recipe
    return neglect_weak_features!.(Ref(recipe), Ref(fp), 1:nmodels(fp))
end
function neglect_weak_features!(recipe::CRecipe{<: QSOGeneric}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    model = getmodel(fp, ith)

    # Disable "nuisance" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    if :NuisanceLines in keys(model)
        for cname in model[:NuisanceLines].list
            isfreezed(model, cname)  &&  continue
            if model[cname].norm.val == 0.
                freeze!(model, cname)
                needs_fitting = true
                println("Disabling $cname (norm. = 0)")
            elseif model[cname].norm.unc / model[cname].norm.val > 3
                model[cname].norm.val = 0.
                    freeze!(model, cname)
                needs_fitting = true
                println("Disabling $cname (unc. / norm. > 3)")
            end
        end
    end
    return needs_fitting
end


function postanalysis(recipe::CRecipe{<: QSOGeneric}, bestfit::GModelFit.ModelSnapshot)
    @track_recipe
    EW = OrderedDict{Symbol, Float64}()

    cnames = Symbol[]
    for i in 1:length(recipe.specs)
        append!(cnames, keys(recipe.specs[i].meta[:lines]))
    end
    cnames = sort(unique(cnames))

    cont = deepcopy(bestfit())
    for cname in cnames
        haskey(bestfit, cname) || continue
        cont .-= bestfit(cname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for cname in cnames
        haskey(bestfit, cname) || continue
        EW[cname] = QSFit.int_tabulated(coords(bestfit.domain),
                                        bestfit(cname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end


include("QSORecipes_Type1.jl")
include("QSORecipes_Type2.jl")

end
