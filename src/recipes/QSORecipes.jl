module QSORecipes

using Printf, DataStructures, Statistics, Dates
using Dierckx
using ..QSFit, ..QSFit.ATL, GModelFit

import GModelFit: domain
import QSFit: init_recipe!, preprocess_spec!, spec2data, line_suffix, line_component, set_lines_dict!, analyze, reduce, Results

abstract type QSOGeneric <: AbstractRecipe end

#=
Note: the following structure has the same name and purpose as the one in GModelFit, the only difference being that the followinf can trigger a multi-model scan and evaluation
=#
struct DataModel{N}
    parent::GModelFit.MEval{N}
    data::Vector{Measures{1}}
    id::Int
    DataModel(meval::GModelFit.MEval{N}, data::Vector{Measures{1}}, id::Int) where N = new{N}(meval, data, id)
end

getmodel(s::DataModel) = s.parent.v[s.id].model
getdomain(s::DataModel) = s.parent.v[s.id].domain
getdata(s::DataModel) = s.data[s.id]

import GModelFit: last_eval
last_eval(s::DataModel) = GModelFit.last_eval(s.parent, s.id)
last_eval(s::DataModel, cname::Symbol) = GModelFit.last_eval(s.parent, s.id, cname)


scan_and_evaluate!(dm::Vector{DataModel{N}}) where N = scan_and_evaluate!(dm[1])
function scan_and_evaluate!(meval::DataModel)
    GModelFit.scan_model!(meval.parent)
    GModelFit.update_eval!(meval.parent)
end

line_component(recipe::CRecipe{<: QSOGeneric}, center::Float64) = recipe.line_component(center)

abstract type BlueWing <: NarrowLine end
line_suffix(recipe::CRecipe{<: QSOGeneric}, ::Type{BlueWing}) = :_bw

function line_component(recipe::CRecipe, center::Float64, ::Type{<: BlueWing})
    @track_recipe
    comp = line_component(recipe, center, NarrowLine)
    comp.fwhm.low, comp.fwhm.val, comp.fwhm.high = 0, 3e3, 5e3
    comp.voff.low, comp.voff.val, comp.voff.high = 0, 0, 2e3
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


fit!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N = fit!(recipe, dm[1])
function fit!(recipe::CRecipe{<: QSOGeneric}, dm::DataModel{N}) where N
    @track_recipe
    scan_and_evaluate!(dm)
    fp = GModelFit.FitProblem(dm.parent, dm.data)
    bestfit, fsumm = GModelFit.fit!(fp, recipe.solver)
    show(fsumm)
    return bestfit[1], fsumm
end



function add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_qso_continuum!.(Ref(recipe), dm)
    scan_and_evaluate!(dm)
end
function add_qso_continuum!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe
    λ = coords(domain(getdata(meval)))
    comp = QSFit.powerlaw(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(getdata(meval))) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5
    comp.alpha.low  = -3
    comp.alpha.high =  1

    getmodel(meval)[:QSOcont] = comp
    push!(getmodel(meval)[:Continuum].list, :QSOcont)
end


function add_host_galaxy!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_host_galaxy!.(Ref(recipe), dm)
    for i in 2:length(dm)
        getmodel(dm[i])[:Galaxy].norm.mpatch = @fd m -> m[1][:Galaxy].norm
    end
    scan_and_evaluate!(dm)
end
function add_host_galaxy!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe
    λ = coords(getdomain(meval))
    if recipe.use_host_template  &&
        (recipe.host_template_range[1] .< maximum(λ))  &&
        (recipe.host_template_range[2] .> minimum(λ))
        getmodel(meval)[:Galaxy] = QSFit.hostgalaxy(recipe.host_template[:template],
                                                library=recipe.host_template[:library],
                                                refwl=recipe.host_template[:ref_wavelength])
        push!(getmodel(meval)[:Continuum].list, :Galaxy)

        # Split total flux between continuum and host galaxy
        refwl = recipe.host_template[:ref_wavelength]
        vv = Dierckx.Spline1D(λ, values(getdata(meval)), k=1, bc="extrapolate")(refwl)
        @assert !isnan(vv) "Predicted L_λ at $(refwl)A is NaN"
        if vv <= 0
            @warn "Predicted L_λ at $(refwl)A is negative, set host galaxy guess value at zero."
            getmodel(meval)[:Galaxy].norm.val   = 0.
        else
            getmodel(meval)[:Galaxy].norm.val   = 1/2 * vv
            getmodel(meval)[:QSOcont].norm.val *= 1/2 * vv / Dierckx.Spline1D(λ, GModelFit.last_eval(meval, :QSOcont), k=1, bc="extrapolate")(refwl)
        end
    end
end


function renorm_cont!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    renorm_cont!.(Ref(recipe), dm)
end
function renorm_cont!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe
    freeze!(getmodel(meval), :QSOcont)
    c = getmodel(meval)[:QSOcont]
    initialnorm = c.norm.val
    if c.norm.val > 0
        println("Cont. norm. (before): ", c.norm.val)
        scaling = 0.99
        while c.norm.val * scaling > c.norm.low
            residuals = (GModelFit.last_eval(meval) - values(getdata(meval))) ./ uncerts(getdata(meval))
            ratio = count(residuals .< 0) / length(residuals)
            (ratio > 0.9)  &&  break
            (c.norm.val < (initialnorm / 5))  &&  break # give up
            c.norm.val *= scaling
            scan_and_evaluate!(meval)
        end
        println("Cont. norm. (after) : ", c.norm.val)
    else
        println("Skipping cont. renormalization")
    end
    scan_and_evaluate!(meval)
end


function guess_norm_factor!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    guess_norm_factor!.(Ref(recipe), dm)
end
function guess_norm_factor!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel, name::Symbol; quantile=0.95)
    @track_recipe
    @assert getmodel(meval)[name].norm.val != 0
    m = GModelFit.last_eval(meval, name)
    c = cumsum(m)
    @assert maximum(c) != 0. "Model for $name evaluates to zero over the whole domain"
    c ./= maximum(c)
    i1 = findfirst(c .> ((1 - quantile)/2))
    i2 = findlast( c .< ((1 + quantile)/2))
    if i1 >= i2
        return #Can't calculate normalization for component
    end
    r = values(getdata(meval)) - GModelFit.last_eval(meval)
    ratio = getmodel(meval)[name].norm.val / sum(m[i1:i2])
    off = sum(r[i1:i2]) * ratio
    getmodel(meval)[name].norm.val += off
    @assert !isnan(off) "Norm. offset is NaN for $name"
    if getmodel(meval)[name].norm.val < 0  # ensure line has positive normalization
        getmodel(meval)[name].norm.val = abs(off)
    end
end


function add_emission_lines!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_emission_lines!.(Ref(recipe), dm)
end
function add_emission_lines!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe
    groups = OrderedDict{Symbol, Vector{Symbol}}()

    # Create model components
    for (cname, line) in recipe.lines
        getmodel(meval)[cname] = line.comp
        haskey(groups, line.group)  ||  (groups[line.group] = Vector{Symbol}())
        push!(groups[line.group], cname)
    end

    # Create reducers for groups
    for (group, lnames) in groups
        @assert group in [:BroadLines, :NarrowLines, :VeryBroadLines] "Unexpected group for emission lines: $group"
        getmodel(meval)[group] = SumReducer(lnames)
        push!(getmodel(meval)[:main].list, group)
    end
    scan_and_evaluate!(meval)

    # Guess normalizations
    for group in [:BroadLines, :NarrowLines, :VeryBroadLines]  # Note: order is important
        if group in keys(groups)
            for cname in groups[group]
                guess_norm_factor!(recipe, meval, cname)
            end
        end
    end
end


function add_nuisance_lines!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    add_nuisance_lines!.(Ref(recipe), dm)
end
function add_nuisance_lines!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe
    (recipe.n_nuisance > 0)  ||  (return nothing)

    # Prepare nuisance line components
    for i in 1:recipe.n_nuisance
        getmodel(meval)[Symbol(:nuisance, i)] = line_component(recipe, 3000., NuisanceLine)
    end
    getmodel(meval)[QSFit.line_group(recipe, NuisanceLine)] = SumReducer([Symbol(:nuisance, i) for i in 1:recipe.n_nuisance])
    push!(getmodel(meval)[:main].list, QSFit.line_group(recipe, NuisanceLine))
    scan_and_evaluate!(meval)
    for j in 1:recipe.n_nuisance
        freeze!(getmodel(meval), Symbol(:nuisance, j))
    end
    scan_and_evaluate!(meval)

    # Set "nuisance" line center wavelength where there is a maximum in
    # the fit residuals, and re-run a fit.
    λ = coords(getdomain(meval))
    λnuisance = Vector{Float64}()
    while true
        (length(λnuisance) >= recipe.n_nuisance)  &&  break
        scan_and_evaluate!(meval)
        Δ = (values(getdata(meval)) - GModelFit.last_eval(meval)) ./ uncerts(getdata(meval))

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
        getmodel(meval)[cname].norm.val = 1.
        getmodel(meval)[cname].center.val  = λ[iadd]

        # Allow to shift by a quantity equal to ...
        @assert recipe.nuisance_maxoffset_from_guess > 0
        getmodel(meval)[cname].center.low  = λ[iadd] * (1 - recipe.nuisance_maxoffset_from_guess / 3e5)
        getmodel(meval)[cname].center.high = λ[iadd] * (1 + recipe.nuisance_maxoffset_from_guess / 3e5)

        # In any case, we must stay out of avoidance regions
        for rr in recipe.nuisance_avoid
            @assert !(rr[1] .< λ[iadd] .< rr[2])
            if rr[1] .>= λ[iadd]
                getmodel(meval)[cname].center.high = min(getmodel(meval)[cname].center.high, rr[1])
            end
            if rr[2] .<= λ[iadd]
                getmodel(meval)[cname].center.low  = max(getmodel(meval)[cname].center.low, rr[2])
            end
        end

        thaw!(getmodel(meval), cname);     scan_and_evaluate!(meval)
        fit!(recipe, meval)
        freeze!(getmodel(meval), cname);   scan_and_evaluate!(meval)
    end
    scan_and_evaluate!(meval)
end


function neglect_weak_features!(recipe::CRecipe{<: QSOGeneric}, dm::Vector{DataModel{N}}) where N
    out = neglect_weak_features!.(Ref(recipe), dm)
    return out
end
function neglect_weak_features!(recipe::CRecipe{<: QSOGeneric}, meval::DataModel)
    @track_recipe

    # Disable "nuisance" lines whose normalization uncertainty is larger
    # than X times the normalization
    needs_fitting = false
    for ii in 1:recipe.n_nuisance
        cname = Symbol(:nuisance, ii)
        isfreezed(getmodel(meval), cname)  &&  continue
        if getmodel(meval)[cname].norm.val == 0.
            freeze!(getmodel(meval), cname)
            needs_fitting = true
            println("Disabling $cname (norm. = 0)")
        elseif getmodel(meval)[cname].norm.unc / getmodel(meval)[cname].norm.val > 3
            getmodel(meval)[cname].norm.val = 0.
            freeze!(getmodel(meval), cname)
            needs_fitting = true
            println("Disabling $cname (unc. / norm. > 3)")
        end
    end
    return needs_fitting
end


function preprocess_spec!(recipe::CRecipe{T}, spec::QSFit.Spectrum) where T <: QSOGeneric
    @track_recipe
    @invoke preprocess_spec!(recipe::CRecipe{<: AbstractRecipe}, spec)
    recipe.spec = spec  # save spectrum for later use

    set_lines_dict!(recipe)

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
    for loop in 1:2
        # The second pass is required to neglect lines whose coverage
        # has been affected by the neglected spectral samples.
        if loop == 2
            println()
            println("Updated coverage:")
        end
        for (cname, line) in get_lines_dict(recipe)
            threshold = get(recipe.min_spectral_coverage, cname, recipe.min_spectral_coverage[:default])
            (λmin, λmax, coverage) = QSFit.spectral_coverage(spec.x[findall(spec.good)],
                                                             spec.resolution, line.comp)
            @printf("Line %-15s coverage: %5.3f on range %10.5g < λ < %10.5g", cname, coverage, λmin, λmax)
            if coverage < threshold
            @printf(", threshold is < %5.3f, neglecting...", threshold)
                spec.good[λmin .<= spec.x .< λmax] .= false
                delete!(recipe.lines, cname)
            end
            println()
        end
    end
    println("Good samples after line coverage filter: ", count(spec.good) , " / ", length(spec.good))

    # Sort lines according to center wavelength, and save list in recipe
    kk = collect(keys(  recipe.lines))
    vv = collect(values(recipe.lines))
    ii = sortperm(getfield.(getfield.(getfield.(vv, :comp), :center), :val))
    recipe.lines = OrderedDict(Pair.(kk[ii], vv[ii]))
end


function reduce(recipe::CRecipe{<: QSOGeneric}, bestfit::GModelFit.ModelSnapshot)
    @track_recipe
    EW = OrderedDict{Symbol, Float64}()

    cont = deepcopy(bestfit())
    for cname in keys(recipe.lines)
        haskey(bestfit, cname) || continue
        cont .-= bestfit(cname)
    end
    @assert all(cont .> 0) "Continuum model is zero or negative"
    for cname in keys(recipe.lines)
        haskey(bestfit, cname) || continue
        EW[cname] = QSFit.int_tabulated(coords(bestfit.domain),
                                        bestfit(cname) ./ cont)[1]
    end
    return OrderedDict{Symbol, Any}(:EW => EW)
end


include("QSORecipes_Type1.jl")
include("QSORecipes_Type2.jl")

end
