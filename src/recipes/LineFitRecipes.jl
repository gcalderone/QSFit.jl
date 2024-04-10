module LineFitRecipes

using Statistics
using ..QSFit, ..QSFit.ATL, GModelFit, Gnuplot

import QSFit: init_recipe!, preprocess_spec!, line_component, analyze

export LineFit, InteractiveLineFit

abstract type LineFit <: AbstractRecipeSpec end
abstract type InteractiveLineFit <: LineFit end

function wait_mouse(sid=Gnuplot.options.default)
    gpexec(sid, "pause mouse")
    return gpvars(sid, "MOUSE")
end


function line_component(recipe::Recipe{<: LineFit}, ::Type{<: QSFit.LineTemplate}, center::Float64)
    if recipe.line_profiles == :gauss
        return QSFit.SpecLineGauss(center)
    elseif recipe.line_profiles == :lorentz
        return QSFit.SpecLineLorentz(center)
    elseif recipe.line_profiles == :voigt
        return QSFit.SpecLineVoigt(center)
    end
end


function init_recipe!(recipe::Recipe{T}) where T <: LineFit
    @invoke init_recipe!(recipe::Recipe{<: AbstractRecipeSpec})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.line_profiles = :gauss
end


function preprocess_spec!(recipe::Recipe{T}, spec::QSFit.Spectrum) where T <: LineFit
    @invoke preprocess_spec!(recipe::Recipe{<: AbstractRecipeSpec}, spec)
    spec.good[findall(spec.x .< recipe.wavelength_range[1])] .= false
    spec.good[findall(spec.x .> recipe.wavelength_range[2])] .= false
end


function analyze(recipe::Recipe{<: LineFit}, spec::QSFit.Spectrum, resid::GModelFit.Residuals)
    resid.mzer.config.ftol = 1.e-6
    model = resid.meval.model
    model[:QSOcont] = QSFit.powerlaw(median(coords(domain(resid.data))))
    model[:QSOcont].norm.val = median(values(resid.data))
    model[:QSOcont].norm.low = median(values(resid.data)) / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    model[:QSOcont].alpha.val  = -1.5
    model[:QSOcont].alpha.low  = -3
    model[:QSOcont].alpha.high =  1

    model[:main] = SumReducer(:QSOcont)
    select_maincomp!(model, :main)

    for (cname, line) in recipe.lines
        model[cname] = line.comp
        if !(line.group in keys(model))
            model[line.group] = SumReducer(cname)
            push!(model[:main].list, line.group)
        else
            push!(model[line.group].list, cname)
        end
    end

    return GModelFit.minimize!(resid)
end


function analyze(_recipe::Recipe{<: InteractiveLineFit}, _spec::Spectrum) where T <: AbstractRecipeSpec
    recipe = deepcopy(_recipe)
    spec = deepcopy(_spec)
    preprocess_spec!(recipe, spec)
    @gp :LineFit spec

    println()
    printstyled("Interactive emission line fit\n", color=:green)
    printstyled("Step 1: ", color=:green);  println("zoom on the desired range")
    println("(use mouse right button to start drawing a region, then left button to zoom)")
    print("Press ENTER to continue...")
    readline()
    vars = gpvars(:LineFit)
    range = [vars.X_MIN, vars.X_MAX]
    recipe.wavelength_range = range

    linetemplates = [ForbiddenLine,
                     NarrowLine,
                     BroadLine,
                     NuisanceLine]
    println()
    printstyled("Step 2: ", color=:green);  println("Identify emission lines by click on the peaks (CTRL+left click to stop)")
    accum_lines = String[]
    while true
        print("Waiting for a click...")
        vars = wait_mouse(:LineFit)
        (vars.MOUSE_CTRL != 0)  &&  break
        print("\r")
        λ = round(vars.MOUSE_X * 1e2) / 1e2
        if length(get_lines_dict(recipe)) == 0
            for i in 1:length(linetemplates)
                println("$(i): ", linetemplates[i])
            end
            println("0: skip this line")
        end
        print("Insert space separated list of line type(s) for the emission line at $λ ", spec.unit_x, ": ")
        i = Int.(Meta.parse.(string.(split(readline()))))
        (0 in i)  &&  continue
        @assert all(1 .<= i .<= length(linetemplates))
        add_line!(recipe, ATL.UnidentifiedTransition(λ), linetemplates[i]...)
        push!(accum_lines, "add_line!(recipe, QSFit.ATL.UnidentifiedTransition($(λ)), " * join(string.(linetemplates[i]), ",") * ")")
    end
    println()
    Gnuplot.quit(:LineFit)

    println()
    printstyled("Reproduce analysis with:\n", color=:green)
    println("recipe.wavelength_range = ", recipe.wavelength_range)
    println("recipe.line_profiles = :", recipe.line_profiles)
    for l in accum_lines
        println(l)
    end
    println()
    return @invoke analyze(recipe::Recipe{<: LineFit}, _spec)
end

end
