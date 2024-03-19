module LineFitRecipes

using Statistics
using ..QSFit, GModelFit, Gnuplot

import QSFit: init_recipe!, preprocess_spec!, line_profile, analyze

export LineFit, InteractiveLineFit

abstract type LineFit <: AbstractRecipeSpec end
abstract type InteractiveLineFit <: LineFit end

function wait_mouse(sid=Gnuplot.options.default)
    gpexec(sid, "pause mouse")
    return gpvars(sid, "MOUSE")
end


line_profile(recipe::Recipe{<: LineFit}, ::Type{<: AbstractLine}, id::Val) = recipe.line_profiles
line_profile(recipe::Recipe{<: LineFit}, ::Type{<: AbstractLine})          = recipe.line_profiles

function init_recipe!(recipe::Recipe{T}) where T <: LineFit
    @invoke init_recipe!(recipe::Recipe{<: AbstractRecipeSpec})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.line_profiles = :gauss
    recipe.lines = LineDescriptor[]
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

    for ld in recipe.:lines
        for T in ld.types
            cname = line_cname(recipe, T, ld.id)
            @assert !(cname in keys(model)) "Duplicate component name: $cname"
            model[cname] = line_component(recipe, T, ld.id)

            group = line_group(recipe, T)
            if !(group in keys(model))
                model[group] = SumReducer(cname)
                push!(model[:main].list, group)
            else
                push!(model[group].list, cname)
            end
        end
    end

    display(model)
    display(resid.meval.model)
    
    return fit!(resid)
end


function analyze(recipe::Recipe{<: InteractiveLineFit}, spec::QSFit.Spectrum, resid::GModelFit.Residuals)
    ld = Vector{LineDescriptor}()
    linetypes = [ForbiddenLine,
                 NarrowLine,
                 BroadLine,
                 VeryBroadLine,
                 NuisanceLine]

    @gp :LineFit spec
    println()
    printstyled("Interactive emission line fit\n", color=:green)
    printstyled("Step 1: ", color=:green);  println("zoom on the desired range")
    println("(use mouse right button to start drawing a region, then left button to zoom)")
    print("Press ENTER to continue...")
    readline()
    vars = gpvars(:LineFit)
    range = [vars.X_MIN, vars.X_MAX]

    println()
    printstyled("Step 2: ", color=:green);  println("Identify emission lines by click on the peaks (CTRL+left click to stop)")
    while true
        print("Waiting for a click...")
        vars = wait_mouse(:LineFit)
        (vars.MOUSE_CTRL != 0)  &&  break
        print("\r")
        if length(ld) == 0
            for i in 1:length(linetypes)
                println("$(i): ", linetypes[i])
            end
            println("0: skip this line")
        end
        λ = round(vars.MOUSE_X * 1e2) / 1e2
        print("Insert space separated list of line type(s) for the emission line at $λ ", spec.unit_x, ": ")
        i = Int.(Meta.parse.(string.(split(readline()))))
        (0 in i)  &&  continue
        @assert all(1 .<= i .<= length(linetypes))
        push!(ld, LineDescriptor(λ, linetypes[i]...))
    end
    println()
    Gnuplot.quit(:LineFit)

    recipe.wavelength_range = range
    recipe.lines = ld

    println()
    printstyled("Step 3: ", color=:green);  println("choose line profile:")
    println("1: Gaussian")
    println("2: Lorentzian")
    println("3: Voigt")
    print("Insert line profile: ")
    i = Int(Meta.parse(readline()))
    @assert (1 <= i <= 3)
    if i == 1
        recipe.line_profiles = :gauss
    elseif i == 2
        recipe.line_profiles = :lorentz
    else
        recipe.line_profiles = :voigt
    end

    println()
    println("recipe = Recipe(LineFit)")
    println("recipe.wavelength_range = ", recipe.wavelength_range)
    println("recipe.line_profiles = :", recipe.line_profiles)
    println("recipe.lines = [")
    for i in 1:length(ld)
        println("     LineDescriptor(", ld[i].id, ", ", join(string.(ld[i].types), ", "), ")")
    end
    println("]")

    return @invoke analyze(recipe::Recipe{<: LineFit}, spec, resid)
end

end
