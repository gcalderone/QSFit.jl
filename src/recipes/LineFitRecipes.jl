module LineFitRecipes

using Statistics
using ..QSFit, GModelFit, Gnuplot

import QSFit: Food, init_recipe!, preprocess_spec!, line_component, analyze

export LineFit, InteractiveLineFit

abstract type LineFit <: AbstractRecipe end
abstract type InteractiveLineFit <: LineFit end

function wait_mouse(sid=Gnuplot.options.default)
    gpexec(sid, "pause mouse")
    return gpvars(sid, "MOUSE")
end

line_component(recipe::CRecipe{<: LineFit}, center::Float64) = recipe.line_component(center)

function init_recipe!(recipe::CRecipe{T}) where T <: LineFit
    @invoke init_recipe!(recipe::CRecipe{<: AbstractRecipe})
    recipe.wavelength_range = [1215, 7.3e3]
    recipe.line_component = QSFit.SpecLineGauss
end


function preprocess_spec!(recipe::CRecipe{T}, spec::QSFit.Spectrum) where T <: LineFit
    @invoke preprocess_spec!(recipe::CRecipe{<: AbstractRecipe}, spec)
    spec.good[findall(spec.x .< recipe.wavelength_range[1])] .= false
    spec.good[findall(spec.x .> recipe.wavelength_range[2])] .= false
end


function analyze(recipe::CRecipe{<: LineFit}, food::Food)
    food.model[:QSOcont] = QSFit.powerlaw(median(coords(domain(food.data))))
    food.model[:QSOcont].norm.val = median(values(food.data))
    food.model[:QSOcont].norm.low = median(values(food.data)) / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    food.model[:QSOcont].alpha.val  = -1.5
    food.model[:QSOcont].alpha.low  = -3
    food.model[:QSOcont].alpha.high =  1

    food.model[:main] = SumReducer(:QSOcont)
    select_maincomp!(food.model, :main)

    for (cname, line) in recipe.lines
        food.model[cname] = line.comp
        if !(line.group in keys(food.model))
            food.model[line.group] = SumReducer(cname)
            push!(food.model[:main].list, line.group)
        else
            push!(food.model[line.group].list, cname)
        end
    end

    solver = GModelFit.cmpfit()
    solver.config.ftol = 1.e-6
    return GModelFit.fit!(food.model, food.data, solver)
end


function analyze(recipe::CRecipe{<: InteractiveLineFit}, food::Food)
    @gp :LineFit food.spec

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
            println()
            println("List of line templates:")
            for i in 1:length(linetemplates)
                println("$(i): ", linetemplates[i])
            end
            println("0: skip this line")
        end
        print("Insert space separated list of line template(s) for the emission line at $λ ", food.spec.unit_x, ": ")
        i = Int.(Meta.parse.(string.(split(readline()))))
        (0 in i)  &&  continue
        @assert all(1 .<= i .<= length(linetemplates))
        add_line!(recipe, λ, linetemplates[i]...)
        push!(accum_lines, "add_line!(recipe, $λ, " * join(string.(linetemplates[i]), ",") * ")")
    end
    println()
    Gnuplot.quit(:LineFit)

    println()
    printstyled("Reproduce analysis with the LineFit recipe:\n", color=:green)
    println("recipe = CRecipe(LineFit, ... keywords ...)")
    println("... set recipe properties ...")
    println("recipe.wavelength_range = ", recipe.wavelength_range)
    for l in accum_lines
        println(l)
    end
    println()

    recipe2 = CRecipe{LineFit}()
    for p in propertynames(recipe)
        setproperty!(recipe2, p, getproperty(recipe, p))
    end
    res = analyze(recipe2, food.origspec)
    return res.bestfit, res.fsumm    
end

end
