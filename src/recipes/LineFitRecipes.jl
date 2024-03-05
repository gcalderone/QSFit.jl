module LineFitRecipes

using Statistics
using ..QSFit, GModelFit, Gnuplot

import QSFit: line_profile, default_options, analyze

export LineFitRecipe, InteractiveLineFitRecipe

abstract type LineFitRecipe <: AbstractRecipe end
abstract type InteractiveLineFitRecipe <: LineFitRecipe end

function wait_mouse(sid=Gnuplot.options.default)
    gpexec(sid, "pause mouse")
    return gpvars(sid, "MOUSE")
end


line_profile(recipe::RRef{<: LineFitRecipe}, ::Type{<: AbstractLine}, id::Val) = recipe.options[:line_profiles]
line_profile(recipe::RRef{<: LineFitRecipe}, ::Type{<: AbstractLine})          = recipe.options[:line_profiles]

function default_options(::Type{T}) where T <: LineFitRecipe
    out = default_options(supertype(T))
    out[:wavelength_range] = [1215, 7.3e3]
    out[:line_profiles] = :gauss
    out[:lines] = LineDescriptor[]
    return out
end


function analyze(recipe::RRef{<: LineFitRecipe}, state::QSFit.State)
    state.spec.good[findall(state.spec.x .< recipe.options[:wavelength_range][1])] .= false
    state.spec.good[findall(state.spec.x .> recipe.options[:wavelength_range][2])] .= false
    QSFit.update_data!(state)

    model = Model()
    model[:QSOcont] = QSFit.powerlaw(median(coords(domain(state.data))))
    model[:QSOcont].norm.val = median(values(state.data))
    model[:QSOcont].norm.low = median(values(state.data)) / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    model[:QSOcont].alpha.val  = -1.5
    model[:QSOcont].alpha.low  = -3
    model[:QSOcont].alpha.high =  1

    model[:main] = SumReducer(:QSOcont)
    select_maincomp!(model, :main)

    for ld in recipe.options[:lines]
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

    mzer = GModelFit.cmpfit()
    mzer.config.ftol = 1.e-6
    return fit(model, state.data, minimizer=mzer)
end


function analyze(recipe::RRef{<: InteractiveLineFitRecipe}, state::QSFit.State)
    ld = Vector{LineDescriptor}()
    linetypes = [ForbiddenLine,
                 NarrowLine,
                 BroadLine,
                 VeryBroadLine,
                 NuisanceLine]

    @gp :LineFit state.spec
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
        print("Insert space separated list of line type(s) for the emission line at $λ ", state.spec.unit_x, ": ")
        i = Int.(Meta.parse.(string.(split(readline()))))
        (0 in i)  &&  continue
        @assert all(1 .<= i .<= length(linetypes))
        push!(ld, LineDescriptor(λ, linetypes[i]...))
    end
    println()
    Gnuplot.quit(:LineFit)

    empty!(recipe.options)
    recipe.options[:wavelength_range] = range
    recipe.options[:lines] = ld

    println()
    printstyled("Step 3: ", color=:green);  println("choose line profile:")
    println("1: Gaussian")
    println("2: Lorentzian")
    println("3: Voigt")
    print("Insert line profile: ")
    i = Int(Meta.parse(readline()))
    @assert (1 <= i <= 3)
    if i == 1
        recipe.options[:line_profiles] = :gauss
    elseif i == 2
        recipe.options[:line_profiles] = :lorentz
    else
        recipe.options[:line_profiles] = :voigt
    end

    println()
    println("recipe = RRef(LineFitRecipe)")
    println("recipe.options[:wavelength_range] = ", recipe.options[:wavelength_range])
    println("recipe.options[:line_profiles] = :", recipe.options[:line_profiles])
    println("recipe.options[:lines] = [")
    for i in 1:length(ld)
        println("     LineDescriptor(", ld[i].id, ", ", join(string.(ld[i].types), ", "), ")")
    end
    println("]")

    return @invoke analyze(recipe::RRef{<: LineFitRecipe}, state)
end

end
