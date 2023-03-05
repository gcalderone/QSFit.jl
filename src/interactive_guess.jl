using REPL.TerminalMenus
export interactive_guess

function interactive_guess(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: DefaultRecipe
    cnames = keys(model)
    menu = RadioMenu(string.(cnames), pagesize=20)
    choice = request("Choose a component:", menu)
    (choice == -1)  &&  (return nothing)
    interactive_guess(source, pspec, model, cnames[choice])
end


function interactive_guess(source::QSO{T}, pspec::PreparedSpectrum, model::Model, cname::Symbol) where T <: DefaultRecipe
    pnames = collect(keys(GModelFit.getparams(model[cname])))
    choices = [string(p.name) * (p.index > 0 ? string(p.index) : "") for p in pnames]
    menu = RadioMenu(choices, pagesize=20)
    choice = request("Choose a parameter of component $cname:", menu)
    (choice == -1)  &&  (return nothing)
    interactive_guess(source, pspec, model, cname, pnames[choice])
end


interactive_guess(source::QSO{T}, pspec::PreparedSpectrum, model::Model, cname::Symbol, pname::Symbol) where T <: DefaultRecipe =
    interactive_guess(source, pspec, model, cname, GModelFit.ParamID(pname, 0))


function interactive_guess(source::QSO{T}, pspec::PreparedSpectrum, model::Model, cname::Symbol, pid::GModelFit.ParamID) where T <: DefaultRecipe
    if pid.index == 0
        par = getproperty(model[cname], pid.name)
    else
        par = getproperty(model[cname], pid.name)[pid.index]
    end
    orig_val = par.val
    println("Component: $cname, Parameter: $(pid.name)", (pid.index > 0  ?  "[" * string(pid.index) * "]"  :  ""))
    println("Initial value: ", par.val)
    println("('q' to accept current value, 'c' to restore initial value)")
    while true
        evaluate!(model)
        @gp    :interactive_guess title=(source.name * ", " * pspec.orig.label) "set grid" :-
        @gp :- :interactive_guess reverse([Gnuplot.recipe(model)..., Gnuplot.recipe(pspec.data)])
        print("Insert new value: ")
        input = readline()
        if input == "q"
            println("Consider adding the following code to your recipe:")
            println()
            println("function guess_emission_lines!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: ", string(typeof(source))[5:end-1])
            println("    guess_emission_lines!(parent_recipe(source), pspec, model)")
            println("    if (source.name == \"", source.name, "\")  &&  (pspec.orig.label == \"", pspec.orig.label, "\")")
            println("        model[:$(cname)].$(pid.name)", (pid.index > 0  ?  "[" * string(pid.index) * "]"  :  ""), ".val = ", par.val)
            println("    end")
            println("end")
            println()
            println("(press ENTER to continue)"); readline()
            break
        elseif input == "c"
            par.val = orig_val
            break
        else
            par.val = Meta.parse(input)
        end
    end
    println("Using value $(par.val)")
    evaluate!(model)
    Gnuplot.quit(:interactive_guess)
end
