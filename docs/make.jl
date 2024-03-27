using Revise, Documenter, GModelFit, Gnuplot, QSFit

makedocs(sitename="QSFit.jl",
         authors = "Giorgio Calderone",
         format = Documenter.HTML(prettyurls = false),  # uncomment for local use, comment for deployment
         modules=[QSFit],
         pages = [
             "QSFit" => "index.md",
             "Components" => "components.md",
             "Recipes" => "recipes.md",
         ])
Gnuplot.quitall()

if !(@isdefined is_compiled)
    is_compiled = true
    error("Re-run with compiled code!")
end
