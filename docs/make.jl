using Revise, Documenter, GModelFit, Gnuplot, QSFit

makedocs(sitename="QSFit.jl",
         authors = "Giorgio Calderone",
         format = Documenter.HTML(prettyurls = false),  # uncomment for local use, comment for deployment
         modules=[QSFit],
         pages = [
             "Home" => "index.md",
         ])
Gnuplot.quitall()

if !(@isdefined is_compiled)
    is_compiled = true
    error("Re-run with compiled code!")
end
