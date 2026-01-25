using GModelFitViewer
import GModelFitViewer: ViewerOpts, _priv_lower

function ViewerOpts(res::Results; kws...)
    # ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    return GModelFitViewer.ViewerOpts(; title=res.spec[:label],
                                      xlabel="Wavelength",
                                      xunit=res.spec[:xunit],
                                      xscale=res.spec[:xscale],
                                      ylabel="Lum. density",
                                      yunit=res.spec[:yunit],
                                      yscale=res.spec[:yscale],
                                      kws...)
end


function _priv_lower(res::Results; kws...)
    extra = Vector{AuxTable}()
    for (tab, dict) in res.post
        @assert isa(dict, AbstractDict)

        if tab == :Data_stats
            e = AuxTable("Data stats")
            add_col!(e, "Quantity", collect(keys(dict)))
            add_col!(e, "Value",    collect(values(dict)))
            push!(extra, e)
        elseif tab == :Issues
            e = AuxTable("Issues")
            add_col!(e, "Component", String[])
            add_col!(e, "Parameter / dependency", String[])
            add_col!(e, "Issue",    String[])

            for (cname, subdict) in dict
                for (par_or_dep, issue) in subdict
                    push!(e.fields[1].data, string(cname))
                    push!(e.fields[2].data, string(par_or_dep))
                    push!(e.fields[3].data, issue)
                end
            end
            push!(extra, e)
        elseif tab == :Continuum_luminosity
            e = AuxTable("Continuum luminosities")
            add_col!(e, "Wavelength [A]", collect(keys(dict)))
            add_col!(e, "lambda L_lambda [10^42 erg/s]", collect(values(dict)))
            push!(extra, e)
        elseif tab == :Equivalent_widths
            e = AuxTable("Equivalent widths")
            add_col!(e, "Emission line", collect(keys(dict)))
            add_col!(e, "[A]", collect(values(dict)))
            push!(extra, e)
        elseif isa(dict, OrderedDict)
            e = AuxTable(string(tab))
            add_col!(e, "Key", collect(keys(dict)))
            add_col!(e, "Value", collect(values(dict)))
            push!(extra, e)
        else
            @warn "Unexpected entry in post-analysis dictionary: $tab"
        end
    end
    out = GModelFitViewer._priv_lower(res.bestfit, res.fsumm, res.data, ViewerOpts(res; kws...), extra)
    return out
end


#=
function ViewerData(multires::MultiResults; kws...)
    N = length(multires.bestfit)
    res = [Results(multires.timestamp,
                   multires.elapsed,
                   multires.spec[i],
                   multires.data[i],
                   multires.bestfit[i],
                   multires.fsumm::GModelFit.FitSummary,
                   multires.post[i]) for i in 1:N]

    meta = [ViewerData_meta(res[i]; kws...) for i in 1:N]
    out = GModelFitViewer.ViewerData(multires.bestfit, multires.fsumm, multires.data, meta=meta)
    for i in 1:N
        out.data[1][i]["extra"] = ViewerData_extra(res[i])
    end
    return out
end
=#
