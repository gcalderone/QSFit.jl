using GModelFitViewer
import GModelFitViewer: ViewerOpts, _priv_lower

function ViewerOpts(spec::AbstractDict; kws...)
    # ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    return GModelFitViewer.ViewerOpts(; title=spec[:label],
                                      xlabel="Wavelength",
                                      xunit=spec[:xunit],
                                      xscale=spec[:xscale],
                                      ylabel="Lum. density",
                                      yunit=spec[:yunit],
                                      yscale=spec[:yscale],
                                      kws...)
end


function format_post(post::AbstractDict)
    extra = Vector{AuxTable}()
    for (tab, dict) in post
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
    return extra
end

function _priv_lower(res::Results; kws...)
    extra = format_post(res.post)
    out = GModelFitViewer._priv_lower(res.bestfit, res.fsumm, res.data, ViewerOpts(res.spec; kws...), extra)
    return out
end


function _priv_lower(res::MultiResults; kws...)
    extra = [format_post(res.post[i]) for i in 1:length(res.post)]
    opts  = [ViewerOpts(res.spec[i]) for i in 1:length(res.spec)]
    out = GModelFitViewer._priv_lower(res.bestfit, res.fsumm, res.data, opts, extra)
end
