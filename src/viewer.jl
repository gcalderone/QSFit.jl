using GModelFitViewer
import GModelFitViewer: ViewerData

function ViewerData_meta(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    return GModelFitViewer.Meta(; title=res.spec[:label],
                                xlabel="Wavelength",
                                xunit=res.spec[:xunit],
                                xscale=res.spec[:xscale],
                                ylabel="Lum. density",
                                yunit=res.spec[:yunit],
                                yscale=res.spec[:yscale],
                                kws...)
end


function ViewerData_extra(res::Results)
    out = OrderedDict{Symbol, Any}()
    for (tab, dict) in res.post
        @assert isa(dict, AbstractDict)

        if tab == :Data_stats
            out[tab] = OrderedDict(:label => "Data stats", :fields => OrderedDict{Symbol, Any}())
            out[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Quantity"),
                                                 :data => collect(keys(dict)))
            out[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Value"),
                                                 :data => collect(values(dict)))
        elseif tab == :Issues
            out[tab] = OrderedDict(:label => "Issues", :fields => OrderedDict{Symbol, Any}())
            out[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Component"),
                                                 :data => Vector{Symbol}())
            out[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Parameter / dependency"),
                                                 :data => Vector{Symbol}())
            out[tab][:fields][:c3] = OrderedDict(:meta => Dict(:desc => "Issue"),
                                                 :data => Vector{String}())
            for (cname, subdict) in dict
                for (par_or_dep, issue) in subdict
                    push!(out[tab][:fields][:c1][:data], cname)
                    push!(out[tab][:fields][:c2][:data], par_or_dep)
                    push!(out[tab][:fields][:c3][:data], issue)
                end
            end
        elseif tab == :Continuum_luminosity
            out[tab] = OrderedDict(:label => "Continuum luminosities", :fields => OrderedDict{Symbol, Any}())
            out[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Wavelength [A]"),
                                                 :data => collect(keys(dict)))
            out[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "lambda L_lambda [10^42 erg/s]"),
                                                 :data => collect(values(dict)))
        elseif tab == :Equivalent_widths
            out[tab] = OrderedDict(:label => "Equivalent widths", :fields => OrderedDict{Symbol, Any}())
            out[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Emission line"),
                                                 :data => collect(keys(dict)))
            out[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "[A]"),
                                                 :data => collect(values(dict)))
        elseif isa(dict, OrderedDict)
            vv = string.(collect(values(dict)))
            out[tab] = OrderedDict(:label => tab, :fields => OrderedDict{Symbol, Any}())
            out[tab][:fields][:c1] = OrderedDict(:meta => Dict(:desc => "Key"),
                                                 :data => collect(keys(dict)))
            out[tab][:fields][:c2] = OrderedDict(:meta => Dict(:desc => "Value"),
                                                 :data => collect(values(dict)))
        else
            @warn "Unexpected entry in post-analysis dictionary: $tab"
        end
    end
    return [out]
end


function ViewerData(res::Results; kws...)
    meta = ViewerData_meta(res; kws...)
    out = GModelFitViewer.ViewerData(res.bestfit, res.fsumm, res.data, meta=meta)
    out.data[1]["extra"] = ViewerData_extra(res)
    return out
end


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
