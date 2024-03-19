using GModelFitViewer
import GModelFitViewer: ViewerData

function ViewerData(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes)))
    keep = string.(keys(res.bestfit))[i]
    meta = GModelFitViewer.Meta(; title=res.spec.label,
                                xlabel="Wavelength",
                                xunit=string(unit(res.spec.unit_x)),
                                xscale=ustrip(res.spec.unit_x),
                                ylabel="Lum. density",
                                yunit=string(unit(res.spec.unit_y)),
                                yscale=ustrip(res.spec.unit_y),
                                keep=keep, kws...)

    d = GModelFitViewer.ViewerData(res.bestfit, res.fitstats, res.data, meta=meta)

    # Add further content
    if haskey(res.reduced, :EW)
        m = OrderedDict{Symbol, Any}()
        m[:EW] = OrderedDict{Symbol, Any}()
        m[:EW][:label] = "Em. lines EW"
        m[:EW][:fields] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label][:meta] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Label][:meta][:desc] = "Em. line"
        m[:EW][:fields][:Label][:data] = collect(keys(res.reduced[:EW]))
        m[:EW][:fields][:Value] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Value][:meta] = OrderedDict{Symbol, Any}()
        m[:EW][:fields][:Value][:meta][:desc] = "EW [A]"
        m[:EW][:fields][:Value][:data] = collect(values(res.reduced[:EW]))
        d.data[1]["extra"] = [m]
    end

    return d
end
