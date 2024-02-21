using GModelFitViewer
import GModelFitViewer: ViewerData

function ViewerData(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes)))
    keep = string.(keys(res.bestfit))[i]
    meta = GModelFitViewer.Meta(; title=res.source.name * ", z=" * string(res.source.z) * ", E(B-V)=" * string(res.source.MW_ebv),
                                xlabel="Rest frame wavelength",
                                xunit=string(QSFit.unit_λ()),
                                xscale=10. ^ QSFit.log10_scale_λ(),
                                ylabel="Lum. density",
                                yunit=string(QSFit.unit_lum_density()),
                                yscale=10. ^ QSFit.log10_scale_lum(),
                                keep=keep, kws...)

    d = GModelFitViewer.ViewerData(res.bestfit, res.fitstats, res.pspec.data, meta=meta)

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
