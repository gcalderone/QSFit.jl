using GModelFitViewer
import GModelFitViewer: viewer

function viewer(res::Results; kws...)
    ctypes = [comptype(res.bestfit, cname) for cname in keys(res.bestfit)]
    i = findall(isnothing.(match.(r"SpecLine", ctypes)))
    keep = string.(keys(res.bestfit))[i]
    meta = GModelFitViewer.Meta(; title=res.source.name * ", z=" * string(res.source.z) * ", E(B-V)=" * string(res.source.mw_ebv),
                                xlabel="Rest frame wavelength",
                                xunit=string(QSFit.unit_λ()),
                                xscale=10. ^ QSFit.log10_scale_λ(),
                                ylabel="Lum. density",
                                yunit=string(QSFit.unit_lum_density()),
                                yscale=10. ^ QSFit.log10_scale_lum(),
                                keep=keep, kws...)
    return viewer(res.bestfit, res.fitstats, res.pspec.data, meta=meta)
end
