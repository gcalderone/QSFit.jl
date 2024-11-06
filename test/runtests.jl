using Test, QSFit, QSFit.QSORecipes, GModelFitViewer

spec = Spectrum(Val(:SDSS_DR10), joinpath(QSFit.qsfit_data(), "test", "spec-0752-52251-0323.fits"))
recipe = CRecipe{Type1}(redshift=0.3806, Av=3.1 * 0.06846)
res = analyze(recipe, spec)
display(res.bestfit)
display(res.fitstats)
rm(GModelFitViewer.serialize_html(res))
@test res.fitstats.ndata == 3309
@test res.fitstats.nfree == 70
# @test abs(res.fitstats.fitstat - 1.201116391196934) < 1e-5


spec = Spectrum(Val(:SDSS_DR10), joinpath(QSFit.qsfit_data(), "test", "spec-2233-53845-0594.fits"))
recipe = CRecipe{Type1}(redshift=0.0999, Av=3.1 * 0.0209587)
res = analyze(recipe, spec)
display(res.bestfit)
display(res.fitstats)
rm(GModelFitViewer.serialize_html(res))
@test res.fitstats.ndata == 3118
@test res.fitstats.nfree == 83
# @test abs(res.fitstats.fitstat - 1.2621845693945035) < 1e-5
