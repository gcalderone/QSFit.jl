using Test, QSFit, QSFit.QSORecipes, QSFit.LineFitRecipes, GModelFitViewer

spec = Spectrum(Val(:SDSS_DR10), joinpath(QSFit.qsfit_data(), "test", "spec-0752-52251-0323.fits"))
recipe = CRecipe{Type1}(redshift=0.3806, Av=3.1 * 0.06846)
res = analyze(recipe, spec)
display(res.bestfit)
display(res.fsumm)
rm(GModelFitViewer.serialize_html(res))
@test res.fsumm.ndata == 3309
@test res.fsumm.nfree == 70


recipe = CRecipe{LineFit}(redshift=0.3806)
recipe.wavelength_range = [4530.90809628009, 5392.50547045952]
recipe.lines = QSFit.SpecLineSet()
QSFit.add_line!(recipe, recipe.lines, 4864.77, NarrowLine, BroadLine)
QSFit.add_line!(recipe, recipe.lines, 5010.88, ForbiddenLine)
res = analyze(recipe, spec)


res = analyze(recipe, [spec, deepcopy(spec)])
display(res.bestfit)
display(res.fsumm)
rm(GModelFitViewer.serialize_html(res))
@test res.fsumm.ndata == 6618
@test res.fsumm.nfree == 156


spec = Spectrum(Val(:SDSS_DR10), joinpath(QSFit.qsfit_data(), "test", "spec-2233-53845-0594.fits"))
recipe = CRecipe{Type1}(redshift=0.0999, Av=3.1 * 0.0209587)
res = analyze(recipe, spec)
display(res.bestfit)
display(res.fsumm)
rm(GModelFitViewer.serialize_html(res))
@test res.fsumm.ndata == 3118
@test res.fsumm.nfree == 83
