#qsfit = QSFit("dd", 0.3806, 0.06846);
#adddata!(qsfit, read_sdss_dr10("spec-0752-52251-0323.fits"));

#using Revise
include("QSFIT.jl")
#using Main.QSFIT

qsfit = QSFit("spec-1959-53440-0066.fits", 0.178064, 0.02);
adddata!(qsfit, read_sdss_dr10("spec-1959-53440-0066.fits"));
(model, bestfit) = run(qsfit);
@info "aaaa"

