
# Components

The **QSFit.jl** package implements the following GModelFit-compatible [components](https://gcalderone.github.io/GModelFit.jl/concepts/):

- `balmercont`: Balmer continuum (λ < 3645Å) and pseudo-continuum (i.e. unresolved Balmer emission lines at λ < 3645Å);
- `cutoff_powerlaw`: Continuum cutoff-powerlaw;
- `gaussconv`: Convolution of a spectrum sampled on a log-regular grid with a Gaussian kernel;
- `hostgalaxy`: Host galaxy templates;
- `interpolator`: Interpolate generic template on a wavelength grid;
- `ironopt`: Iron complex emission lines at optical wavelengths (from Véron-Cetty et al. 2004);
- `ironuv`: Iron complex emission lines at UV wavelengths (from Vestergaard & Wilkes 2001);
- `powerlaw`: Continuum powerlaw;
- `sbpl`: Continuum smoothly broken powerlaw;
- `SpecLineAsymmGauss`: Emission line with asymmetric Gaussian profile;
- `SpecLineGauss`: Emission line with Gaussian profile;
- `SpecLineLorentz`: Emission line with Lorentz profile;
- `SpecLineVoigt`: Emission line with Voigt profile;

To be written...
