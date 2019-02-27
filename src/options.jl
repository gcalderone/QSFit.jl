@with_kw mutable struct QSFitOptions
    # The minimum wavelength used during fit.  Smaller wavelengths are
    # ignored.
    min_wavelength = 1210.

    # Fraction of negative residuals after continuum re-normalization
    cont_negative_fraction = 0.9

    # Value for the continuum.alpha1 parameter when z<=alpha1_fixed_max_z
    alpha1_fixed_value = -1.7
        
    # Max redshift to keep continuum.alpha1 fixed.  Beyond this
    # redshift the parameter is free to vary.  This functionality
    # allows to avoid degeneracy with the host galaxy template.
    alpha1_fixed_z = 0.6

    # Flag to use the host galaxy component
    use_galaxy = true

    # Name of host galaxy template
    galaxy_template = "Ell5"

    # Max redshift to use the galaxy template
    galaxy_enabled_z = 0.8

    # Flag to use the Balmer continuum component
    use_balmer = true

    # Min redshift to keep the Balmer component fixed.
    balmer_fixed_z = 1.1

    # Flag to use the iron UV component
    use_ironuv = true

    # Flag to use the iron optical component
    use_ironopt = true

    # If true use Lorentzian (rather than Gaussian) profiles for
    # emission lines
    lorentzian = false
    
    # If true tie the FWHM of broad component to be larger than the
    # FWHM of the associated narrow line
    bn_Fwhmtied = false

    # If true use a further emission line component for the blue
    # wing of the [OIII]5007 lne.
    oiii5007_bluewing = true

    # The number of unknown lines whose center wavelength is not
    # a-priori assigned: they are placed (after all other emission
    # lines have been fitted) at wavelengths where the largest fitting
    # residuals occur.
    unkLines = 10
end
