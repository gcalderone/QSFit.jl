# Reference: http://www.astro.wisc.edu/~dolan/constants.html
@with_kw_noshow struct gpc
    # Dimensionless or common to both CGS and MKS
    jd_mjd    = 2400000.5              # Offset to change from Julian date to MJD
    NA        = 6.0221367    * 1e23    # Avagadro's number
    α         = 7.29735308   * 1e-3    # Fine structure constant
    day       = 86400.                 # days to sec fator
    month     = 2592000.               # months to sec fator
    year      = 31536000.              # year to sec fator

    # CGS         
    c         = 2.99792458  * 1e10     * u"cm" / u"s"                          # Vacuum speed of light
    h         = 6.6260755   * 1e-27    * u"erg" * u"s"                         # Planck's constant
    ħ         = 1.05457266  * 1e-27    * u"erg" * u"s"                         # Reduced Planck's constant
    G         = 6.67259     * 1e-8     * u"cm"^3 * u"g"^-1 * u"s"^-2           # Gravitational constant 
    e         = 4.8032068   * 1e-10                                            # Electron charge (esu)
    me        = 9.1093897   * 1e-28    * u"g"                                  # Mass of electron
    mp        = 1.6726231   * 1e-24    * u"g"                                  # Mass of proton
    mn        = 1.6749286   * 1e-24    * u"g"                                  # Mass of neutron
    mH        = 1.6733      * 1e-24    * u"g"                                  # Mass of hydrogen
    amu       = 1.6605402   * 1e-24    * u"g"                                  # Atomic mass unit
    k         = 1.380658    * 1e-16    * u"erg" / u"K"                         # Boltzmann constant
    eV        = 1.6021772   * 1e-12    * u"erg" / u"eV"                        # Electron volt to erg factor
    a         = 7.5646      * 1e-15    * u"erg" * u"cm"^-3 * u"K"^-4           # Radiation density constant
    σ         = 5.67051     * 1e-5     * u"erg" * u"cm"^-2 * u"K"^-4 * u"s"^-1 # Stefan-Boltzmann constant
    R_inf     = 1.097373    * 1e5      * u"cm"^-1                              # R_infinity
    au        = 1.496       * 1e13     * u"cm" / u"AU"                         # Astronomical unit to cm
    pc        = 3.086       * 1e18     * u"cm" / u"pc"                         # Parsec to cm factor
    ly        = 9.463       * 1e17     * u"cm" / u"ly"                         # Light year to cm factor
    ld        = 2.5902      * 1e15                                             # Light day to cm factor (cm ly^-1)
    Msun      = 1.99        * 1e33     * u"g"                                  # Solar mass
    Rsun      = 6.96        * 1e10     * u"cm"                                 # Solar radius
    Lsun      = 3.9         * 1e33     * u"erg" / u"s"                         # Solar luminosity
    Tsun      = 5.780       * 1e3      * u"K"                                  # Solar temperature
    Mearth    = 5.974       * 1e27     * u"g"                                  # Earth mass
    Rearth    = 6372.8      * 1e5      * u"cm"                                 # Earth mean radius
    g         = 981.52                 * u"cm" * u"s"^-2                       # Acceleration at Earth surface
    thom      = 0.66524616  * 1e-24    * u"cm"^2                               # Thomson cross section
    jansky    = 1e-23                  * u"erg" * u"cm"^-2 * u"s"^-1 * u"Hz"^-1# Flux density (keV cm^-2 s^-1 keV^/1)
    wien      = 2.82 * k / h           
    deg       = pi/180                 * u"rad" * u"°"^-1
    arcsec    = pi/180/3600                                                    # rad arcsec^-1
    mas       = pi/180/3600/1000                                               # rad milliarcsec^-1
    edd       = 1.26e38                * u"erg" * u"s"^-1 * u"Msun"^-1         # Eddington luminosity (erg s^-1 M_sun^-1)
    r_cm      = 2 * pi^2 * me * e^4 / (h^3) / c                                # Rydberg constant (cm^-1)
    r2_cm     = r_cm / (1 + me/mp)                                             # Rydberg constant (cm^-1, reduced mass)
    A         = 1.e-8                  * u"cm" / UnitfulAstro.angstrom
end

function qsfitversion()
    return v"0.0.1"
end

function qsfitpath()
    #dirname(pathof(QSFIT))
    (file, line) = functionloc(qsfitversion, ())
    return dirname(file)
end


function gauss(x, μ, σ)
    return exp.(-0.5 .* ((x .- μ) ./ σ).^2) ./ sqrt(2pi) ./ σ
end

function planck(λ, T)
    c = gpc()
    b = 2 * c.h.val * c.c.val^2. ./ (λ.^5.)
    d = c.h.val * c.c.val ./ (λ .* c.k.val .* T)
    return b ./ (exp.(d) .- 1)
end
  
function convol(v, _k)
    k = reverse(_k)
    nk = length(k)
    @assert nk < length(v)
    @assert mod(nk, 2) == 1
    r = div(nk, 2)
    out = fill(0., length(v))
    for i in r+1:length(v)-r
        out[i-r:i+r] .+= v[i] .* k
    end
    return out ./ sum(abs.(k))
end

function interpol(y, x, X)
    out = fill(0., length(X))
    itp = interpolate((x,), y, Gridded(Linear()))
    ii = findall(minimum(x) .<= X .<= maximum(x))
    out[ii] .= itp(X[ii])
    return out
end
interpol(y, x, X::Number) = interpol(y, x, [X])[1]

function interpol1(y, x, X)
    out = fill(0., length(X))
    itp = interpolate((x,), y, Gridded(Linear()))
    ii = findall(minimum(x) .<= X .<= maximum(x))
    out[ii] .= itp(X[ii])
    etp = extrapolate(itp, Interpolations.Line())
    ii = findall(minimum(x) .> X)
    out[ii] .= etp(X[ii])
    ii = findall(maximum(x) .< X)
    out[ii] .= etp(X[ii])
    return out
end
interpol1(y, x, X::Number) = interpol1(y, x, [X])[1]



function smooth(y, n)
    out = y .* 1.
    @assert mod(n, 2) == 1
    @assert n >= 3
    h = div(n-1, 2)
    for i in 1+h:length(y)-h
        out[i] = mean(y[i-h:i+h])
    end
    return out
end

#=
function boole_int(x, f)
    @assert length(x) == length(f)
    N = 5
    Xmin = minimum(x)
    Xmax = maximum(x)
    X = range(Xmin, stop=Xmax, length=N)
    h = X[2] - X[1]
    
    itp = interpolate((x,), f, Gridded(Linear()))
    F = collect(itp(X))
    @assert length(F) == N
    #@gp    x f "w lp" :-
    #@gp :- X F "w lp"

    #t1 = range(0., stop=1, length=length(x))
    #A = hcat(x,f)
    #itp = scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))), NoInterp())), t1, 1:2)
    #t2 = range(0., stop=1, length=N)
    #X, F = [itp(t,1) for t in t2], [itp(t,2) for t in t2]
    #@assert length(F) == N
    #@gp :- X F "w lp"

    return 2*h / 45 * (7  * F[1] +
                       32 * F[2] +
                       12 * F[3] + 
                       32 * F[4] +
                       7  * F[5])
end
=#
