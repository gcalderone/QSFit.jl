function read_sdss_dr10(file::AbstractString)
    f = FITS(file)
    λ = 10 .^read(f[2], "loglam")
    flux = float.(read(f[2], "flux"))
    ivar = float.(read(f[2], "ivar"))
    mask = read(f[2], "and_mask")
    close(f)
    
    ndrop = 100
    λ = λ[ndrop:end-ndrop]
    flux = flux[ndrop:end-ndrop]
    ivar = ivar[ndrop:end-ndrop]
    mask = mask[ndrop:end-ndrop]

    igood =  findall((mask .== 0)  .&
                     (ivar .> 0)   .&
                     (flux .> 0))

    d = QSFitData(λ, flux, sqrt.(1 ./ ivar), igood, label="SDSS-DR10: " * file)
    return d
end
