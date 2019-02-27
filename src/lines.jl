@quasiabstract struct Line
    label::String
    λ::Float64
    enabled::Bool
end

@quasiabstract struct EmissionLine <: Line
    fwhm::Float64
    fwhm_limits::NTuple{2, Float64}
    voff_limits::NTuple{2, Float64}
end

@quasiabstract mutable struct NarrowLine <: EmissionLine
    NarrowLine(label, λ) =
        new(label,        λ, true,   500, (100, 2.0e3), 1000 .* (-1,1))
end

@quasiabstract mutable struct CombinedNarrowLine <: EmissionLine
    CombinedNarrowLine(label, λ) =
        new(label,        λ, true,   500, (100, 1.0e3), 1000 .* (-1,1))
end

@quasiabstract mutable struct BroadLine <: EmissionLine
    BroadLine(label, λ) =
        new(label,        λ, true,  5000, (900, 1.5e4), 3000 .* (-1,1))
end

@quasiabstract mutable struct CombinedBroadLine <: EmissionLine
    CombinedBroadLine(label, λ) =
        new(label,        λ, true,  5000, (900, 1.5e4), 3000 .* (-1,1))
end

@quasiabstract mutable struct UnknownLine <: Line
    fwhm::Float64
    fwhm_limits::NTuple{2, Float64}
    λ_limits::NTuple{2, Float64}
    UnknownLine(λ) =
        new("unknown",    λ, false, 5000, (600, 1.0e4),  λ .+ 100 .* (-1,1))
end

@quasiabstract mutable struct AbsorptionLine <: UnknownLine
    AbsorptionLine(λ) =
        new("absorption", λ, false, 1000, (200, 3.0e4),  100 .* (-1,1))
end


function qsfit_lines()
    out = Vector{Line}()
    l = CombinedNarrowLine("na_Lyb"         , 1026.0  );  push!(out, l)
    l = CombinedBroadLine( "br_Lyb"         , 1026.0  );  push!(out, l)
    l = CombinedNarrowLine("na_Lya"         , 1215.24 );  push!(out, l)
    l = CombinedBroadLine( "br_Lya"         , 1215.24 );  push!(out, l)
    l = NarrowLine(        "na_NV_1241"     , 1240.81 );  push!(out, l)
    l = BroadLine(         "br_OI_1306"     , 1305.53 );  push!(out, l)
    l = BroadLine(         "br_CII_1335"    , 1335.31 );  push!(out, l)
    l = BroadLine(         "br_SiIV_1400"   , 1399.8  );  push!(out, l)
    l = BroadLine(         "br_CIV_1549"    , 1549.48 );  push!(out, l)
    # l = BroadLine(       "br_HeII"        , 1640.4  );  push!(out, l)
    # l = BroadLine(       "br_OIII"        , 1665.85 );  push!(out, l)
    # l = BroadLine(       "br_AlIII"       , 1857.4  );  push!(out, l)
    l = BroadLine(         "br_CIII_1909"   , 1908.734);  push!(out, l)
    l = BroadLine(         "br_CII"         , 2326.0  );  push!(out, l)
    # l = BroadLine(       "br_F2420"       , 2420.0  );  push!(out, l)
    l = BroadLine(         "br_MgII_2798"   , 2799.117);  l.fwhm_limits = 1000. .* (-1,1);    push!(out, l)
    # l = NarrowLine(      "na_NeV"N        , 3346.79 );  push!(out, l)
    l = NarrowLine(        "na_NeVI_3426"   , 3426.85 );  push!(out, l)
    l = NarrowLine(        "na_OII_3727"    , 3729.875);  push!(out, l)
    l = NarrowLine(        "na_NeIII_3869"  , 3869.81 );  push!(out, l)
    l = BroadLine(         "br_Hd"          , 4102.89 );  push!(out, l)
    l = BroadLine(         "br_Hg"          , 4341.68 );  push!(out, l)
    # l = BroadLine(       "br_HeII"        , ????    );  push!(out, l)
    l = CombinedNarrowLine("na_Hb"          , 4862.68 );  push!(out, l)
    l = CombinedBroadLine( "br_Hb"          , 4862.68 );  push!(out, l)
    l = NarrowLine(        "na_OIII_4959"   , 4960.295);  push!(out, l)
    l = NarrowLine(        "na_OIII_5007"   , 5008.240);  push!(out, l)
    l = NarrowLine(        "na_OIII_5007_bw", 5008.240);  l.fwhm = 500; l.fwhm_limits = (1e2, 1e3); l.voff_limits = (0, 2e3); push!(out, l)
    l = BroadLine(         "br_HeI_5876"    , 5877.30 );  push!(out, l)
    l = NarrowLine(        "na_NII_6549"    , 6549.86 );  push!(out, l)
    l = CombinedNarrowLine("na_Ha"          , 6564.61 );  push!(out, l)
    l = CombinedBroadLine( "br_Ha"          , 6564.61 );  push!(out, l)
    l = CombinedBroadLine( "br_Ha_base"     , 6564.61 );  l.fwhm = 2e4; l.fwhm_limits = (1e4, 3e4); push!(out, l)
    l = NarrowLine(        "na_NII_6583"    , 6585.27 );  push!(out, l)
    l = NarrowLine(        "na_SII_6716"    , 6718.29 );  push!(out, l)
    l = NarrowLine(        "na_SII_6731"    , 6732.67 );  push!(out, l)
end


function qsfit_lines_linkparameters(model::DataFitting.UI{Model})
    model.na_OIII_4959.voff.expr = "na_OIII_5007_voff"
end


function addEmLine(model, line::EmissionLine)
    cname = Symbol(line.label)
    addcomp!(model, cname => emline(line.λ))
    model[cname].fwhm.val  = line.fwhm
    model[cname].fwhm.low  = line.fwhm_limits[1]
    model[cname].fwhm.high = line.fwhm_limits[2]
    model[cname].voff.val  = 0.
    model[cname].voff.low  = line.voff_limits[1]
    model[cname].voff.high = line.voff_limits[2]
    return cname
end

function addEmLine(model, line::UnknownLine)
    cname = Symbol(line.label)
    addcomp!(model, cname => emline(line.λ))
    model[cname].fwhm.val  = line.fwhm
    model[cname].fwhm.low  = line.fwhm_limits[1]
    model[cname].fwhm.high = line.fwhm_limits[2]
    model[cname].center.low  = line.λ_limits[1]
    model[cname].center.high = line.λ_limits[2]
    model[cname].center.fixed = false
    model[cname].voff.fixed = true
    model[cname].voff.val = 0.
    return cname
end
