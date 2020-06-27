@quasiabstract struct SpectralLine
    label::String
    λ::Float64
    enabled::Bool
    fwhm::Float64
    fwhm_limits::NTuple{2, Float64}
end

@quasiabstract struct KnownLine <: SpectralLine
    voff_limits::NTuple{2, Float64}
end

@quasiabstract mutable struct UnknownLine <: SpectralLine
    λ_limits::NTuple{2, Float64}
end

@quasiabstract struct KnownEmLine  <: KnownLine; end
@quasiabstract struct KnownAbsLine <: KnownLine; end

@quasiabstract mutable struct NarrowLine         <: KnownEmLine;  end
@quasiabstract mutable struct BroadLine          <: KnownEmLine;  end
@quasiabstract mutable struct CombinedNarrowLine <: KnownEmLine;  end
@quasiabstract mutable struct CombinedBroadLine  <: KnownEmLine;  end
@quasiabstract mutable struct AbsorptionLine     <: KnownAbsLine; end

NarrowLine(label, λ)         = NarrowLine(label,         λ, true,   500, (100, 2.0e3), 1000 .* (-1,1))
CombinedNarrowLine(label, λ) = CombinedNarrowLine(label, λ, true,   500, (100, 1.0e3), 1000 .* (-1,1))
BroadLine(label, λ)          = BroadLine(label,          λ, true,  5000, (900, 1.5e4), 3000 .* (-1,1))
CombinedBroadLine(label, λ)  = CombinedBroadLine(label,  λ, true,  5000, (900, 1.5e4), 3000 .* (-1,1))
AbsorptionLine(label, λ)     = AbsorptionLine(label,     λ, false, 1000, (200, 3.0e4),  100 .* (-1,1))
UnknownLine(label, λ)        = UnknownLine(label,        λ, false, 5000, (600, 1.0e4),  λ .+ 100 .* (-1,1))


function emline(line::SpectralLine)
    comp = emline(line.λ)
    comp.fwhm.val  = line.fwhm
    comp.fwhm.low  = line.fwhm_limits[1]
    comp.fwhm.high = line.fwhm_limits[2]
    if isa(line, KnownLine)
        comp.center.free = false
        comp.voff.free = true
        comp.voff.val  = 0.
        comp.voff.low  = line.voff_limits[1]
        comp.voff.high = line.voff_limits[2]
    end
    if  isa(line, UnknownLine)
        comp.center.free = true
        comp.center.low  = line.λ_limits[1]
        comp.center.high = line.λ_limits[2]
        comp.voff.free = false
        comp.voff.val = 0.
    end
    return comp
end
