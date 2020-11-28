@quasiabstract mutable struct SpectralLine
    λ::Float64
    fwhm::Float64
    fwhm_limits::NTuple{2, Float64}
end

@quasiabstract mutable struct KnownLine <: SpectralLine; end

@quasiabstract mutable struct KnownEmLine <: KnownLine
    both_brna::Bool
    voff_limits::NTuple{2, Float64}
end

@quasiabstract mutable struct NarrowLine   <: KnownEmLine;  end
@quasiabstract mutable struct BroadLine    <: KnownEmLine;  end

mutable struct CombinedLine <: KnownEmLine
    λ::Float64
    na::NarrowLine
    br::BroadLine
end

NarrowLine(λ) = NarrowLine(λ,  500, (100, 2.0e3), false, 1000 .* (-1,1))
BroadLine(λ)  = BroadLine( λ, 5000, (900, 1.5e4), false, 3000 .* (-1,1))
CombinedLine(λ) = CombinedLine(λ,
    NarrowLine( λ,  500, (100, 1.0e3), true, 1000 .* (-1,1)),
    BroadLine(  λ, 5000, (900, 1.5e4), true, 3000 .* (-1,1)))

@quasiabstract mutable struct UnknownLine <: SpectralLine
    λ_limits::NTuple{2, Float64}
end
UnknownLine() = UnknownLine(5e3, 5000, (600, 1.0e4), (0., Inf))


@quasiabstract mutable struct AbsorptionLine <: SpectralLine;
    λ_limits::NTuple{2, Float64}
end
AbsorptionLine(λ) = AbsorptionLine(0, 1000, (200, 3.0e4), 100 .* (-1,1))

function SpecLineGauss(line::T) where T <: Union{NarrowLine, BroadLine}
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = line.fwhm
    comp.fwhm.low  = line.fwhm_limits[1]
    comp.fwhm.high = line.fwhm_limits[2]
    comp.center.fixed = true
    comp.voff.fixed = false
    comp.voff.val  = 0.
    comp.voff.low  = line.voff_limits[1]
    comp.voff.high = line.voff_limits[2]
    return comp
end

function SpecLineGauss(line::UnknownLine)
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = line.fwhm
    comp.fwhm.low  = line.fwhm_limits[1]
    comp.fwhm.high = line.fwhm_limits[2]
    comp.center.fixed = false
    comp.center.low  = line.λ_limits[1]
    comp.center.high = line.λ_limits[2]
    comp.voff.fixed = true
    comp.voff.val = 0.
    return comp
end


function SpecLineGauss(line::AbsorptionLine)
    comp = SpecLineGauss(line.λ)
    comp.fwhm.val  = line.fwhm
    comp.fwhm.low  = line.fwhm_limits[1]
    comp.fwhm.high = line.fwhm_limits[2]
    comp.center.fixed = false
    comp.center.low  = line.λ_limits[1]
    comp.center.high = line.λ_limits[2]
    comp.voff.fixed = true
    comp.voff.val = 0.
    return comp
end


function tocomps(lines::OrderedDict{Symbol, SpectralLine}, prefix="")
    comps = OrderedDict{Symbol, SpecLineGauss}()
    for (cname, line) in lines
        if isa(line, CombinedLine)
            comps[Symbol(prefix, :na_, cname)] = SpecLineGauss(line.na)
            comps[Symbol(prefix, :br_, cname)] = SpecLineGauss(line.br)
        else
            comps[Symbol(prefix, cname)] = SpecLineGauss(line)
        end
    end
    return comps
end
