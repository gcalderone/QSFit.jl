using CMPFit, GFit, Gnuplot
using QSFit, DataStructures, Statistics, Dierckx

import QSFit: known_spectral_lines, add_qso_continuum!, LineComponent

export T2Recipe

abstract type T2Recipe <: DefaultRecipe end

function default_options(::Type{T}) where T <: T2Recipe
    out = OrderedDict{Symbol, Any}()
    out[:host_template] = "Ell5"
    out[:use_host_template] = true
    out[:use_balmer] = false
    out[:use_ironuv] = false
    out[:use_ironopt] = false
    out[:use_lorentzian_profiles] = false
    out[:n_unk] = 2
    out[:unk_avoid] = [4863 .+ [-1,1] .* 50, 
                       6565 .+ [-1,1] .* 150,
                       5008 .+ [-1,1] .* 25]
    return out
end

function known_spectral_lines(source::QSO{T}) where T <: T2Recipe
    list = [
        NarrowLine(                      :Lyb                       ),
        NarrowLine(                      :Lya                       ),
        NarrowLine(                      :NV_1241                   ),
        NarrowLine(                      :CIV_1549                  ),
        NarrowLine(                      :CIII_1909                 ),
        NarrowLine(                      :MgII_2798                 ),
        NarrowLine(                      :NeV_3426                  ),
        NarrowLine(                      :OII_3727                  ),
        NarrowLine(                      :NeIII_3869                ),
        NarrowLine(                      :Hg                        ),
        NarrowLine(                      :Hb                        ),
        NarrowLine(                      :OIII_4959                 ),
        NarrowLine(                      :OIII_5007                 ),
        AsymmTailLine(                   :OIII_5007      ,     :blue),
        NarrowLine(                      :OI_6300                   ),
        NarrowLine(                      :OI_6364                   ),
        NarrowLine(                      :NII_6549                  ),
        NarrowLine(                      :Ha                        ),
        NarrowLine(                      :NII_6583                  ),
        NarrowLine(                      :SII_6716                  ),
        NarrowLine(                      :SII_6731                  )]
    return list
end

function add_qso_continuum(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: T2Recipe
    λ = domain(model)[:]

    comp = QSFit.powerlaw(median(λ))
    comp.alpha.val = -1.8

    model[:qso_cont] = comp
    push!(model[:Continuum].list, :qso_cont)
    evaluate!(model)
end

function LineComponent(source::QSO{T}, line::NarrowLine, multicomp::Bool) where T <: T2Recipe
    lc = LineComponent(parent_recipe(source), line, multicomp) # invoke parent recipe
    lc.comp.fwhm.low  = 10 
    lc.comp.fwhm.high = 1000
    lc.comp.voff.high = 500
    @info line.tid lc.comp.fwhm.val
    return lc
end

function LineComponent(source::QSO{T}, line::AsymmTailLine, multicomp::Bool) where T <: T2Recipe
    comp = LineComponent(source, GenericLine(line.tid), multicomp).comp
    comp.fwhm.val  = 500
    comp.fwhm.high = 1e3
    comp.voff.low  = line.side == :blue ?    0  :  720
    comp.voff.high = line.side == :blue ?  720  :    0
    return LineComponent(line, comp, multicomp)
end

function default_unk_line(source::QSO{T}) where T <: T2Recipe
    comp = SpecLineGauss(5e3)
    comp.norm.val = 0.
    comp.center.fixed = false
    comp.center.low = 0
    comp.center.high = Inf
    comp.fwhm.val  = 500
    comp.fwhm.low  = 10
    comp.fwhm.high = 2e3
    comp.voff.fixed = true
    return comp
end

function add_patch_functs!(source::QSO{T}, pspec::PreparedSpectrum, model::Model) where T <: T2Recipe
    # Patch parameters
    @try_patch! begin
        # model[:OIII_4959].norm = model[:OIII_5007].norm / 3
        model[:OIII_4959].voff = model[:OIII_5007].voff
    end
    @try_patch! begin
        model[:OIII_5007_bw].voff += model[:OIII_5007].voff
        model[:OIII_5007_bw].fwhm += model[:OIII_5007].fwhm
        model[:OIII_5007_bw].norm += model[:OIII_5007].norm # is this correct way round?
    end
    @try_patch! begin
        # model[:OI_6300].norm = model[:OI_6364].norm / 3
        model[:OI_6300].voff = model[:OI_6364].voff
    end
    @try_patch! begin
        # model[:NII_6549].norm = model[:NII_6583].norm / 3
        model[:NII_6549].voff = model[:NII_6583].voff
    end
    @try_patch! begin
        # model[:SII_6716].norm = model[:SII_6731].norm / 1.5
        model[:SII_6716].voff = model[:SII_6731].voff
    end
end