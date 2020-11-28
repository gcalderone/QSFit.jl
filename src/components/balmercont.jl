function eval_balmer_pseudocont(Temp, Ne, fwhm)
    function readfile(file, Ne)
        lower = Vector{Float64}()
        upper = Vector{Float64}()
        norm  = Vector{Float64}()    
        density = NaN
        current_upper = 0
        for l in readlines(file)
            l = strip(l)
            while occursin(" =", l); l = replace(l, " =" => "="); end
            while occursin("= ", l); l = replace(l, "= " => "="); end
            s = split(l)
            if s[1] == "DENS"
                density = NaN
            elseif !isfinite(density)
                density = parse(Float64, s[1])
            elseif density == Ne
                if (length(s[1]) >= 4)  &&  (s[1][1:4] == "E_NU")
                    current_upper = parse(Int, split(s[1],'=')[2])
                elseif occursin(r"\D+", s[1])
                    current_upper = 0
                elseif (occursin(r"\d+", s[1]))  &&  (current_upper > 0)
                    @assert mod(length(s), 2) == 0
                    for ii in 1:2:length(s)
                        push!(upper, current_upper)
                        push!(lower, parse(Int    , s[ii]))
                        push!(norm , parse(Float64, s[ii+1]))
                    end
                end
            end
        end
        norm ./= maximum(norm)
        return (lower, upper, norm)
    end

    # List of available temperatures in SH95 [K/100]
    Tavail = [5, 10, 30, 50, 75, 100, 125, 150, 200, 300]
    iTemp = argmin(abs.(Temp./100. .- Tavail))
    if Temp != Tavail[iTemp]*100 
        @info "Balmer pseudo-continuum: T (requested)=$Temp, T (available) = $(Tavail[iTemp]*100)"
    end

    # List of available electron densities in SH95 [log (Ne / cm^-3)]
    logNeAvail = [2., 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    iNe = argmin(abs.(log10(Ne) .- logNeAvail))
    if Ne != 10^logNeAvail[iNe]
        @info "Balmer pseudo-continuum: Ne (requested)=$Ne, Ne (available) = $(10^logNeAvail[iNe])"
    end
    
    # Read appropriate file
    path = qsfit_data() * "/SH1995/"
    file = path * @sprintf("r1b%04d.d", Tavail[iTemp])
    electronDensity = 10^logNeAvail[iNe]

    (lower, upper, norm) = readfile(file, electronDensity)
    
    ryd = 911.2671808036101 # Rydberg in A
    edge = 4 * ryd          # Balmer (2^2=4) edge in Angstrom

    wave = ryd ./ (1 ./ lower.^2 - 1 ./ upper.^2)

    ii = findall((lower .== 2)  .&  (upper .> 6))
    wave = wave[ii]
    norm = norm[ii]

    # Prepare output variables
    λ = 10 .^range(log10(3400.), stop=log10(4300), length=500) # Angstrom
    cont = fill(0., length(λ))       # Lum. density (arbitrary units)
    
    # Broadening of high order Balmer lines ==> pseudo-continuum
    σ = (fwhm / 3.e5) .* wave ./ 2.35
    for i in 1:length(wave)
        cont += norm[i] * gauss.(λ, wave[i], σ[i])
    end
    cont ./= interpol(cont, λ, edge)
    return (λ, cont)
end


function eval_balmer_continuum(Temp, Tau, fwhm)
    ryd = 911.2671808036101 # Rydberg in A
    edge = 4 * ryd          # Balmer (2^2=4) edge in Angstrom

    λ = 10 .^range(log10(912.), stop=log10(5000), length=600) .* 1e-8 # cm
    
    # Planck function
    b = planck(λ, Temp)
  
    # Take into account optical depth
    λ .*= 1e8
    cont = b .* (1. .- exp.(-(Tau .* (λ ./ edge).^3.)))
    cont[findall(λ .>= edge)] .= 0.
    cont ./= maximum(cont)

    # Broadening
    σ = fwhm / 3.e5 / 2.35 * edge
    σ = 2.51705 * 6 # TODO: this is to match IDL results, above formula is the correct one
    kernel = gauss(λ, mean(λ), σ)
    kernel = kernel[findall(kernel .> maximum(kernel) / 1e3)]
    (mod(length(kernel), 2) == 0)  &&  (push!(kernel, 0.))
    cont = convol(cont, kernel)

    # Normalize Balmer continuum at 3000A
    cont ./= interpol1(cont, λ, 3000.)
    contAtEdge = interpol1(cont, λ, edge)
    return (λ, cont, contAtEdge)
end


# ____________________________________________________________________
# balmercont
#
mutable struct balmercont <: AbstractComponent
    norm::Parameter
    ratio::Parameter
    function balmercont(norm::Number, ratio::Number)
        out = new(Parameter(norm), Parameter(ratio))
        out.norm.low = 0
        out.ratio.low = 0        
        return out
    end
end

mutable struct balmercont_cdata
    c1::Vector{Float64}
    c2::Vector{Float64}
end

function compeval_cdata(comp::balmercont, domain::Domain_1D)
    T = 15000.
    Ne = 1e9
    Tau = 1.
    fwhm = 5040.8 # 5000.
    (λ1, c1, contAtEdge) = eval_balmer_continuum(T, Tau, fwhm)
    (λ2, c2) = eval_balmer_pseudocont(T, Ne, fwhm)
    c2 .*= contAtEdge
    out = balmercont_cdata(interpol1(c1, λ1, domain[1]),
                           interpol1(c2, λ2, domain[1]))
    out.c1[findall(out.c1 .< 0.)] .= 0.
    out.c2[findall(out.c2 .< 0.)] .= 0.
    return out
end
compeval_array(comp::balmercont, domain::Domain_1D) = fill(NaN, length(domain))

function evaluate(c::CompEval{balmercont, Domain_1D},
                  norm, ratio)
    c.buffer .= norm .* (c.cdata.c1 .+ ratio .* c.cdata.c2)
end



#=
m = Model(:c => balmercont(0.1, 0.5))
n = DataFitting.wrappee(m)
add_dom!(m, 1000.:4500)
addexpr!(m, :c)
evaluate!(m)

@gp "set grid" yr=(0, 1.1) :-
@gp :- xlab="Wavelength [A]" ylab="Lum. density [arb.units]" :-
@gp :- 3645.07.*[1,1] [0, 1.1] "w l dt 4 lc rgb 'black'" :-
m.c.norm.val = 1
m.c.ratio.val = 1
evaluate!(m)
@gp :- m(:domain) m() "w l lw 2 lc rgb 'red'"
m.c.ratio.val = 0.5
evaluate!(m)
@gp :- m(:domain) m() "w l lw 2 lc rgb 'blue'"
=#
