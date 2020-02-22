function poly(x, c)
    out = x .* 0. .+ c[1]
    for i in 2:length(c)
        out .+= c[i] .* x.^(i-1)
    end
    return out
end


function ccm_unred(wave, ebv; R_V=3.1)
    x = 10000 ./ wave    # Convert to inverse microns 
    npts = length( x )
    a = fill(0., npts)  
    b = fill(0., npts)  

    good = findall((x .>  0.3)  .&   (x .< 1.1)) #Infrared
    a[good] =  0.574 .* x[good].^(1.61)
    b[good] = -0.527 .* x[good].^(1.61)

    good = findall((x .>= 1.1)  .&   (x .< 3.3)) #Optical/NIR
    #Use new constants from O'Donnell (1994)
    y = x[good] .- 1.82
#     c1 = [ 1. , 0.17699, -0.50447, -0.02427,  0.72085,      #Original
#                 0.01979, -0.77530,  0.32999 ]               #coefficients
#     c2 = [ 0.,  1.41338,  2.28305,  1.07233, -5.38434,      #from CCM89
#                -0.62251,  5.30260, -2.09002 ]
      c1 = [ 1. , 0.104,   -0.609,    0.701,  1.137,         #New coefficients
                 -1.718,   -0.827,    1.647, -0.505 ]        #from O'Donnell
      c2 = [ 0.,  1.952,    2.908,   -3.989, -7.985,         #(1994)
                 11.102,    5.491,  -10.805,  3.347 ]
     a[good] = poly(y, c1)
     b[good] = poly(y, c2)

    good = findall((x .>= 3.3)  .&   (x .< 8.0)) #Mid-UV
    y = x[good]
    F_a = fill(0., length(good))
    F_b = fill(0., length(good))
    good1 = findall(y .> 5.9)
    y1 = y[good1] .- 5.9
    F_a[good1] = -0.04473 .* y1.^2 .- 0.009779 .* y1.^3
    F_b[good1] =   0.2130 .* y1.^2 .+   0.1207 .* y1.^3
    
    a[good] =  1.752 .- 0.316.*y .- (0.104 ./ ( (y.-4.67).^2 .+ 0.341 )) .+ F_a
    b[good] = -3.090 .+ 1.825.*y .+ (1.206 ./ ( (y.-4.62).^2 .+ 0.263 )) .+ F_b


    good = findall((x .>= 8.0)  .&   (x .< 11.)) #Far-UV
    y = x[good] .- 8.
    c1 = [ -1.073, -0.628,  0.137, -0.070 ]
    c2 = [ 13.670,  4.257, -0.420,  0.374 ]
    a[good] = poly(y, c1)
    b[good] = poly(y, c2)


    # Now apply extinction correction to input flux vector
    A_V = R_V * ebv
    A_lambda = A_V .* (a .+ b./R_V)
    return 10 .^(0.4.*A_lambda)
 end                               
