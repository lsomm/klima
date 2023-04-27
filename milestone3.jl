import Plots
import DelimitedFiles
include("milestone2_solution_julia.jl")

geo_dat = read_geography("input/The_World128x65.dat")

function calc_area(n_latitude, n_longitude)

    delta_theta = pi/(n_latitude-1)
    area = zeros(Float64, n_latitude)
    area[1] = 0.5 * (1 - cos(0.5 * delta_theta))
    area[n_latitude] = area[1]

    for j=2:n_latitude-1
        theta_j = delta_theta * (j - 1)
        area[j] = sin(0.5 * delta_theta) * sin(theta_j) / n_longitude
    end

    return area
end

function calc_mean(field, area, n_latitude, n_longitude)
    if size(field)[1] !== n_latitude || size(field)[2] !== n_longitude
      error("field and area sizes do not match")
    end
  
    # Initialize mean with the values at the poles
    mean = area[1]*field[1,1] + area[end]*field[end,end]
  
    for j in 2:n_latitude-1
        for i in 1:n_longitude
            mean += area[j] * field[j,i]
        end 
    end
  
    return mean #/(4*pi*r(Radius der Erde Si-Einheit?))
end

function calc_radiative_cooling_co2(c ,c0 = 315.0)

    return 210.3 − 5.35 * log(c/c0)

end

function compute_equilibrium(T_0,c,C,S_sol, n_latitude, n_longitude,B=2.15)
    A = calc_radiative_cooling_co2(c)
    C_quer = calc_mean(C, calc_area(n_latitude,n_longitude),n_latitude,n_longitude)
    #S_sol_quer = calc_mean(S_sol,calc_area(n_latitude,n_longitude),n_latitude,n_longitude)
    # forward Euler method

    T = zeros(48)
    T[1] = T_0
    # secpy = 3.15576e7
    # true_longitude = readdlm("input/True_Longitude.dat")
    # delta_t = secpy/lengths(true_longitude)
    # timesteps = [delta_t*k for k in eachindex(true_longitude)]
    
    function f(T,t)
        S_sol_quer_t = calc_mean(S_sol[:,:,t],calc_area(n_latitude,n_longitude),n_latitude,n_longitude)
        return (S_sol_quer_t-A-B*T)/C_quer
    end

    for k = 1:47
        T[k+1] = T[k] + F(T[k],k)
    end
    return T

end