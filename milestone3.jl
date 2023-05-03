import Plots
import DelimitedFiles
using Infiltrator
include("milestone2.jl")

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

function compute_equilibrium(T_0,c,C,S_sol_average, n_latitude, n_longitude, B=2.15)
    A = calc_radiative_cooling_co2(c)
    C_quer = calc_mean(C, calc_area(n_latitude,n_longitude),n_latitude,n_longitude)
    #@info "" C_quer A

    T = zeros(49)
    T[1] = T_0
    T_implicit_euler = zeros(49)
    T_implicit_euler[1] = T_0
    
    function f(T,t,S_sol_average)
        
        result = (S_sol_average[t]-A-B*T)/C_quer
        #@infiltrate
        return result
    end

    function f_implicit(T,t,S_sol_average)

        result = (T + (S_sol_average[t] - A )/(48*C_quer) ) / (1+B/(48*C_quer) )
        return result
    end

    #@info "" f(T[2],2)
    #error()
    
    # forward Euler method
    for k = 1:48
        T[k+1] = T[k] + (1/48) *f(T[k],k,S_sol_average)
        T_implicit_euler[k+1] =  f_implicit(T_implicit_euler[k],k,S_sol_average)
    end
    return T#T_implicit_euler
     
   

end

function compare_temperatures(T_initial, steps,eps, c, C, S_sol_average, n_latitude, n_longitude)
    T_1 = zeros(49) 
    T_1[end] = T_initial
    for i in 1:steps
        T_2 = compute_equilibrium(T_1[end], c, C, S_sol_average, n_latitude, n_longitude)
        if maximum(abs.(T_1-T_2)) < eps #norm(T_1-T_2,inf)
            print("Reached equilibrium at step $(i) and temp $(T_1)!")
            print(T_2)
            return T_2
        end
        T_1 = copy(T_2)
    end
end

c = 315 #400
C = calc_heat_capacity(geo_dat)



S_sol = calc_solar_forcing(calc_albedo(geo_dat), read_true_longitude("input/True_Longitude.dat"))
#S_sol = readdlm("solar_forcing_averages.txt")
S_sol_quer = zeros(48)
for t = 1:48    
    S_sol_quer[t] = calc_mean(S_sol[:,:,t],calc_area(n_latitude,n_longitude),n_latitude,n_longitude)
end
n_latitude = 65
n_longitude = 128
T_initial = 5.0
steps = 50
eps = 10^(-5)
T_end = compare_temperatures(T_initial, steps, eps, c, C, S_sol_quer, n_latitude, n_longitude)
plot(T_end)

