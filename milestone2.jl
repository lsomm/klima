using DelimitedFiles
using Plots
pythonplot()

# Creates an nlatitude x nlongitude = 65 x 128 array with integer digits decoding the geography.
function read_geography(filepath)
    return readdlm(filepath)
end

function robinson_projection(nlatitude, nlongitude)
    function x_fun(lon, lat)
        return lon / pi * (0.0379 * lat^6 - 0.15 * lat^4 - 0.367 * lat^2 + 2.666)
    end

    function y_fun(_, lat)
        return 0.96047 * lat - 0.00857 * sign(lat) * abs(lat)^6.41
    end

    # Longitude goes from -pi to pi (not included), latitude from pi/2 to -pi/2.
    # Remove endpoint of longitude to avoid overlap.
    x_lon = LinRange(-pi, pi, nlongitude + 1)[1:(end - 1)]

    # Latitude goes backwards because our plot starts at the bottom,
    # but the matrix is read from the top.
    y_lat = LinRange(pi / 2, -pi / 2, nlatitude)

    x = [x_fun(lon, lat) for lat in y_lat, lon in x_lon]
    y = [y_fun(lon, lat) for lat in y_lat, lon in x_lon]

    return x, y
end

function plot_geo(geo_dat)
    nlatitude, nlongitude = size(geo_dat)
    x, y = robinson_projection(nlatitude, nlongitude)

    # Unfortunately, this is a bit buggy. When passing a `cgrad` with `categorical=true`,
    # one would expect to only get the colors of this discrete `cgrad`, but for some reason,
    # that is not the case.
    # Instead, we got the colors that we want by experimenting with the levels.
    # We tried to make the range for ocean only slightly larger than the others to avoid
    # a weird looking colorbar.
    plot = contourf(x, y, geo_dat,
                    levels=[0.5, 1.7, 2.9, 4.1, 5.5],
                    clims=(1, 5),
                    aspect_ratio=1,
                    title="Earth Geography",
                    c=cgrad([:darkgreen, :lightsteelblue, :lavender, :navy]),
                    colorbar_ticks=([1.1, 2.3, 3.5, 4.8],
                                    ["land", "sea ice", "snow cover", "ocean"]),
                    axis=([], false),
                    dpi=300)

    return plot
end

function calc_albedo(nlatitude,nlongitude)
    surface_type = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    lat = LinRange(pi / 2, -pi / 2, nlatitude)
    #lon = LinRange(-pi, pi, nlongitude + 1)[1:(end - 1)]
    # 1 - land
    # 2 - sea ice
    # 3 - snow
    # 5 - ocean
    function p2(i)
        return  (3*(sin(lat[i])^2)-1)/2
    end
    function albedo(i,j,type)
        if type[i,j]==1
            alb = 0.3 + 0.12*p2(i)
        elseif type[i,j]==2
            alb = 0.6
        elseif type[i,j] ==3
            alb = 0.75
        elseif type[i,j] == 5
            alb = 0.29 + 0.12*p2(i)
        end
        return alb
    end
    #print(albedo(17,2,surface_type))


    return [albedo(i,j,surface_type) for i in 1:nlatitude, j in 1:nlongitude]
end

function calc_heatcapacity(nlatitude,nlongitude)
    surface_type = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    #[atm,land, sea_ice,snow,ocn_mxd_lyr]
    rho = [1.225,1350,917,400,1030]
    C_spec = [1000,750,2000,880,4000]
    d = [3850,1.0,1.5,0.5,70]
    C = rho .* C_spec .* d
    C_atm = C[1]
    function heatcapacity(i,j,type)
        if type[i,j]==1
            cap = C[2]
        elseif type[i,j]==2
            cap = C[3]
        elseif type[i,j]==3
            cap=C[4]
        elseif type[i,j]==5
            cap= C[5]
        end
        return cap
    end
    C_matrix = [heatcapacity(i,j,surface_type) for i in 1:nlatitude, j in 1:nlongitude]
    C_matrix = C_matrix .+ C_atm
    C_matrix = C_matrix ./ (3.15576 * 10^7)
    return C_matrix
end

function read_true_longitude()
    return readdlm(joinpath(@__DIR__, "input", "True_Longitude.dat"))
end

function insolation(lat, S_0, eccen, obliq, prec_dist, true_lon)
    delta = asin(sin(obliq)*sin(true_lon))
    rho = ((1-eccen*cos(true_lon-prec_dist))/(1-eccen^2))^2
    if -(pi/2-abs(delta))<lat<pi/2-abs(delta)
        H_0 = acos(-tan(lat)tan(delta))
        return S_0*rho/pi*(H_0*sin(lat)*sin(delta)+cos(lat*cos(delta)*sin(H_0)))
    elseif lat*delta >0
        return S_0*rho*sin(lat)*sin(delta)
    else
        return 0
    end
end

function calc_solar_forcing(albedo, true_lon_list, ntimesteps, nlatitude, nlongitude, S_0=1371.685, eccen = 0.01674, obliq = 0.409253, prec_dist = 1.783037)
    lat = LinRange(pi / 2, -pi / 2, nlatitude)
    return [insolation(lat[i],S_0, eccen, obliq, prec_dist, true_lon_list[t]) * (1-albedo[i,j]) for i in 1:nlatitude, j in 1:nlongitude, t in 1:ntimesteps]
end

function plot_data(data)
    nlatitude, nlongitude = size(data)
    x, y = robinson_projection(nlatitude, nlongitude)
    plot = contourf(x,y,data, levels = 200)
    return plot
end
# Run code
function milestone1()
    geo_dat = read_geography(joinpath(@__DIR__, "input", "The_World128x65.dat"))
    plot_geo(geo_dat)
end

#read_geography("/home/luca/klima/input/The_World128x65.dat")
albedo = calc_albedo(65,128)
heat_capacity = calc_heatcapacity(65,128)
solar_forcing = calc_solar_forcing(albedo,read_true_longitude(),48,65,128)


#plot_data(albedo)
#plot_data(heat_capacity)
#plot_data(solar_forcing[:,:,1])

#anim = @animate for t in 1:48
#    plot_data(solar_forcing[:,:,t])
#end
#gif(anim, "anim_fps15.gif", fps = 15)
