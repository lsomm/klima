using PlotlyJS
using DelimitedFiles

include("milestone2.jl")
# TODO
# animation/gif/slider/schieberegler
# In Klimakoffer einbauen / Temperatur/etc plots damit machen
# Grid anzeigen?
# Welt höher aufgelöst
# Dokumentation
# Scientific programming standard (s. Website)
# Himmel/Sterne/Weltall anzeigen
# Solar mechanics ?

function parametrisierung(surface_data, radius=10)

    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)

    return X,Y,Z

end

function plot_earth(X,Y,Z,surface_data)
    #=
    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)
    =#

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        colorscale = [
            [0, "rgb(60,149,33)"], #land, green
            [0.25, "rgb(165,255,255)"], #sea ice, icy blue
            [0.5,"rgb(255,255,255)"], #snow, white
            [1,"rgb(26,100,212)" ], #ocean, blue
        ],
        colorbar = (
            autotick = false, 
            tick0 = 0,
            dtick = 1,
            tickcolor = 888,
            tickvals = [1,2,3,4,5],
            ticktext = [
                "land",
                "sea ice",
                "snow cover",
                "lakes",
                "ocean"
            ],
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
            )
        )
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        contours = attr(
            geo = attr(
                show = true,
                size = 0.05
            )
        ),
        showscale = false,
        opacity = 0.0,
    )
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black"
            )


    PlotlyJS.plot([fig1,fig2], layout)   
end

function plot_data_field(X,Y,Z, data_field, surface_data)
    # function plot_data_field(data_field, surface_data, radius=10)
    #=
    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)
    =#

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = data_field',
        showscale = true,
        opacity = 0.5,
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black"
            )


    PlotlyJS.plot([fig1,fig2], layout)   

end

function get_outlines(geo_dat)
    outlines = zeros(size(geo_dat))
    for i in 2:64 # Edge cases can be ignored because there is no continent 
        for j in 2:127
            if geo[i,j] == 1
                if ((geo[i,j-1] != 1) || (geo[i,j+1] != 1) || (geo[i-1,j] != 1) || (geo[i+1,j] != 1))
                    outlines[i,j] = 1
                end
            end
        end
    end
    return round.(Int, outlines)
end

function plot_albedo_3d(X,Y,Z, albedo, surface_data)
    # function plot_albedo_3d(albedo, surface_data, radius=10)

    #=
    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)
    =#

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = albedo',
        showscale = true,
        opacity = 0.5,
        colorscale = [
            [0, "rgb(0,0,0)"], # black, nothing reflected
            [1,"rgb(255,255,255)" ], # white, everything reflected
        ],
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
                ),
            title="albedo",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 15
                )   
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black"
            )


    PlotlyJS.plot([fig1,fig2], layout)   

end

function plot_heatcapacity_3d(X,Y,Z, heat_capacity, surface_data)
    # function plot_heatcapacity_3d(heat_capacity, surface_data, radius=10)

    #=
    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)
    =#

    heat_capacity_log = log10.(heat_capacity)

    vmin = minimum(heat_capacity_log)
    vmax = maximum(heat_capacity_log)

    cbar_ticks = LinRange(vmin, vmax, 10)
    cbar_tick_values = [@sprintf("%.3f", 10^tick) for tick in cbar_ticks]

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        showscale = false,
        colorscale = [
            [0, "rgb(255,255,255)"], # no border
            [1,"rgb(0,0,0)" ], # continent border
        ],
    )
    
    fig2 = PlotlyJS.surface(
        x = 1.01 .* X,
        y = 1.01 .* Y,
        z = 1.01 .* Z,
        surfacecolor = heat_capacity_log',
        showscale = true,
        opacity = 0.8,
        colorscale =[
            [0, "rgb(255,255,255)" ], # white
            [1, "rgb(139,0,0)"], # red 
        ],
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 15
                ),
            tickvals = cbar_ticks,
            ticktext = cbar_tick_values,
            title="heat capacity",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 15
                )   
        )
    )
    
    layout = Layout(
                scene = attr(
                    xaxis = attr(
                        visible = false
                    ),
                    yaxis = attr(
                        visible = false
                    ),
                    zaxis = attr(
                        visible = false
                    ),
                ),
                paper_bgcolor = "black"
            )


    PlotlyJS.plot([fig1,fig2], layout)   

end

geo = readdlm("input/The_World128x65.dat")

albedo = calc_albedo(geo)
heat_capacity = calc_heat_capacity(geo)

#plot_data_field(fake_temps, get_outlines(geo),10)

# Parametrisierung
X,Y,Z = parametrisierung(geo)
display(plot_earth(X,Y,Z,geo)) # earth
display(plot_albedo_3d(X,Y,Z, albedo, get_outlines(geo))) # albedo
plot_heatcapacity_3d(X,Y,Z, heat_capacity, get_outlines(geo)) # heat capacity

