using PlotlyJS
using DelimitedFiles


# TODO
# animation/gif/slider/schieberegler
# Nur Ränder plotten als Unterlage
# In Klimakoffer einbauen / Temperatur/etc plots damit machen
# Grid anzeigen?
# Welt höher aufgelöst
# Dokumentation

function plot_earth(surface_data, radius)

    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)


    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = geo',
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
            ]
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
            )


    PlotlyJS.plot([fig1,fig2], layout)   
    end

geo = readdlm("input/The_World128x65.dat")
plot_earth(geo, 10)