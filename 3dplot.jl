using PlotlyJS
using DelimitedFiles



function plot_sphere(surface_data, radius)

    lat_resolution, long_resolution = size(surface_data)

    latitude = range(0, pi, lat_resolution)
    longitude = range(0, 2*pi, long_resolution)

    lat_grid = latitude'.*ones(long_resolution)
    long_grid = ones(lat_resolution)'.*longitude

    X = radius*sin.(lat_grid).*cos.(long_grid)
    Y = radius*sin.(lat_grid).*sin.(long_grid)
    Z = radius*cos.(lat_grid)


    fig = PlotlyJS.surface(x=X, y=Y, z=Z, surfacecolor=geo')
    
    layout = Layout(
                scene=attr(
                    xaxis=attr(
                        visible=false
                    ),
                    yaxis=attr(
                        visible=false
                    ),
                    zaxis=attr(
                        visible=false
                    )
                )
            )
    
    plot(fig, layout)   
    end

geo = readdlm("input/The_World128x65.dat")
plot_sphere(geo, 10)