using PlotlyJS
using DelimitedFiles
using LaTeXStrings

include("milestone2.jl")
# TODO
# Welt höher aufgelöst 
# Dokumentation
# Scientific programming standard (s. Website)
# Himmel/Sterne/Weltall anzeigen
# Solar mechanics ?
# Überschriften 
# difussion coefficient plot
# Bedienungsanleitung (insb. wir können andere unser benutzen?)
# Diskrete Skala für Erde
# Am Anfang einmal transponieren, als immer
# Unterliegende Erde funktion schreiben
# Opacity vereinheitlichen?
# Latex strings?

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

    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = surface_data',
        colorscale = [
            [0, "rgb(60,149,33)"], #land, green
            [0.25, "rgb(165,255,255)"], #sea ice, icy blue
            [0.5, "rgb(255,255,255)"], #snow, white
            [1, "rgb(26,100,212)" ], #ocean, blue
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
                paper_bgcolor = "black",
                title = attr(
                    text = "The Earth",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
    )

    return PlotlyJS.plot(fig1, layout)   
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
            #=    
            title="albedo",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 15
            ) =#  
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
                paper_bgcolor = "black",
                title = attr(
                    text = "Surface albedo",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
            )


    return PlotlyJS.plot([fig1,fig2], layout)   

end

function plot_heatcapacity_3d(X,Y,Z, heat_capacity, surface_data)

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
        opacity = 0.7,
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
            title="[ J / m² / K ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
                ),
            titleside = "right"   
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
                paper_bgcolor = "black",
                title = attr(
                    text = "Heat capacity",
                    y=0.95,
                    x=0.5
                ),
                titlefont = (
                    color = "rgb(255,255,255)",
                    size = 40
                ),
    )

    return PlotlyJS.plot([fig1,fig2], layout)   

end

function plot_solar_forcing_3d_anim(X,Y,Z,solar_forcing, surface_data)

    cmin, cmax = extrema(solar_forcing)
    n_frames = size(solar_forcing, 3) #Länge von true longs, aber dartauf können wir hier nicht zugreifen, deswegen solar
    
    fig1 = PlotlyJS.surface(
        x = X,
        y = Y,
        z = Z,
        surfacecolor = get_outlines(surface_data)',
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
        surfacecolor = solar_forcing[:,:,1]',
        showscale = true,
        opacity = 0.8,
        colorscale = "Hot",
        cmax = cmax,
        cmin = cmin,
        colorbar = (
            autotick = false, 
            tickcolor = 888,
            tickfont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            title="[ W / m² ]",
            titlefont = (
                color = "rgb(255,255,255)",
                size = 20
            ),
            titleside = "right"   
        )
    )

    

    sliders = [
        attr(
            steps = [
                attr(
                    method = "animate",
                    args= [
                        [
                            "fr$k"
                        ],                           
                        attr(
                            mode = "immediate",
                            frame = attr(
                                duration = 40,
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            )
                        )
                    ],
                    label = "$k"
                )
            for k in 1:n_frames], 
            active = 17,
            transition = attr(
                duration = 0
            ),
            x=0, # slider starting position  
            y=0, 
            currentvalue = attr(
                font = attr(
                        size = 12
                    ), 
                prefix = "Step: ", 
                visible = true, 
                xanchor = "center"
            ),  
            len = 1.0 #slider length
        )    
    ];
    
    updatemenus = [
        attr(
            type = "buttons", 
            active = 0,
            y = 0.0,  #(x,y) button position 
            x = 1,
            buttons = [
                attr(
                    label = "Play",
                    method = "animate",
                    args = [
                        nothing,
                        attr(
                            frame = attr(
                                duration = 5, 
                                redraw = true
                            ),
                            transition = attr(
                                duration = 0
                            ),
                            fromcurrent = true,
                            mode = "immediate"
                        )
                    ]
                )
            ]
        )
    ];



    frames  = Vector{PlotlyFrame}(undef, n_frames)
    for k in 1:n_frames
        day = (round(Int, (k - 1) / n_frames * 365) + 80) % 365
        frames[k] = PlotlyJS.frame(
                        data = [
                            attr(
                                surfacecolor = solar_forcing[:,:,k]',            
                                colorbar = attr(
                                    #title = "Day $day",
                                )
                            )
                        ],
                        layout = attr(
                            title_text = "Solar forcing - Day $day"
                        ),
                        name="fr$k",
                        traces=[1],
                    )
    end    

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
        paper_bgcolor = "black",
        sliders = sliders,
        title = attr(
            text = "Solar forcing",
            y=0.95,
            x=0.5
        ),
        titlefont = (
            color = "rgb(255,255,255)",
            size = 40
            ),
        updatemenus = updatemenus
    )

    return Plot([fig1,fig2], layout, frames)
end

geo = readdlm("input/The_World128x65.dat")

albedo = calc_albedo(geo)
heat_capacity = calc_heat_capacity(geo)
true_lon = read_true_longitude("input/True_Longitude.dat")
solar_forcing = calc_solar_forcing(albedo,true_lon)
X,Y,Z = parametrisierung(geo)
#display(plot_earth(X,Y,Z,geo)) # earth
#display(plot_albedo_3d(X,Y,Z, albedo, get_outlines(geo)))# albedo
#display(plot_heatcapacity_3d(X,Y,Z, heat_capacity, get_outlines(geo))) # heat capacity



display(plot_solar_forcing_3d_anim(X,Y,Z,solar_forcing,get_outlines(geo)))


