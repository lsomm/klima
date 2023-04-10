using Plots
using DelimitedFiles
pythonplot()

function read_geography(filepath)
    geo_matrix = readdlm(filepath)
    return geo_matrix'
end

geo_matrix = read_geography("./input/The_World128x65.dat")

# 1 - land
# 2 - sea ice
# 3 - snow
# 5 - ocean

function robinson_projection(nlatitude, nlongitude)
    #i = 1:nlongitude
    #j = 1:nlatitude
    longitude(i) = -pi + (i - 1) * 2.8125 * pi / 180
    latitude(j) = pi / 2 - (j - 1) * 2.8125 * pi / 180
    x(phi, theta) = phi / pi * (0.0379 * theta^6 - 0.15 * theta^4 - 0.367 * theta^2 + 2.666)
    y(phi, theta) = 0.96947 * theta - 0.00857 * sign(theta) * abs(theta)^6.41
    X = [x(longitude(i), latitude(j)) for i in 1:nlongitude, j in 1:nlatitude]
    Y = [y(longitude(i), latitude(j)) for i in 1:nlongitude, j in 1:nlatitude]
    return X, Y
end

function plot_geo(geo_dat)
    # Example to demonstrate plotting against points that don't lie on a rectangular grid.
    # This does not work with other backends. We need PythonPlot for that.
    x, y = robinson_projection(65, 128)
    print(x)
    z = [1.0 2.0
         3.0 4.0]

    plot = contourf(x, y, geo_matrix, aspect_ratio=1)

    return plot
end

# Run code
plot_geo(nothing)
