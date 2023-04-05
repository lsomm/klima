using Plots
using DelimitedFiles
pythonplot()

function read_geography(filepath)
     geo_matrix = readdlm(filepath)
end

read_geography("./input/The_World128x65.dat")

# 1 - land
# 2 - sea ice
# 3 - snow
# 5 - ocean



function robinson_projection(nlatitude, nlongitude)
     # Example data to demonstrate plotting against points that don't lie on a rectangular grid.
     x = [1.0 2.0
          1.0 1.8]

     y = [1.0 1.0
          2.0 1.8]

     return x, y
end

function plot_geo(geo_dat)
     # Example to demonstrate plotting against points that don't lie on a rectangular grid.
     # This does not work with other backends. We need PythonPlot for that.
     x, y = robinson_projection(2, 2)

     z = [1.0 2.0
          3.0 4.0]

     plot = contourf(x, y, z, aspect_ratio=1)

     return plot
end

# Run code
plot_geo(nothing)
