using ImageFiltering
using Plots

F = [1 1 0 0 0; 1 1 0 0 0; 0 0 1 1 1; 0 0 1 0 0]
heatmap(F, color=[:white, :blue])
imgg = imfilter(F, Kernel.gaussian(1))
heatmap(imgg, color=[:white, :blue])
imgl = imfilter(F, Kernel.Laplacian())
heatmap(imgl, color=[:white, :blue])
kern = centered([1/3, 1/3, 1/3])
im = imfilter(F, kern)
