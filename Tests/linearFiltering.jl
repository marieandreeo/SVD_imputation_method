using ImageFiltering
using Plots

# Testing linear filtering of images
Y = [1.0 1.0 0.0 0.0 0.0; 1.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 1.0 0.0 0.0]
heatmap(Y, color=[:white, :blue])
imgg = imfilter(Y, Kernel.gaussian(1))
heatmap(imgg, color=[:white, :blue])
imgl = imfilter(Y, Kernel.Laplacian())
heatmap(imgl, color=[:white, :blue])
kern = centered([1/3, 1/3, 1/3])
im = imfilter(Y, kern)
Y

# Reproducing Stock et al. example
F = zeros(size(Y))
a = [1/4, 1/4, 1/4, 1/4];
for i in 1:size(Y,1)
    for j in 1:size(Y,2)
        F[i,j] = (a[1]*Y[i,j]) + (a[2]/size(Y,1)*sum(Y[:,j])) + (a[3]/size(Y,2)*sum(Y[i,:])) + (a[4]/length(Y)*sum(Y))
    end
end
heatmap(F, color=[:white, :blue])
F
