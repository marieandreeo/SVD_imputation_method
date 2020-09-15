using LinearAlgebra;
using Plots
using EcologicalNetworks


A = convert(BipartiteNetwork, web_of_life("A_HP_001"))



data = web_of_life("M_PL_001")
species(data)
interactions(data)
convert(Bool, data)
heatmap(data)
data_SVD = svd(data);
data_SVD.U
data_SVD.S
Diagonal(data_SVD.S)
data_SVD.V
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
epsilon = sum(data .- new_data)
plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[7,7] = new_data[7,7]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[7,7] = new_data[7,7]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[7,7] = new_data[7,7]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[7,7] = new_data[7,7]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[7,7] = new_data[7,7]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)




v = randn(5,5)
v[sortperm(v[:,1]),:]
heatmap(v)
