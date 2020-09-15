using LinearAlgebra;
using Plots

data = [1.0 1.0 1.0 0.0 0.0 0.0; 1.0 1.0 1.0 0.0 0.0 0.0 ; 1.0 0.0 1.0 0.0 0.0 0.0 ; 0.0 0.0 0.0 1.0 1.0 1.0 ; 0.0 0.0 0.0 1.0 1.0 1.0 ; 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 1.0 1.0 1.0]
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
data[3,2] = new_data[3,2]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[3,2] = new_data[3,2]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[3,2] = new_data[3,2]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[3,2] = new_data[3,2]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
data[3,2] = new_data[3,2]
data
data_SVD = svd(data)
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)

plot(data_SVD.S)
data_SVD.S[3:end] .= 0
data_SVD.S
new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
heatmap(new_data)
