using CSV
using DataFrames
using DataFramesMeta
using LightGraphs
using GraphPlot
using GraphDataFrameBridge
using EcologicalNetworks
using Distributions
using Plots
using LinearAlgebra


df = CSV.read("virionette.csv")
df.host_species
df.virus_genus
unique(df.host_species)
unique(df.virus_genus)

# Creating a Dictionary and Attributing host ID's
Dictionary_host = Dict(unique(df.host_species)[h] => h for h in 1:length(unique(df.host_species)))
# Test to make sure of the ID's attribution
Dictionary_host["Rhinolophus euryale"]

# Creating a Dictionary and Attributing virus ID's
Dictionary_virus = Dict(unique(df.virus_genus)[v-length(unique(df.host_species))] => v for v in length(unique(df.host_species))+1:length(unique(df.virus_genus))+length(unique(df.host_species)))
# Test to make sure of the ID's attribution
Dictionary_virus["Alphacoronavirus"]

##Create an interactions matrix
#create matrix of zeros
Matrix_interactions = zeros(Float64, length(unique(df.host_species))+1, length(unique(df.virus_genus))+1)
#assign names for rows
for r in 2:length(unique(df.host_species))+1
    Matrix_interactions[r,1] = Dictionary_host[unique(df.host_species)[r-1]]
end


#assign names for cols
for c in 2:length(unique(df.virus_genus))+1
    Matrix_interactions[1,c] = Dictionary_virus[unique(df.virus_genus)[c-1]]
end

# Fill matrix with ones where there are interactions
for interactions in 1:length(df.host_species)
    row_ID = Dictionary_host[df.host_species[interactions]]
    col_ID = Dictionary_virus[df.virus_genus[interactions]]
    Matrix_interactions[row_ID + 1, col_ID - length(unique(df.host_species)) + 1] = 1.0
end
heatmap(Matrix_interactions[2:end, 2:end])

df.host_species[1]
Dictionary_host[df.host_species[1]]
Dictionary_virus[df.virus_genus[1]]


Matrix_interactions_data = Matrix_interactions[2:end, 2:end]
new_Matrix = copy(Matrix_interactions_data)


function imputation(data, rank, rows, cols)
    data[rows, cols] == 0.1
    for iterations in 1:5
        SVD = svd(data)
        SVD.S[rank:end] .= 0.0
        data = SVD.U*Diagonal(SVD.S)*transpose(SVD.V)
    end
    return data[rows, cols]
end


for rows in 1:length(unique(df.host_species))
    for cols in 1:length(unique(df.virus_genus))
        if Matrix_interactions_data[rows, cols] == 0.0
            new_Matrix[rows, cols] = imputation(Matrix_interactions_data, 5, rows, cols)
        end
    end
end

heatmap(new_Matrix)
## SVD
Matrix_interactions_data_SVD = svd(Matrix_interactions_data)
plot(Matrix_interactions_data_SVD.S, seriestype = :scatter)
new_data = Matrix_interactions_data_SVD.U*Diagonal(Matrix_interactions_data_SVD.S)*transpose(Matrix_interactions_data_SVD.V)
epsilon = sum(Matrix_interactions_data .- new_data)

Matrix_interactions_data_SVD.S[5:end] .= 0.0
new_data = Matrix_interactions_data_SVD.U*Diagonal(Matrix_interactions_data_SVD.S)*transpose(Matrix_interactions_data_SVD.V)
heatmap(new_data)
Matrix_interactions_data[1,1] = new_data[1,1]
heatmap(Matrix_interactions_data)


function imputation(data, rank)
    data_SVD.S[rank:end] .= 0
    new_data = data_SVD.U*Diagonal(data_SVD.S)*transpose(data_SVD.V)
    #heatmap(new_data)
    return new_data[7,3]
end


# 1 technic to find top10
matrix = [10 2 4 10; 7 9 12 3; 14 2 5 7]
minValue = findmin(matrix)
minValue[2]
maxValues = zeros(5,3)
maxValues[1] = minValue[1]
maxValues[2] = minValue[2,1]
for w in eachindex(matrix)
    if matrix[w] > findmin(maxValues)[1]
        maxValues[findmin(maxValues)[2]] = matrix[w]
    end
end
