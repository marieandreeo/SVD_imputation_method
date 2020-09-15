using CSV
using DataFrames
using Plots
using LinearAlgebra
include("functions.jl")

# Reading the data file
df = CSV.read("virionette.csv");
# Creating a Dictionary and Attributing host ID's
Dictionary_host = Dict(unique(df.host_species)[h] => h for h in 1:length(unique(df.host_species)));

# Creating a Dictionary and Attributing virus ID's
Dictionary_virus = Dict(unique(df.virus_genus)[v-length(unique(df.host_species))] => v for v in length(unique(df.host_species))+1:length(unique(df.virus_genus))+length(unique(df.host_species)));
## Creating an interactions matrix
# Creating a matrix of zeros
Matrix_interactions = zeros(Float64, length(unique(df.host_species))+1, length(unique(df.virus_genus))+1);

# Assigning IDs for rows
for r in 2:length(unique(df.host_species))+1
    Matrix_interactions[r,1] = Dictionary_host[unique(df.host_species)[r-1]]
end

# Assigning IDs for cols
for c in 2:length(unique(df.virus_genus))+1
    Matrix_interactions[1,c] = Dictionary_virus[unique(df.virus_genus)[c-1]]
end

Matrix_interactions

# Fill matrix with ones where there are interactions
for interactions in 1:length(df.host_species)
    row_ID = Dictionary_host[df.host_species[interactions]]
    col_ID = Dictionary_virus[df.virus_genus[interactions]]
    Matrix_interactions[row_ID + 1, col_ID - length(unique(df.host_species)) + 1] = 1.0
end

# Selecting only the data in the matrix
const Matrix_interactions_data = Matrix_interactions[2:end, 2:end];

# Visualizing the interactions
heatmap(Matrix_interactions_data , xlabel = "Virus", ylabel="Host")

# Creating a copy of the matrix data
new_Matrix = copy(Matrix_interactions_data);

# Iterating through each missing interactions for the imputation
for rows in 1:length(unique(df.host_species))
    for cols in 1:length(unique(df.virus_genus))
        if Matrix_interactions_data[rows, cols] == 0.0
            new_Matrix[rows, cols] = imputation(Matrix_interactions_data, 5, rows, cols)
        end
    end
end

# Visualizing the new interactions
heatmap(new_Matrix, xlabel = "Virus", ylabel="Host")
