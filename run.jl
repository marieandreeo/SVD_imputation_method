using CSV
using DataFrames
using Plots
using LinearAlgebra
using Statistics
using ProgressMeter

include("functions.jl")

# Reading the data file
df = CSV.read("virionette.csv");

# Make a sorted list of unique hosts and viruses
const hosts = sort(unique(df.host_species))
const viruses = sort(unique(df.virus_genus))

# Use this information to fill the matrix
M = zeros(Float64, (length(viruses), length(hosts)))
for interaction in eachrow(df)
    host_idx = findfirst(hosts .== interaction.host_species)
    virus_idx = findfirst(viruses .== interaction.virus_genus)
    M[virus_idx, host_idx] = 1.0
end

# Visualizing the interactions
heatmap(M, xlabel = "Hosts", ylabel = "Viruses", c = :Greys, leg = false, frame = :box)

# Do the imputation for every zero value
positions_to_impute = findall(M .== 0.0)
K = copy(M)

@showprogress for p in positions_to_impute
    K[p] = imputation(M, p, mean(M))
end

# Visualizing the interactions
heatmap(K, xlabel = "Hosts", ylabel = "Viruses", c = :Greys, leg = false, frame = :box)
