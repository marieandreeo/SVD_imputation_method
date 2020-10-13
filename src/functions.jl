"""
    getData()

Reads the CSV data file and return the dataframe, hosts and viruses.
"""
function getData()

    # Reading the data file
    df = CSV.read("./data/virionette.csv");
    #df=df[df.host_order .== "Chiroptera",:]
    #df = CSV.read("./data/TestData.csv");

    # Make a sorted list of unique hosts and viruses
    hosts = sort(unique(df.host_species))
    viruses = sort(unique(df.virus_genus))
    return df, hosts, viruses
end

"""
    buildInteractionMatrix()

Prepares the interaction matrix with initial data.
"""
function buildInteractionMatrix()

    (df, hosts, viruses) = getData()

    # Use this information to fill the matrix
    interaction_matrix = zeros(Float64, (length(viruses), length(hosts)))
    for interaction in eachrow(df)
        host_idx = findfirst(hosts .== interaction.host_species)
        virus_idx = findfirst(viruses .== interaction.virus_genus)
        interaction_matrix[virus_idx, host_idx] = 1.0
    end
    return interaction_matrix, hosts, viruses, df
end

"""
    crossValidation(interaction_matrix, targetedValue, initialValue, rank)

Computes the imputation for a given matrix if the targeted value is 0
Computes the leave one out validation for a given matrix if the targeted value is 1
"""
function crossValidation(targetedValue, initialValue, rank)
    println("Targeted Value: $(targetedValue) ",
        "Initial Value: $(round(initialValue, digits=2)) ", "Rank: $(rank)")

    (interaction_matrix, hosts, viruses, df) = buildInteractionMatrix()
    # Do the imputation for every targeted value
    positions_to_impute = findall(interaction_matrix .== targetedValue)
    output_matrix = copy(interaction_matrix)

    @showprogress for position in positions_to_impute
        output_matrix[position] = imputation(interaction_matrix,
            position, initialValue, rank)
    end
    # Visualizing the interactions
    display(plot(generateHeatmap("Initial Matrix \n(targeted value: $(targetedValue), rank: $(rank))", interaction_matrix),
        generateHeatmap("Output matrix", output_matrix,)))

    println("Variation: $(calculateVariation(interaction_matrix, output_matrix))%")
    getTopInteractions(10, interaction_matrix, output_matrix, hosts, viruses, df)
end

"""
    imputation(matrix, position, initialValue, rank; tolerance=1e-2, maxiter=50)

This function returns the imputed value at a given position (expressed in
CartesianCoordinates) in the matrix, seeded from a value initialeValue. The imputation
is done by iterating an SVD at a given rank, and stops when the iteration
difference is smaller than the absolute tolerance, or after maxiter steps
have been done.
"""
function imputation(matrix, position, initialValue, rank; tolerance = 1e-2, maxiter = 50)
    tempMatrix = copy(matrix)
    tempMatrix[position] = initialValue
    Δ = 1.0
    iter = 1
    while Δ > tolerance
        iter += 1
        approx_matrix = lowrank(tempMatrix, rank)
        # The change in value is the absolute difference between the
        # next and current iteration
        Δ = abs(approx_matrix[position] - tempMatrix[position])
        tempMatrix[position] = approx_matrix[position]
        # We stop if there are more than a set number of iterations
        iter ≥ maxiter && break
    end

    return tempMatrix[position]
end

"""
    lowrank(matrix, svd_rank)

Returns the low rank appproximation of the matrix, after a SVD.
"""
function lowrank(matrix, svd_rank)
    @assert svd_rank ≤ rank(matrix)
    factorization = svd(matrix)
    factorization.S[(svd_rank + 1):end] .= zero(eltype(factorization.S))
    # Note that the SVD method returns Vt, which is the transposed
    # right singular matrix
    return factorization.U * Diagonal(factorization.S) * factorization.Vt
end

"""
    generateHeatMap(title, matrix)

Returns a heatmap based on the given matrix.
"""
function generateHeatmap(title, matrix)

    return heatmap(matrix, xlabel = "Hosts", ylabel = "Viruses", c = :Greys,
        leg = true, frame = :box, title = title)

end

"""
    calculateVariation(initial_matrix, output_matrix)

Calculates the variation between two matrices.
"""
function calculateVariation(initial_matrix, output_matrix)

    delta = sum(abs.(initial_matrix .- output_matrix))
    return round((delta/length(initial_matrix)) * 100, digits = 2);

end

"""
    getTopInteractions(top, initial_matrix, output_matrix, hosts, viruses)

Returns the imputed interactions with the highest probability of occurrence.

"""
function getTopInteractions(top, initial_matrix, output_matrix, hosts, viruses, df)

    minValue = findmin(output_matrix)
    #Creating a buffer to stock the maximum probability values and their indexes
    maxValues = CircularBuffer{Tuple{Float64,Int64,Int64}}(top)
    # initializing the maximum value array to the minimum value
    push!(maxValues, (minValue[1],minValue[2][1],minValue[2][2]));
    for r in 1:size(output_matrix,1)
        for c in 1:size(output_matrix,2)
            # making sure that the present interaction was initially missing, and is now a max
            if initial_matrix[r,c] == 0 && occursin("Betacoronavirus", viruses[r]) &&  first(df[(df.host_species .== hosts[c]),:]).host_order == "Chiroptera" && output_matrix[r,c] > findmin(maxValues)[1][1]
            #if initial_matrix[r,c] == 0 && occursin("Betacoronavirus", viruses[r]) && output_matrix[r,c] > findmin(maxValues)[1][1]
                # sorting the maximums
                sort!(maxValues)
                # adding the new value at the end of tthe array and overwriting the smallest one
                push!(maxValues, (output_matrix[r,c],r,c))
            end
        end
    end
    # sorting the array in decreasing order
    sort!(maxValues, rev=true)
    println("The top $(top) predicted interactions are:")
    for t in 1:top
        println("Virus: ", viruses[maxValues[t][2]], "  Host: ", hosts[maxValues[t][3]], " With: ", round(maxValues[t][1]*100, digits=1),"%")
    end
end
