"""
    getData()

Reads the CSV data file and return the dataframe, hosts and viruses.
"""
function getData()

    # Reading the data file
    df = CSV.read("./data/virionette.csv");
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
    crossValidation(interaction_matrix, targetedValue, initialValueMatrix, rank, suspected_newData, unlikely_newData)

Computes the imputation for a given matrix if the targeted value is 0
Computes the leave one out validation for a given matrix if the targeted value is 1
"""
function crossValidation(targetedValue, initialValueMatrix, α, rank, suspected_newData, unlikely_newData)
    #println("Targeted Value: $(targetedValue) ",
        #"Initial Value: $(round(initialValueMatrix, digits=2)) ", "Rank: $(rank)")
    println("Rank: $(rank)")
    println("Alpha: $(α)")

    (interaction_matrix, hosts, viruses, df) = buildInteractionMatrix()
    # Do the imputation for every targeted value
    positions_to_impute = findall(interaction_matrix .== targetedValue)
    output_matrix = copy(interaction_matrix)

    @showprogress for position in positions_to_impute
        output_matrix[position] = imputation(interaction_matrix,
            position, initialValueMatrix, rank)
    end
    # Visualizing the interactions
    display(plot(generateHeatmap("Initial Matrix \n(targeted value: $(targetedValue), rank: $(rank))", interaction_matrix),
        generateHeatmap("Output matrix", output_matrix,)))
    #println("Variation: $(calculateVariation(interaction_matrix, output_matrix))%")
    maxValues = getTopInteractions(10, interaction_matrix, output_matrix, hosts, viruses, df)
    generateResultsTable(maxValues, α, rank, suspected_newData, unlikely_newData)
end

"""
    imputation(matrix, position, initialValueMatrix, rank; tolerance=1e-2, maxiter=50)

This function returns the imputed value at a given position (expressed in
CartesianCoordinates) in the matrix, seeded from a value in initialeValueMatrix. The imputation
is done by iterating an SVD at a given rank, and stops when the iteration
difference is smaller than the absolute tolerance, or after maxiter steps
have been done.
"""
function imputation(matrix, position, initialValueMatrix, rank; tolerance = 1e-2, maxiter = 50)
    tempMatrix = copy(matrix)
    tempMatrix[position] = initialValueMatrix[position]
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
            # selecting only the betacoronaviruses and bat hosts
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
    # Printing the top10 interactions
    println("The top $(top) predicted interactions are:")
    for t in 1:top
        println("Virus: ", viruses[maxValues[t][2]], "  Host: ", hosts[maxValues[t][3]])
        #println("Virus: ", viruses[maxValues[t][2]], "  Host: ", hosts[maxValues[t][3]], " With: ", round(maxValues[t][1]*100, digits=1),"%")
    end
    println("The highest scoring interaction is:")
    println("Virus: ", viruses[maxValues[1][2]], "  Host: ", hosts[maxValues[1][3]], " With a Δ of: ", maxValues[1][1] - initial_matrix[maxValues[1][2],maxValues[1][3]])
    # Returning the array containing the top10
    return(maxValues)
end


"""
    calculateInitialeValues(Y, α)
Calculates the initial values to be attributed by linear filtering.

"""
function calculateInitialeValues(Y, α)
    # Convert ones and zeros in boolean
    matrix_bool = convert(Array{Bool}, Y.== 1)
    # Convert matrix in bipartite network
    B = BipartiteNetwork(matrix_bool)
    # Apply the linear filter
    F = linearfilter(B, α = α)
    #return the adjacency matrix
    return(F.A)
end


"""
    generateResultsTable(top10, α, rank, suspected_newData, unlikely_newData)

Generates an automated results table containing the combinations of alphas and
ranks, and the number of species predicted in the top10 that are also included
in the new dataset.
"""

function generateResultsTable(top10, α, rank, suspected_newData, unlikely_newData)
    suspected_count = 0;
    unlikely_count = 0;
    for i in 1:length(top10)
        if findfirst(suspected -> suspected == top10[i], suspected_newData) != nothing
            suspected_count += 1;
            println("Suspected count", suspected_count)
        end
        if findfirst(unlikely -> unlikely == top10[i], unlikely_newData) != nothing
            unlikely_count += 1;
            println("unlikely count", unlikely_count)
        end
    end
    results = DataFrame(Rank = [rank], Alpha = [α], Suspected = [suspected_count], Unlikely = [unlikely_count])
    CSV.write("Results.csv", results, append = true, delim =';')
end
