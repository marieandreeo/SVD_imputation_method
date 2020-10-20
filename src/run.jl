include("./required.jl")
include("./functions.jl")

# Test the imputation
α = [0.0, 0.0, 1.0, 0.0];
(matrix, hosts, viruses, df) = buildInteractionMatrix();
crossValidation(0.0, calculateInitialeValues(matrix, α), 5);

# Test LOO
#crossValidation(1.0, 0.9, 10);
