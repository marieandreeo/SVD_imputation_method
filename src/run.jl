include("./required.jl")
include("./functions.jl")

# Test the impuration
α = [0.0, 0.0, 0.0, 1.0];
(matrix, hosts, viruses, df) = buildInteractionMatrix();
crossValidation(0.0, calculateInitialeValues(matrix, α), 1);

# Test LOO
#crossValidation(1.0, 0.9, 10);
