include("./required.jl")
include("./functions.jl")

# Testing initial value = mean
matrix = buildInteractionMatrix();
crossValidation(0.0, mean(matrix), 3)

# Test the imputation
#crossValidation(0.0, 0.035, 3)
# Test LOO
#crossValidation(1.0, 0.9, 3)
