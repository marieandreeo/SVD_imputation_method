include("./required.jl")
include("./functions.jl")

# Test the impuration with initial value = mean
matrix = buildInteractionMatrix();
#crossValidation(0.0, mean(matrix), 4);

# Test LOO
crossValidation(1.0, 0.9, 1);
