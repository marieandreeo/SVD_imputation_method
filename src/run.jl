include("./src/required.jl")
include("./src/functions.jl")


Y = [1.0 1.0 0.0 0.0 0.0; 1.0 1.0 0.0 0.0 0.0; 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 1.0 0.0 0.0]

F = linearFilter([1/4 1/4 1/4 1/4], Y)
heatmap(F, color=[:white, :blue])

# Test the impuration with initial value = mean
(matrix, hosts, viruses, df) = buildInteractionMatrix();
crossValidation(0.0, mean(matrix), 2);

# Test LOO
#crossValidation(1.0, 0.9, 10);
