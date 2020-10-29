include("./src/required.jl")
include("./src/functions.jl")

# Test the imputation
header_list = DataFrame(Rank = [], Alpha = [], Suspected = [], Unlikely = [], Delta = []);
CSV.write("Results.csv", header_list, header = true, delim =';');
α_list = [[0.0, 0.0, 0.0, 1.0], [0.0, 1.0, 1.0, 0.0], [0.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]];
rank_list = [1, 2, 3, 4, 5];
suspected_newData = ["Artibeus jamaicensis","Desmodus rotundus", "Hipposideros larvatus","Hipposideros pomona", "Megaerops kusnotoi", "Myonycteris angolensis", "Myotis pequinius", "Nanonycteris veldkampii", "Nycteris macrotis", "Pipistrellus deserti", "Plecotus auritus", "Pteropus lylei", "Scotophilus heathii", "Scotophilus kuhlii"];
unlikely_newData = ["Carollia sowelli", "Hipposideros gigas", "Hipposideros lekaguli", "Macroglossus minimus", "Myotis horsfieldii", "Pipistrellus coromandra", "Tadarida teniotis"];

for i in 1:length(α_list)[1]
    for j in 1:length(rank_list)
        (matrix, hosts, viruses, df) = buildInteractionMatrix();
        crossValidation(0.0, calculateInitialeValues(matrix, α_list[i]), α_list[i], rank_list[j], suspected_newData, unlikely_newData);
    end
end

# Test LOO
#crossValidation(1.0, 0.9, 10);
