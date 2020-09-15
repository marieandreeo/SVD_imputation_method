# Creating a function to process the imputation
function imputation(data, rank, rows, cols)
    data[rows, cols] = 0.1
    for iterations in 1:1
        SVD = svd(data)
        SVD.S[rank:end] .= 0.0
        data = SVD.U*Diagonal(SVD.S)*transpose(SVD.V)
    end
    return data[rows, cols]
end
