"""
    lowrank(N, r)

Returns the low rank appproximation of the matrix N in rank r, after a SVD.
"""
function lowrank(N, r)
    @assert r ≤ rank(N)
    F = svd(N)
    F.S[(r + 1):end] .= zero(eltype(F.S))
    # Note that the SVD method returns Vt, which is the transposed right singular matrix
    return F.U * Diagonal(F.S) * F.Vt
end

"""
    imputation(N, position, i0; rank=3, tol=1e-2, maxiter=50)

This function returns the imputed value at a given position (expressed in
CartesianCoordinates) in the matrix N, seeded from a value i0. The imputation
is done by iterating an SVD at a given rank, and stops when the iteration
difference is smaller than the absolute tolerance, of after maxiter steps
have been done.
"""
function imputation(N, position, i0; rank = 3, tol = 1e-2, maxiter = 50)
    N[position] = i0
    Δ = 1.0
    iter = 1
    while Δ > tol
        iter += 1
        Z = lowrank(N, rank)
        # The change in value is the absolute difference between the next and current iteration
        Δ = abs(Z[position] - N[position])
        N[position] = Z[position]
        # We stop if there are more than a set number of iterations
        iter ≥ maxiter && break
    end
    return N[position]
end
