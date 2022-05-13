function lastline(io)
    local line
    for l in eachline(io)
        line = l
    end
    line
end

"""
        _numpairs2N(x)::Integer
    solve choose(n,k)=x for n
    for numbers around a trillion use a BigInt for x
"""
function _numpairs2N(x)::Integer
    # solve choose(n,k)=x for n
    n = (1 + sqrt(1 + 8x)) / 2
    # if n is whole, return n
    n % 1 == 0 ? n : throw(ArgumentError("Your vector has the wrong number of pairs to be pairwise combinations"))
    # n not whole means that x was not a binomial number
    return n
end

using Combinatorics: combinations
"""
    asdistancematrix(pairs::Vector; defaultval=zeros)
take the columnwise vectorized lower diagonal of distance matrix and remake a symetric distance matrix.
"""
function asdistancematrix(pairsvec; defaultval=zeros)
    n = _numpairs2N(length(pairsvec))
    distmtx = defaultval(n, n)
    for (k, (i,j)) in enumerate(combinations(1:n, 2))
        distmtx[i,j] = pairsvec[k]
        distmtx[j,i] = pairsvec[k]
    end
    return distmtx
end