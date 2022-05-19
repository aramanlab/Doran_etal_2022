using SparseArrays
using Combinatorics: combinations
using LinearAlgebra: Symmetric

# d = [
#     0 5  9  9 8
#     5 0 10 10 9
#     9 10 0  8 7
#     9 10 8  0 3
#     8  9 7  3 0
# ]

# fastNJ(d)
# @time regNJ(Symmetric(rand(1000,1000)))
# @time fastNJ(Symmetric(rand(10000,10000)))

function fastNJ(d::AbstractMatrix{<:Number}) 
    n = size(d, 1)
    n == size(d,2) || ArgumentError("d must be a square matrix")
    merges = Vector()
    heights = Vector()
    sd = spzeros(eltype(d), 2n-1, 2n-1)
    sd[1:n, 1:n] .= d
    currentindices = vcat(ones(Bool,n), zeros(Bool, n-1))
    fastNJ!(merges, heights, sd, currentindices)
    return vcat(merges...), vcat(heights...)
end

function fastNJ!(merges::AbstractVector, heights::AbstractVector, d::AbstractSparseMatrix{<:Number}, currentindices::AbstractVector{<:Bool})
    n = sum(currentindices)
    mergestep = 1
    while (sum(currentindices) > 1) && (mergestep < n)
        for (i, j) in _Qlist(d, currentindices)
            push!(merges, [mergeidx(i, n) mergeidx(j, n)])
            push!(heights, _distancetoparent(d, currentindices, i, j))
            
            currentindices[i] = false
            currentindices[j] = false

            for c in findall(currentindices)
                d[c, n+mergestep] = _distancetonewnode(d, i, j, c)
                d[n+mergestep, c] = d[c, n+mergestep]
            end
            currentindices[n+mergestep] = true
            mergestep += 1
        end
    end
end

# julia> regNJ(d)
# (Any[[-2 -1], [1 -3], [-5 -4], [3 2]], Any[[3.0 2.0], [3.0 4.0], [1.0 2.0], [1.0 1.0]])

function regNJ(d::AbstractMatrix{<:Number})
    n = size(d, 1)
    n == size(d,2) || ArgumentError("d must be a square matrix")
    merges = Vector()
    heights = Vector()
    sd = spzeros(eltype(d), 2n-1, 2n-1)
    sd[1:n, 1:n] .= d
    currentindices = vcat(ones(Bool,n), zeros(Bool, n-1))
    regNJ!(merges, heights, sd, currentindices)
    return vcat(merges...), vcat(heights...)
end


function regNJ!(merges::AbstractVector, heights::AbstractVector, d::AbstractSparseMatrix{<:Number}, currentindices::AbstractVector{<:Bool})
    n=sum(currentindices)
    for mergestep in 1:n-1
        tmpd = @view d[currentindices, currentindices]
        Qij = spzeros(eltype(d), size(d))
        Qij[currentindices, currentindices] .= _Q(tmpd)
        
        idxs = argmin(Qij)
        
        push!(merges, [mergeidx(idxs[1], n) mergeidx(idxs[2], n)])
        push!(heights, _distancetoparent(d, currentindices, idxs[1], idxs[2]))
        
        currentindices[idxs[1]] = false
        currentindices[idxs[2]] = false
        for c in findall(currentindices)
            d[c, n+mergestep] = _distancetonewnode(d, idxs[1], idxs[2], c)
            d[n+mergestep, c] = d[c, n+mergestep]
        end
        currentindices[n+mergestep] = true
    end
end

mergeidx(i, n) = i ≤ n ? -i : i-n

function _Q(d::AbstractMatrix{<:Number})
    n = size(d,1)
    Qij = spzeros(eltype(d), n,n)
    marginsums = mapslices(sum, d, dims=1)[:]
    for (j,i) in combinations(1:n, 2)
            Qij[j, i] = (n-2)*d[i,j] - marginsums[i] - marginsums[j]
    end
    return Qij
end

function _Qlist(d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool})
    indxs = (1:length(currentindices))[currentindices]
    n = sum(currentindices)
    Qk = zeros(eltype(d), 3, binomial(length(indxs), 2))
    marginsums = mapslices(sum, d, dims=1)[:]
    for (k, (j, i)) in enumerate(combinations(indxs, 2))
            Qk[1, k] = i
            Qk[2, k] = j
            Qk[3, k] = (n-2)*d[i,j] - marginsums[i] - marginsums[j]
    end
    Qk = Qk[:, sortperm(Qk[3,:])]
    _get_independent_merges(Qk)
end

function _get_independent_merges(q::AbstractMatrix{<:Number})
    s = Set{eltype(q)}()
    sl = length(s)
    indxs = Vector{Vector{Int64}}()
    for v in eachcol(q)
        push!(s, v[1])
        push!(s, v[2])
        length(s) == sl+2 || break
        sl = length(s)
        push!(indxs, v[1:2])
    end
    return indxs
end

function _distancetoparent(d::AbstractMatrix{<:Number}, currentindices::AbstractVector{<:Bool}, i::Integer, j::Integer)
    n = sum(currentindices)
    sum_i = sum(d[currentindices,i])
    sum_j = sum(d[currentindices,j])
    β = n > 2 ? (1/(2*(n-2)))*(sum_i - sum_j) : 0.0
    δa = (0.5) * d[i, j] + β
    δb = d[i, j] - δa
    return [δa  δb]
end

function _distancetonewnode(d::AbstractMatrix{<:Number}, a::Integer, b::Integer, c::Integer)
    return (1/2) * (d[a,c] + d[b,c] - d[a, b])
end

