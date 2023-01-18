using SparseArrays

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

function match_column_order(mtx::Matrix{<:Number}, cnames_src, cnames_dst)
    # make destination matrix
    dstmtx = zeros(eltype(mtx), size(mtx, 1), length(cnames_dst))

    # find matching columns
    rawidxs = indexin(cnames_src, cnames_dst)
    mask = .!isnothing.(rawidxs)
    matchedcols = filter(x->.!isnothing.(x), rawidxs)
    
    # set values where we have matched columns
    dstmtx[:, matchedcols] .= mtx[:, mask]
    return dstmtx
end

onehotencode(seqs::AbstractVector{<:AbstractString}) = onehotencode(_stringcolumntocharmtx(seqs))
function onehotencode(chardf::AbstractMatrix{<:Char})
    ohemtx = Vector()
    for col in eachcol(chardf)
        push!(ohemtx, indicatormat(col)')
    end
    return sparse(hcat(ohemtx...))
end

function _stringcolumntocharmtx(seqs)
    Matrix(DataFrame(reduce(hcat, collect.(seqs)) |> permutedims, [string(i) for i in 1:length(first(seqs))]))
end

function writefasta(path, ids, seqs)
    FASTA.Writer(open(path, "w")) do io
        for (id, seq) in zip(ids, seqs) 
            write(io, FASTA.Record(id, seq))
        end
    end
end

function readfasta(path) 
    FASTA.Reader(open(path)) do io
        recs = collect(io)
        DataFrame(
            :label => String.(identifier.(recs)), 
            :sample_id => String.(description.(recs)), 
            :sequence => String.(sequence.(recs)),
        )
    end
end