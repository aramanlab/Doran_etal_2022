using DataFrames

"""
    readphylip(fn::String)

Read phylip alignment file, return dataframe of IDs and Sequences
"""
function readphylip(fn::String)
    smps, seqs = open(fn, "r") do reader
        nsmp, nfeat = parse.(Int, split(readline(reader)))
        smps = Vector()
        seqs = Vector()
        for l in 1:nsmp
            smp, seq = split(readline(reader))
            push!(smps, smp)
            push!(seqs, seq)
        end
        all(s->length(s)==nfeat, seqs) || throw(DimensionMismatch("Length of all sequences must be the equal"))

        ord = sortperm(smps)
        return string.(smps[ord]), seqs[ord]
    end
    return DataFrame("ids" => string.(smps), "seqs" => string.(seqs))
end # read phylip



onehotencode(seqs::Vector{<:AbstractString}) = onehotencode(_stringcolumntochardf(seqs))
function onehotencode(df::D) where D<:AbstractDataFrame
    ohedf = DataFrame()
    origsize=size(df,2)
    for i in 1:origsize
        ohedf = hcat(ohedf, onehotencode(df, string(i)))
    end
    ohedf
end 
function onehotencode(df::D, col, cate = sort(unique(df[!, col])); outnames = Symbol.(:ohe_,col,cate)) where D<:AbstractDataFrame
    select(df, @. col => ByRow(isequal(cate)) .=> outnames)
end

function _stringcolumntochardf(seqs)
    DataFrame(reduce(hcat, collect.(seqs)) |> permutedims, [string(i) for i in 1:length(first(seqs))])
end