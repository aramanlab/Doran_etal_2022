using DrWatson
@quickactivate "Doran_etal_2022"

using CSV, DataFrames
using SparseArrays
using Muon, MAT
using Combinatorics: combinations
using Statistics

researchdir() = joinpath(projectdir(), "_research")

colnames = ["query", "reference", "ANI", "matched_orthologs", "total_sequences"]
fastanidf = CSV.read(joinpath(researchdir(), "FastANI", "fullFastANI.tsv"), DataFrame; delim="\t", header=colnames)
fastanidf[!, :query_id] = first.(split.(basename.(fastanidf.query), "_"))
fastanidf[!, :reference_id] = first.(split.(basename.(fastanidf.reference), "_"))
fastanidf = sort(fastanidf, [:reference_id, :query_id])
allANIids = union(unique(fastanidf.query_id), unique(fastanidf.reference_id))


ANIdistdf = unstack(fastanidf, :reference_id, :query_id, :ANI; fill=0.0, allowduplicates=true)

ANImtx = Matrix(ANIdistdf[:, 2:end])
ANImtx = ANImtx[sortperm(ANIdistdf.reference_id), sortperm(names(ANIdistdf)[2:end])]

for (j, i) in combinations(1:length(allANIids), 2)
    m = max(ANImtx[j, i], ANImtx[i, j])
    ANImtx[j, i] = m
    ANImtx[i, j] = m
end
ANImtx - ANImtx' ≈ zeros(8852, 8852)

all(sort(ANIdistdf.reference_id) .== sort(names(ANIdistdf)[2:end]) .== sort(allANIids)) 

matwrite(joinpath(researchdir(), "FastANI", "fastANIdistmtx.mat"), Dict(
    "mtx"=>ANImtx,
    "ids"=>sort(allANIids),
), compress=true)