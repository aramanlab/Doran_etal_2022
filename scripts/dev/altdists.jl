using DrWatson
@quickactivate "Doran_etal_2022"

using Distances, Clustering
using LinearAlgebra
using CSV, DataFrames, Muon, MAT
using SPI: nwstr
using Gotree_jll

const NBOOT = 100
const RANK = 10

## LOAD DATA ##
@info "Reading in data"
uniprot = readh5ad(joinpath(datadir(), "exp_pro", "UP7047", "2020_02_UP7047.h5ad"))
oggs = uniprot.X[:, :]
usv = SVD(uniprot.obsm["LSVs"][:,:], uniprot.uns["SVs"][:], uniprot.varm["RSVs"][:,:]');


## COMPUTE REFERENCE TREES ##
@info "computing reference trees" 
Deuclidean = pairwise(Euclidean(), oggs');
Dcityblock = pairwise(Cityblock(), oggs');
Dtop10euclidean = pairwise(Euclidean(), usv.U[:,1:10]');
Dtop10cityblock = pairwise(Cityblock(), usv.U[:,1:10]');

fastANIresults = matread(projectdir("_research", "FastANI", "fastANIdistmtx.mat"))
indxmapping = indexin(uniprot.obs_names, fastANIresults["ids"])
DfastANI = zeros(size(uniprot.obsp["SPI_distances"]))
DfastANI[findall(!isnothing, indxmapping), findall(!isnothing, indxmapping)] .= fastANIresults["mtx"][filter(!isnothing, indxmapping), filter(!isnothing, indxmapping)]
DfastANI = 100.0 .- DfastANI


uniprot.obsp["ogg_euclidean"] = Deuclidean
uniprot.obsp["ogg_cityblock"] = Dcityblock
uniprot.obsp["svd_euclidean"] = Dtop10euclidean
uniprot.obsp["svd_cityblock"] = Dtop10cityblock;
uniprot.obsp["fastANI_dist"] = DfastANI;

outdir = joinpath(projectdir(), "_research", "UP7047altdists")
mkpath(outdir)

## WRITING REFERENCE TREES ##
@info "writing reference trees..."
hc = hclust(Deuclidean, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false)
uniprot.uns["ogg_euclidean_reftreestring"] = nws
open(io->println(io, nws), joinpath(outdir, "UP7047_ogg_euclidean-tree.nw"), "w")

hc = hclust(Dcityblock, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false)
uniprot.uns["ogg_cityblock_reftreestring"] = nws
open(io->println(io, nws), joinpath(outdir, "UP7047_ogg_cityblock-tree.nw"), "w")

hc = hclust(Dtop10euclidean, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false)
uniprot.uns["svd_euclidean_reftreestring"] = nws
open(io->println(io, nws), joinpath(outdir, "UP7047_svd_euclidean-tree.nw"), "w")

hc = hclust(Dtop10cityblock, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false)
uniprot.uns["svd_cityblock_reftreestring"] = nws
open(io->println(io, nws), joinpath(outdir, "UP7047_svd_cityblock-tree.nw"), "w")

hc = hclust(DfastANI, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false)
uniprot.uns["svd_cityblock_reftreestring"] = nws
open(io->println(io, nws), joinpath(outdir, "UP7047_seq_fastANI-tree.nw"), "w")

@info "writing updates to h5ad file"
writeh5ad(joinpath(datadir(), "exp_pro", "UP7047", "2020_02_UP7047.h5ad"), uniprot)

# STARTING BOOTSTRAP ## 
@info "starting bootstrap for svd_euclidean..."
@time begin
nwss = Vector{String}()
Threads.@threads for i in 1:NBOOT
    m = size(oggs,2)
    rusv = svd(oggs[:, rand(1:m, m)])   
    Dij = pairwise(Euclidean(), rusv.U[:,1:RANK]')
    hc = hclust(Dij, linkage=:average, branchorder=:optimal)
    push!(nwss, nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false))
    print("\r$(i)/100 done")
end
println()
open(joinpath(outdir, "UP7047_svd_euclidean-boottrees.nw"), "w") do io
    for nws in nwss 
        println(io, nws)
    end
end
end#time

@info "starting bootstrap for svd_cityblock..."
@time begin
nwss = Vector{String}()
Threads.@threads for i in 1:NBOOT
    m = size(oggs,2)
    rusv = svd(oggs[:, rand(1:m, m)])
    Dij = pairwise(Cityblock(), rusv.U[:,1:RANK]')
    hc = hclust(Dij, linkage=:average, branchorder=:optimal)
    push!(nwss, nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false))
    print("\r$(i)/100 done")
end
println()
open(joinpath(outdir, "UP7047_svd_cityblock-boottrees.nw"), "w") do io
    for nws in nwss 
        println(io, nws)
    end
end
end#time

@info "starting bootstrap for ogg_euclidean..."
@time begin
nwss = Vector{String}()
Threads.@threads for i in 1:NBOOT
    m = size(oggs,2)
    Dij = pairwise(Euclidean(), oggs[:, rand(1:m, m)]')
    hc = hclust(Dij, linkage=:average, branchorder=:optimal)
    push!(nwss, nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false))
    print("\r$(i)/100 done")
end
println()
open(joinpath(outdir, "UP7047_ogg_euclidean-boottrees.nw"), "w") do io
    for nws in nwss 
        println(io, nws)
    end
end
end#time

@info "starting bootstrap for ogg_cityblock..."
@time begin
nwss = Vector{String}()
Threads.@threads for i in 1:NBOOT
    m = size(oggs,2)
    tmpmtx = copy(oggs[:, rand(1:m, m)]')
    Dij = pairwise(Cityblock(), tmpmtx)
    hc = hclust(Dij, linkage=:average, branchorder=:optimal)
    push!(nwss, nwstr(hc, uniprot.obs_names.vals; labelinternalnodes=false))
    print("\r$(length(nwss))/100 done  ")
end
println()
open(joinpath(outdir, "UP7047_ogg_cityblock-boottrees.nw"), "w") do io
    for nws in nwss 
        println(io, nws)
    end
end
end#time

## RUNNING BOOSTER ## 
@info "calculating TBE support values (ogg euclidean)..."
@time run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(outdir, "UP7047_ogg_euclidean-tree.nw")) \
            -b $(joinpath(outdir, "UP7047_ogg_euclidean-boottrees.nw")) \
            -o $(joinpath(outdir, "UP7047_ogg_euclidean-supporttree.nw"))`,
            stderr=joinpath(outdir, "booster.log")))
@info "calculating TBE support values (ogg cityblock)..."
@time run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(outdir, "UP7047_ogg_cityblock-tree.nw")) \
            -b $(joinpath(outdir, "UP7047_ogg_cityblock-boottrees.nw")) \
            -o $(joinpath(outdir, "UP7047_ogg_cityblock-supporttree.nw"))`,
            stderr=joinpath(outdir, "booster.log")))
@info "calculating TBE support values (svd euclidean)..."
@time run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(outdir, "UP7047_svd_euclidean-tree.nw")) \
            -b $(joinpath(outdir, "UP7047_svd_euclidean-boottrees.nw")) \
            -o $(joinpath(outdir, "UP7047_svd_euclidean-supporttree.nw"))`,
            stderr=joinpath(outdir, "booster.log")))
@info "calculating TBE support values (svd cityblock)..."
@time run(pipeline(`$(gotree()) compute support tbe --silent \
            -i $(joinpath(outdir, "UP7047_svd_cityblock-tree.nw")) \
            -b $(joinpath(outdir, "UP7047_svd_cityblock-boottrees.nw")) \
            -o $(joinpath(outdir, "UP7047_svd_cityblock-supporttree.nw"))`,
            stderr=joinpath(outdir, "booster.log")))

@info "Finished!"
