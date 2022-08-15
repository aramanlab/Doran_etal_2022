using DrWatson
@quickactivate "Doran_etal_2022"

using Muon
using SPI
using CSV, DataFrames
using NewickTreeTools
using CategoricalArrays: categorical
using Random: shuffle
using Statistics
using Distances
using StatsPlots
using Glob

NCUTS = 100
NPERMS = 5
TODAY = "2022-08-12"
taxaranklabels = [:Phylum, :Class, :Order, :Family, :Genus, :Species]
taxarankcolors = [:red :pink :orange :lightblue :green :aqua];

function clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)
    minmax = extrema(mapinternalnodes(distfun, tree, tree))
    cuts = range(0, minmax[2], length=ncuts)
    clusts = [cuttree(distfun, tree, cut) for cut in cuts]
    clustmappings = map(c->getleafnames.(c), clusts)
    clustersmps = [vcat(clustmapping...) for clustmapping in clustmappings]
    clusterids = [Int.(vcat([zeros(length(c)) .+ j for (j, c) in enumerate(clustmapping)]...)) for clustmapping in clustmappings];
    return clusterids, clustersmps
end

function pairedMIagainstmetacolumn(metacolumns, IDS, clusterids, clustersmps; doshuffle=false)
    tstat_MI = zeros(length(clusterids), size(metacolumns, 2))
    for (i, mcol) in enumerate(eachcol(metacolumns))
        # cat = levelorder(categorical(mcol))
        # pcat = cat .== cat'
        pcat = mcol .== permutedims(replace(mcol, ""=>"missing"))
        tstat_MI[:, i] .= collectMI_across_treedepth(clusterids, clustersmps, IDS, pcat; doshuffle)
    end
    DataFrame(tstat_MI, names(metacolumns)) |> stack |> df->rename!(df,["taxaID","MI"]);
end

function collectMI_across_treedepth(clusterids, clustersmps, IDS, ptax; doshuffle=false)
    uppertriangle = triu(trues(length(IDS), length(IDS)), 1);
    map(clusterids, clustersmps) do cids, smps
        clustorder = indexin(IDS,smps)
        pcids = cids[clustorder] .== cids[clustorder]'
        pcids = doshuffle ? shuffle(pcids[uppertriangle]) : pcids[uppertriangle]
        empiricalMI(ptax[uppertriangle], pcids)
    end
end


uniprot = readh5ad(joinpath(datadir(), "exp_pro", "UP7047", "2020_02_UP7047.h5ad"))
rowmeta = CSV.read(datadir("exp_raw", "UP7047", "UProwmeta.csv"), DataFrame);
UPtaxa = rowmeta[:, [:proteomeID, :Phylum, :Class, :Order, :Family, :Genus, :Species]];
UPtaxa = coalesce.(UPtaxa, "")
UPIDS = UPtaxa.proteomeID
# close(uniprot.file)

SPImtx = uniprot.obsp["SPI_distances"][:, :]


notclassmissing = UPtaxa[:, [:Phylum, :Class, :Order, :Family, :Genus, :Species]] |>
    df -> mapslices(x->!any(x.==""), Matrix(df), dims=2) |> vec
SPImtx[notclassmissing, notclassmissing]
hc = hclust(SPImtx[notclassmissing, notclassmissing], linkage=:average, branchorder=:optimal)
spi_tree = readnw(SPI.nwstr(hc, String.(UPIDS[notclassmissing]); labelinternalnodes=false))
string.(UPIDS[notclassmissing])

## Caclulate MI curves for SPI tree ##
# spi_tree = readnw(read(joinpath(projectdir(), "_research", "runSPIonUP7047rows", "2020_02_UP7047-supporttree.nw"), String));
as_polytomy!(spi_tree, fun=n->NewickTree.support(n)<0.5)
as_polytomy!(spi_tree, fun=n->NewickTree.distance(n)<1e-8)

@info "calculate tree cuts for SPI tree..."
@time clusterids, clustersmps = clusters_per_cutlevel(network_distance, spi_tree, NCUTS);

@info "calculate test MI for SPI tree..."
# @time tstatdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps)
@time tstatdf = pairedMIagainstmetacolumn(UPtaxa[notclassmissing,2:end], UPIDS[notclassmissing], clusterids, clustersmps)

tstatdf

plot(title="SPI tree", ylabel="MI", xlabel="Tree depth", legend=:topleft)
minmax = extrema(mapinternalnodes(network_distance, spi_tree, spi_tree))
cuts = range(0, minmax[2], length=NCUTS)

for (tlab, tcol) in zip(string.(taxaranklabels), taxarankcolors)
    df = filter(:taxaID=> ==(tlab), tstatdf)
    @df df plot!(cuts, scaledcumsum(:MI), label=tlab, c=tcol, lw=3)
end
plot!()

@info "calculate perm MI for SPI tree..."
for i in 1:NPERMS
    print("\rworking on permutation $i   ")
    tmpdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps; doshuffle=true)
    tstatdf[!, "MI_perm$i"] = tmpdf.MI
end
@info "write csv results SPI tree..."
CSV.write(joinpath(datadir(), "exp_pro", "UP7047", "$(TODAY)_MI-spitree_treedepth-by-taxa.csv"), tstatdf)

## Caclulate MI curves for altdist trees ##
@info "starting on Alt Dist trees"
cd(projectdir())
altdisttreefiles = glob(joinpath("_research", "UP7047altdists", "*-tree.nw"))
altdistnames = replace.(first.(split.(basename.(altdisttreefiles), "-")), "_"=>"-")
altdisttrees = readnw.(read.(altdisttreefiles, String));
as_polytomy!.(altdisttrees, fun=n->NewickTree.support(n)<0.5)
as_polytomy!.(altdisttrees, fun=n->NewickTree.distance(n)<1e-8)

@time begin
for (nm, tree) in zip(altdistnames, altdisttrees)
    @info "calculate tree cuts for $nm..."
    @time clusterids, clustersmps = clusters_per_cutlevel(network_distance, tree, NCUTS);

    @info "calculate test MI for $nm..."
    @time tstatdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps);

    # @info "calculate perm MI for $nm..."
    # for i in 1:NPERMS
    #     print("\rworking on permutation $i   ")
    #     tmpdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps; doshuffle=true)
    #     tstatdf[!, "MI_perm$i"] = tmpdf.MI
    # end
    @info "writing csv results for $nm..."
    CSV.write(joinpath(datadir(), "exp_pro", "UP7047", "$(TODAY)_MI-$(nm)_treedepth-by-taxa.csv"), tstatdf)
end
end # time

## Calculate Max depth of trees

maxdepthdf = DataFrame()
minmax = extrema(mapinternalnodes(network_distance, spi_tree, spi_tree))
maxdepthdf = vcat(maxdepthdf, DataFrame("treename"=>"SPI", "maxdepth"=>minmax[2]))

for (nm, tree) in zip(altdistnames, altdisttrees)
    minmax = extrema(mapinternalnodes(network_distance, tree, tree))
    maxdepthdf = vcat(maxdepthdf, DataFrame("treename"=>nm, "maxdepth"=>minmax[2]))
end

maxdepthdf

CSV.write(datadir("exp_pro", "UP7047", "$(TODAY)_treedepths.csv"), maxdepthdf)

