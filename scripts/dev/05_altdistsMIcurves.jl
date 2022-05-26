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
using Glob

const NCUTS = 100
const NPERMS = 5
taxaranklabels = [:Phylum, :Class, :Order, :Family, :Genus, :Species]


function clusters_per_cutlevel(distfun::Function, tree::Node, ncuts::Number)
    minmax = extrema(mapinternalnodes(distfun, tree, tree))
    cuts = range(0, minmax[2], length=ncuts)
    clusts = [cuttree(distfun, tree, cut) for cut in cuts]
    clustmappings = map(c->getleafnames.(c), clusts)
    clustersmps = [vcat(clustmapping...) for clustmapping in clustmappings]
    clusterids = [Int.(vcat([zeros(length(c)) .+ j for (j, c) in enumerate(clustmapping)]...)) for clustmapping in clustmappings];
    return clusterids, clustersmps
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

function pairedMIagainstmetacolumn(metacolumns, IDS, clusterids, clustersmps; doshuffle=false)
    tstat_MI = zeros(length(clusterids), size(metacolumns, 2))
    for (i, mcol) in enumerate(eachcol(metacolumns))
        cat = levelorder(categorical(mcol))
        pcat = cat .== cat'
        tstat_MI[:, i] .= collectMI_across_treedepth(clusterids, clustersmps, IDS, pcat; doshuffle)
    end
    DataFrame(tstat_MI, names(metacolumns)) |> stack |> df->rename!(df,["taxaID","MI"]);
end


uniprot = readh5ad(joinpath(datadir(), "exp_pro", "UP7047", "2022-02-22_UP7047.h5ad"))
UPtaxa = uniprot.obs[:, [:proteomeID, :Phylum, :Class, :Order, :Family, :Genus, :Species]];
UPIDS = uniprot.obs_names.vals
close(uniprot.file)

## Caclulate MI curves for SPI tree ##
spi_tree = readnw(read(joinpath(projectdir(), "_research", "runSPIonUP7047rows", "2022-02-22_UP7047-supporttree.nw"), String));
as_polytomy!(spi_tree, fun=n->NewickTree.support(n)<0.5)
as_polytomy!(spi_tree, fun=n->NewickTree.distance(n)<1e-8)

@info "calculate tree cuts for SPI tree..."
@time clusterids, clustersmps = clusters_per_cutlevel(network_distance, spi_tree, NCUTS);

@info "calculate test MI for SPI tree..."
@time tstatdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps)

@info "calculate perm MI for SPI tree..."
for i in 1:NPERMS
    print("\rworking on permutation $i   ")
    tmpdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps; doshuffle=true)
    tstatdf[!, "MI_perm$i"] = tmpdf.MI
end
@info "write csv results SPI tree..."
CSV.write(joinpath(datadir(), "exp_pro", "UP7047", "2022-05-23_MI-spitree_treedepth-by-taxa.csv"), tstatdf)

## Caclulate MI curves for altdist trees ##
@info "starting on Alt Dist trees"
cd(projectdir())
altdisttreefiles = glob(joinpath("_research", "UP7047altdists", "*-supporttree.nw"))
altdistnames = replace.(first.(split.(basename.(altdisttreefiles), "-")), "_"=>"-")
altdisttrees = readnw.(read.(altdisttreefiles, String));
as_polytomy!.(altdisttrees, fun=n->NewickTree.support(n)<0.5)
as_polytomy!.(altdisttrees, fun=n->NewickTree.distance(n)<1e-8)


for (nm, tree) in zip(altdistnames, altdisttrees)
    @info "calculate tree cuts for $nm..."
    @time clusterids, clustersmps = clusters_per_cutlevel(network_distance, tree, NCUTS);

    @info "calculate test MI for $nm..."
    @time tstatdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps);

    @info "calculate perm MI for $nm..."
    for i in 1:NPERMS
        print("\rworking on permutation $i   ")
        tmpdf = pairedMIagainstmetacolumn(UPtaxa[!,2:end], UPIDS, clusterids, clustersmps; doshuffle=true)
        tstatdf[!, "MI_perm$i"] = tmpdf.MI
    end
    @info "writing csv results for $nm..."
    CSV.write(joinpath(datadir(), "exp_pro", "UP7047", "2022-05-23_MI-$(nm)_treedepth-by-taxa.csv"), tstatdf)
end




