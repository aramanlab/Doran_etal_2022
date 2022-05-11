using DrWatson
@quickactivate "Doran_etal_2022"
using NewickTreeTools
using CSV, DataFrames, Glob
using Statistics
include(joinpath(srcdir(), "helpers.jl"))

const MAXSECS = 36 * 60 * 60 # 36 hours = max time jobs were allowed to run 
const SUPPORT_THRESHOLD = .90

metricdf = CSV.read(joinpath(datadir(), "exp_pro", "MSAs", "MSAs-metrics.csv"), DataFrame)
metricdf[!, :msaname] .= first.(split.(basename.(metricdf.msafile), "."))
rundir = pwd()
cd(projectdir())

#=
## Calculate fscore, precision, & recall of runSPI

Recall is calculated as what fraction of bipartitions in the source tree are 
present in the predicted tree. Precision is what fraction of bipartitions
in the predicted tree are present in the original tree. F-score combines those
fractions as (1 + β) * (precision * recall) / (β * precision + recall) where β
of one equally weights precision and recall, β less than one weights recall more
than precision and β greater than one weights precision more highly. For this 
calculation we hold β at 1.0. 
=#
### SPI METRICS
@info "Finding SPI trees..."
predtreesfiles = glob(joinpath("_research", "runSPI", "**", "*-supporttree.nw"))
origtreesfiles = joinpath.("data", "sims", "trees", replace.(basename.(predtreesfiles), r"-l[0-9].*"=>".nw"))
predtreesmsaname = replace.(basename.(predtreesfiles), r"-supporttree.nw"=>"")
all(metricdf.msaname .== predtreesmsaname) || ErrorException("origtrees is incorrectly ordered")

@info "Collecting SPI elapsed times..."
timefiles = joinpath.("_research", "runSPI", predtreesmsaname, "runSPI.out");
elapsedsecs = parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>""))
metricdf[!,:SPI_elapsedsecs] .= elapsedsecs

@info "Reading original source trees for SPI..."
origtrees = readnw.(read.(origtreesfiles, String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

@info "Reading in SPI predicted trees..."
predtrees = readnw.(read.(predtreesfiles, String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n)<SUPPORT_THRESHOLD);
@info "Calculate SPI Fscore Precision & Recall..."
proformencemetrics = fscore_precision_recall.(origtrees, predtrees)
proformencemetrics = hcat(vcat.(proformencemetrics...)...)
metricdf[!,:SPI_fscore] .= proformencemetrics[:, 1]
metricdf[!,:SPI_precision] .= proformencemetrics[:, 2]
metricdf[!,:SPI_recall] .= proformencemetrics[:, 3]

# sort(getleafnames(origtrees[1])) .== sort(getleafnames(predtrees[1]))
# keys(tally_tree_bifurcations(origtrees[1]))

@info "Calculating SPI Branch Depth"
# recalledtrees = deepcopy(predtrees)
# as_polytomy!(recalledtrees, fun=)
mmmdepths = map(predtrees) do tr
    dists = mapinternalnodes(NewickTreeTools.network_distance, tr, tr)
    return [mean(dists), median(dists), maximum(dists)]
end
mmmdepths = hcat(mmmdepths...)'
metricdf[!,:SPI_meandepth] .= mmmdepths[:, 1]
metricdf[!,:SPI_mediandepth] .= mmmdepths[:, 2]
metricdf[!,:SPI_maxdepth] .= mmmdepths[:, 3]

### FASTTREE METRICS
@info "Finding FastTree trees..."
predtreesfiles = glob(joinpath("_research", "runFastTree", "**", "*-supporttree.nw"))
origtreesfiles = joinpath.("data", "sims", "trees", replace.(basename.(predtreesfiles), r"-l[0-9].*"=>".nw"))
predtreesmsaname = replace.(basename.(predtreesfiles), r"-supporttree.nw"=>"")
all(metricdf.msaname .== predtreesmsaname) || ErrorException("origtrees is incorrectly ordered")

@info "Collecting FastTree elapsed times..."
timefiles = joinpath.("_research", "runFastTree", predtreesmsaname, "runFastTree.out");
elapsedsecs = parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>""))
metricdf[!,:FastTree_elapsedsecs] .= elapsedsecs

@info "Reading original source trees for FastTree..."
origtrees = readnw.(read.(origtreesfiles, String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

@info "Reading in FastTree predicted trees..."
predtrees = readnw.(read.(predtreesfiles, String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n)<SUPPORT_THRESHOLD);
@info "Calculate FastTree Fscore Precision & Recall..."
proformencemetrics = fscore_precision_recall.(origtrees, predtrees)
proformencemetrics = hcat(vcat.(proformencemetrics...)...)
metricdf[!,:FastTree_fscore] .= proformencemetrics[:, 1]
metricdf[!,:FastTree_precision] .= proformencemetrics[:, 2]
metricdf[!,:FastTree_recall] .= proformencemetrics[:, 3]

@info "Calculating FastTree Branch Depth"
mmmdepths = map(predtrees) do tr
    dists = mapinternalnodes(NewickTreeTools.network_distance, tr, tr)
    return [mean(dists), median(dists), maximum(dists)]
end
mmmdepths = hcat(mmmdepths...)'
metricdf[!,:FastTree_meandepth] .= mmmdepths[:, 1]
metricdf[!,:FastTree_mediandepth] .= mmmdepths[:, 2]
metricdf[!,:FastTree_maxdepth] .= mmmdepths[:, 3]

### FASTME METRICS
@info "Finding FastME trees..."
predtreesfiles = glob(joinpath("_research", "runFastME", "**", "*-supporttree.nw"))
origtreesfiles = joinpath.("data", "sims", "trees", replace.(basename.(predtreesfiles), r"-l[0-9].*"=>".nw"))
predtreesmsaname = replace.(basename.(predtreesfiles), r"-supporttree.nw"=>"")
all(metricdf.msaname .== predtreesmsaname) || ErrorException("origtrees is incorrectly ordered")

@info "Collecting FastME elapsed times..."
timefiles = joinpath.("_research", "runFastME", predtreesmsaname, "runFastME.out");
elapsedsecs = parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>""))
metricdf[!,:FastME_elapsedsecs] .= elapsedsecs

@info "Reading original source trees for FastME..."
origtrees = readnw.(read.(origtreesfiles, String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

@info "Reading in FastME predicted trees..."
predtrees = readnw.(read.(predtreesfiles, String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n)<SUPPORT_THRESHOLD);
@info "Calculate FastME Fscore Precision & Recall..."
proformencemetrics = fscore_precision_recall.(origtrees, predtrees)
proformencemetrics = hcat(vcat.(proformencemetrics...)...)
metricdf[!,:FastME_fscore] .= proformencemetrics[:, 1]
metricdf[!,:FastME_precision] .= proformencemetrics[:, 2]
metricdf[!,:FastME_recall] .= proformencemetrics[:, 3]

@info "Calculating FastME Branch Depth"
mmmdepths = map(predtrees) do tr
    dists = mapinternalnodes(NewickTreeTools.network_distance, tr, tr)
    return [mean(dists), median(dists), maximum(dists)]
end
mmmdepths = hcat(mmmdepths...)'
metricdf[!,:FastME_meandepth] .= mmmdepths[:, 1]
metricdf[!,:FastME_mediandepth] .= mmmdepths[:, 2]
metricdf[!,:FastME_maxdepth] .= mmmdepths[:, 3]

### RAxML METRICS
@info "Finding RAxML trees..."
predtreesfolders = glob(joinpath("_research", "runRAxML", "*"))
predtreesmsaname = basename.(predtreesfolders)
predtreesfiles = joinpath.("_research", "runRAxML", predtreesmsaname, "RAxML-supporttree.nw")
origtreesfiles = joinpath.("data", "sims", "trees", replace.(predtreesmsaname, r"-l[0-9].*"=>".nw"))
predictionexists = isfile.(predtreesfiles)
all(metricdf.msaname .== predtreesmsaname) || ErrorException("origtrees is incorrectly ordered")

@info "Reading original source trees for RAxML..."
origtrees = readnw.(read.(origtreesfiles[predictionexists], String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

@info "Reading in RAxML predicted trees..."
predtrees = readnw.(read.(predtreesfiles[predictionexists], String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n) < (SUPPORT_THRESHOLD * 100));
@info "Calculate RAxML Fscore Precision & Recall..."
proformencemetrics = fscore_precision_recall.(origtrees, predtrees)
proformencemetrics = hcat(vcat.(proformencemetrics...)...)
tmpmtx = Array{Union{Missing, Float64}}(missing, length(predtreesfiles), 3)
tmpmtx[predictionexists, :] .= proformencemetrics
metricdf[!,:RAxML_fscore] .= tmpmtx[:, 1]
metricdf[!,:RAxML_precision] .= tmpmtx[:, 2]
metricdf[!,:RAxML_recall] .= tmpmtx[:, 3]

@info "Calculating RAxML Branch Depth"
mmmdepths = map(predtrees) do tr
    dists = mapinternalnodes(NewickTreeTools.network_distance, tr, tr)
    return [mean(dists), median(dists), maximum(dists)]
end
mmmdepths = hcat(mmmdepths...)'
tmpmtx = Array{Union{Missing, Float64}}(missing, length(predtreesfiles), 3)
tmpmtx[predictionexists, :] .= mmmdepths
metricdf[!,:RAxML_meandepth] .= tmpmtx[:, 1]
metricdf[!,:RAxML_mediandepth] .= tmpmtx[:, 2]
metricdf[!,:RAxML_maxdepth] .= tmpmtx[:, 3]

@info "Collecting RAxML elapsed times..."
timefiles = joinpath.("_research", "runRAxML", predtreesmsaname, "runRAxML.out");
elapsedsecs =  
    parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>"", r"\[.*"=>"$MAXSECS")) 
metricdf[!,:RAxML_elapsedsecs] .= elapsedsecs

### PhyML METRICS
@info "Finding PhyML trees..."
predtreesfolders = glob(joinpath("_research", "runPhyML", "*"))
predtreesmsaname = basename.(predtreesfolders)
predtreesfiles = joinpath.("_research", "runPhyML", predtreesmsaname, predtreesmsaname .* ".phy-supporttree.txt")
origtreesfiles = joinpath.("data", "sims", "trees", replace.(predtreesmsaname, r"-l[0-9].*"=>".nw"))
predictionexists = isfile.(predtreesfiles)
all(metricdf.msaname .== predtreesmsaname) || ErrorException("origtrees is incorrectly ordered")

@info "Reading original source trees for PhyML..."
origtrees = readnw.(read.(origtreesfiles[predictionexists], String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

@info "Reading in PhyML predicted trees..."
predtrees = readnw.(read.(predtreesfiles[predictionexists], String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n)<SUPPORT_THRESHOLD);
@info "Calculate PhyML Fscore Precision & Recall..."
proformencemetrics = fscore_precision_recall.(origtrees, predtrees)
proformencemetrics = hcat(vcat.(proformencemetrics...)...)
tmpmtx = Array{Union{Missing, Float64}}(missing, length(predtreesfiles), 3)
tmpmtx[predictionexists, :] .= proformencemetrics
metricdf[!,:PhyML_fscore] .= tmpmtx[:, 1]
metricdf[!,:PhyML_precision] .= tmpmtx[:, 2]
metricdf[!,:PhyML_recall] .= tmpmtx[:, 3]

@info "Calculating PhyML Branch Depth"
mmmdepths = map(predtrees) do tr
    dists = mapinternalnodes(NewickTreeTools.network_distance, tr, tr)
    return [mean(dists), median(dists), maximum(dists)]
end
mmmdepths = hcat(mmmdepths...)'
tmpmtx = Array{Union{Missing, Float64}}(missing, length(predtreesfiles), 3)
tmpmtx[predictionexists, :] .= mmmdepths
metricdf[!,:PhyML_meandepth] .= tmpmtx[:, 1]
metricdf[!,:PhyML_mediandepth] .= tmpmtx[:, 2]
metricdf[!,:PhyML_maxdepth] .= tmpmtx[:, 3]

@info "Collecting PhyML elapsed times..."
timefiles = joinpath.("_research", "runPhyML", predtreesmsaname, "runPhyML.out");
elapsedsecs =  
    parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>"", r"\[.*"=>"$MAXSECS")) 
metricdf[!,:PhyML_elapsedsecs] .= elapsedsecs

CSV.write(joinpath(datadir(), "exp_pro", "MSAs", "MSAs-metrics-support=$(trunc(Int, SUPPORT_THRESHOLD * 100))-withphyml.csv"), metricdf)

cd(rundir)
