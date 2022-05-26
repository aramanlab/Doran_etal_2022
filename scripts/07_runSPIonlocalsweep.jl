using DrWatson
@quickactivate "Doran_etal_2022"

using Glob
using NewickTreeTools
using Statistics

cd(projectdir())
inputfiles = glob(joinpath("data", "sims", "localsweep", "MSAs", "balancedtree*"))
simnames = first.(split.(basename.(inputfiles), "."))
outputdirs = joinpath.("_research", "localsweep", "runSPI", simnames)

for (inputfile, simname, outputdir) in zip(inputfiles, simnames, outputdirs)
    mkpath(outputdir)
    run(pipeline(`julia -t 4 scripts/slurm/runners/runSPI.jl \
        -i $inputfile \
        -o $outputdir \
        --nboot 100`, 
        stdout=joinpath(outputdir, "runSPI.out")))
end

sourcetrees = joinpath.(datadir(), "sims", "localsweep", "trees", replace.(basename.(inputfiles), r"-l[0-9].*"=>".nw"))

metricdf = DataFrame(
    :treefile=>sourcetrees,
    :msafile=>inputfiles,
    :msaname=>simnames,
    :outputdir=>outputdirs,
    :ntaxa=>parse.(Int, replace.(simnames, r"(.*)((?<=-t)[0-9]*)(.*)"=>s"\2")),
    :nfeatures=>parse.(Int, replace.(simnames, r"(.*)((?<=-l)[0-9]*)(.*)"=>s"\2")),
    :nbits=>parse.(Int, replace.(simnames, r"(.*)((?<=-b)[0-9]*)(.*)"=>s"\2")),
)
metricdf = transform(metricdf, :nbits => (x-> x==20 ? "AA" : "DNA") => :chartype)

### SPI METRICS

@info "Collecting SPI elapsed times..."
timefiles = joinpath.("_research", "localsweep", "runSPI", metricdf.msaname, "runSPI.out");

read.(timefiles, String) |> 
    vs->replace(vs, r"\\n"=>s"") |>
    vs->replace(vs, r"(.*)((?<=total                      1)[\s0-9\.]*)(.*)"=>s"\2")[1]


elapsedsecs = parse.(Int64, replace.(open.(lastline, timefiles, "r"), r"time_elapsed:\s*"=>""))
metricdf[!,:SPI_elapsedsecs] .= elapsedsecs

@info "Reading original source trees for SPI..."
origtrees = readnw.(read.(metricdf.treefile, String));
as_polytomy!.(origtrees, fun=n->distance(n)<1e-8);

predtreesfiles = joinpath.("_research", "localsweep", "runSPI", metricdf.msaname, metricdf.msaname .* "-supporttree.nw")
@info "Reading in SPI predicted trees..."
predtrees = readnw.(read.(predtreesfiles, String));
as_polytomy!.(predtrees, fun=n->NewickTree.support(n)<.5);
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

mkpath(joinpath(datadir(), "exp_pro", "localsweep"))
CSV.write(joinpath(datadir(), "exp_pro", "localsweep", "metricsSPI.csv"), metricdf)