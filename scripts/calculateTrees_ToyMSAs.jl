using DrWatson
@quickactivate "Doran_etal_2022"

using Gotree_jll
using SPI
using Distances
using NewickTreeTools
using StatsPlots

cd(projectdir())

projdir = projectdir()
inputfile = datadir("exp_raw", "toyMSAs", "8x14MSA.phy")
sourcetreefile = datadir("exp_raw", "toyMSAs", "8x14sourcetree.nw")
outputdir = projectdir("_research", "toyMSA_8x14predictions")
mkpath(outputdir)
pdir = plotsdir("toyMSA_8x14_predictions")
mkpath(pdir)


run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
stdin = sourcetreefile,
stdout = joinpath(pdir, "8x14sourcetree.svg")
))


## SPI
outdir = joinpath(outputdir, "SPI")
mkpath(outdir)
run(pipeline(`julia scripts/slurm/runners/runSPI.jl \
    -i $inputfile \
    -o $outdir \
    --nboot 100`, stdout=joinpath(outdir, "SPIresults.out")))
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "8x14MSA-supporttree.nw"),
    stdout = joinpath(pdir, "8x14MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "8x14MSA-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "8x14MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "8x14MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "8x14MSA-supporttree_50pct.svg")
))


## Hamming Distances
outdir = joinpath(outputdir, "Hamming")
mkpath(outdir)

inputdf = SPI.readphylip(inputfile)
toymtx = SPI._stringcolumntocharmtx(inputdf.seqs)
dij = pairwise(Hamming(), permutedims(toymtx))
hc = hclust(dij, linkage=:average, branchorder=:optimal)
nws = SPI.nwstr(hc, inputdf.ids; labelinternalnodes=false)
open(joinpath(outdir, "hamming-tree.nw"), "w") do f
        println(f, nws)
end

boottrees = String[]
for N in 1:100
        smpcols = rand(1:14, 14)
        dij = pairwise(Hamming(), permutedims(toymtx[:, smpcols]))
        hc = hclust(dij, linkage=:average, branchorder=:optimal)
        push!(boottrees, SPI.nwstr(hc, inputdf.ids; labelinternalnodes=false))
end
open(joinpath(outdir, "hamming-boottrees.nw"), "w") do io
        for bt in boottrees
                println(io, bt)
        end
end
run(pipeline(`$(gotree()) compute support tbe --silent \
        -i $(joinpath(outdir, "hamming-tree.nw")) \
        -b $(joinpath(outdir, "hamming-boottrees.nw")) \
        -o $(joinpath(outdir, "hamming-supporttree.nw"))`,
        stderr=joinpath(outdir, "booster.log")))

run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "hamming-tree.nw"),
    stdout = joinpath(pdir, "hamming_8x14MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "hamming-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "8x14MSA-hamming-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "8x14MSA-hamming-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "8x14MSA-hamming-supporttree_50pct.svg")
))

## FastME
outdir = joinpath(outputdir, "FastME") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runFastME.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runFastME.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "8x14MSA-supporttree.nw"),
        stdout = joinpath(pdir, "FastME_8x14MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "8x14MSA-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "FastME_8x14MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "FastME_8x14MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "FastME_8x14MSA-supporttree_50pct.svg")
))

## PhyML
outdir = joinpath(outputdir, "PhyML") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runPhyML.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runPhyML.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "8x14MSA.phy-supporttree.txt"),
        stdout = joinpath(pdir, "PhyML_8x14MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "8x14MSA.phy-supporttree.txt"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "PhyML_8x14MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "PhyML_8x14MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "PhyML_8x14MSA-supporttree_50pct.svg")
))

## RAxML
outdir = joinpath(outputdir, "RAxML") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runRAxML.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runRAxML.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "RAxML-supporttree.nw"),
        stdout = joinpath(pdir, "RAxML_8x14MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "RAxML-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "RAxML_8x14MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "RAxML_8x14MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "RAxML_8x14MSA-supporttree_50pct.svg")
))

## MrBayes
outdir = joinpath(outputdir, "MrBayes") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runMrBayes.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runMrBayes.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "8x14MSA-supporttree.nw"),
        stdout = joinpath(pdir, "MrBayes_8x14MSA-supporttree.svg")
))


#######################
## 18x9 data
#######################
projdir = projectdir()
inputfile = datadir("exp_raw", "toyMSAs", "18x9MSA.phy")
sourcetreefile = datadir("exp_raw", "toyMSAs", "18x9sourcetree.nw")
outputdir = projectdir("_research", "toyMSA_18x9predictions")
mkpath(outputdir)
pdir = plotsdir("toyMSA_18x9_predictions")
mkpath(pdir)


run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
stdin = sourcetreefile,
stdout = joinpath(pdir, "18x9sourcetree.svg")
))


## SPI
outdir = joinpath(outputdir, "SPI")
mkpath(outdir)
run(pipeline(`julia scripts/slurm/runners/runSPI.jl \
    -i $inputfile \
    -o $(joinpath(outputdir, "SPI")) \
    --nboot 100`, stdout=joinpath(outdir, "SPIresults.out")))
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "18x9MSA-supporttree.nw"),
    stdout = joinpath(pdir, "18x9MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "18x9MSA-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "18x9MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "18x9MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "18x9MSA-supporttree_50pct.svg")
))

## Hamming Distances
outdir = joinpath(outputdir, "Hamming")
mkpath(outdir)
inputdf = SPI.readphylip(inputfile)
toymtx = SPI.onehotencode(inputdf.seqs)[:, 2:2:18]
dij = pairwise(Hamming(), permutedims(toymtx))
hc = hclust(dij, linkage=:average, branchorder=:optimal)
nws = SPI.nwstr(hc, inputdf.ids, labelinternalnodes=false)
open(joinpath(outdir, "hamming-tree.nw"), "w") do f
        println(f, nws)
end
heatmap(dij, ratio=1, c=:davos)
savefig(joinpath(pdir, "hamming_distance_heatmap.pdf"))

boottrees = String[]
for N in 1:100
        smpcols = rand(1:9, 9)
        dij = pairwise(Hamming(), permutedims(toymtx[:, smpcols]))
        hc = hclust(dij, linkage=:average, branchorder=:optimal)
        push!(boottrees, SPI.nwstr(hc, inputdf.ids; labelinternalnodes=false))
end
open(joinpath(outdir, "hamming-boottrees.nw"), "w") do io
        for bt in boottrees
                println(io, bt)
        end
end
run(pipeline(`$(gotree()) compute support tbe --silent \
        -i $(joinpath(outdir, "hamming-tree.nw")) \
        -b $(joinpath(outdir, "hamming-boottrees.nw")) \
        -o $(joinpath(outdir, "hamming-supporttree.nw"))`,
        stderr=joinpath(outdir, "booster.log")))

run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "hamming-tree.nw"),
    stdout = joinpath(pdir, "hamming_18x9-MSA-supporttree.svg")
))

nwtree = readnw(read(joinpath(outdir, "hamming-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "18x9MSA-hamming-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "18x9MSA-hamming-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "18x9MSA-hamming-supporttree_50pct.svg")
))


## FastME
outdir = joinpath(outputdir, "FastME") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runFastME.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runFastME.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "18x9MSA-supporttree.nw"),
        stdout = joinpath(pdir, "FastME_18x9MSA-supporttree.svg")
))
nwtree = readnw(read(joinpath(outdir, "18x9MSA-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "FastME_18x9MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "FastME_18x9MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "FastME_18x9MSA-supporttree_50pct.svg")
))

## PhyML
outdir = joinpath(outputdir, "PhyML") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runPhyML.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runPhyML.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "18x9MSA.phy-supporttree.txt"),
        stdout = joinpath(pdir, "PhyML_18x9MSA-supporttree.svg")
))
nwtree = readnw(read(joinpath(outdir, "18x9MSA.phy-supporttree.txt"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "PhyML_18x9MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "PhyML_18x9MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "PhyML_18x9MSA-supporttree_50pct.svg")
))

## RAxML
outdir = joinpath(outputdir, "RAxML") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runRAxML.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runRAxML.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "RAxML-supporttree.nw"),
        stdout = joinpath(pdir, "RAxML_18x9MSA-supporttree.svg")
))
nwtree = readnw(read(joinpath(outdir, "RAxML-supporttree.nw"), String))
as_polytomy!(nwtree, fun=n->NewickTree.support(n)<0.5)
open(joinpath(outdir, "RAxML_18x9MSA-supporttree_50pct.nw"), "w") do io
        println(io, nwstr(nwtree))
end
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
    stdin = joinpath(outdir, "RAxML_18x9MSA-supporttree_50pct.nw"),
    stdout = joinpath(pdir, "RAxML_18x9MSA-supporttree_50pct.svg")
))

## MrBayes
outdir = joinpath(outputdir, "MrBayes") 
mkpath(outdir)
run(pipeline(`julia $projdir/scripts/slurm/runners/runMrBayes.jl \
        --inputfile $inputfile \
        --outputdir $outdir \
        --model JC69 \
        --nboot 100`, 
        stdout=joinpath(outdir, "runMrBayes.out")))
run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400`,
        stdin = joinpath(outdir, "18x9MSA-supporttree.nw"),
        stdout = joinpath(pdir, "MrBayes_18x9MSA-supporttree.svg")
))
