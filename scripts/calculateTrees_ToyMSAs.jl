using DrWatson
@quickactivate "Doran_etal_2022"

using Gotree_jll

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


## Hamming Distances
outdir = joinpath(outputdir, "Hamming")
mkpath(outdir)
using SPI
using Distances
inputdf = SPI.readphylip(inputfile)
toymtx = SPI._stringcolumntocharmtx(inputdf.seqs)
dij = pairwise(Hamming(), permutedims(toymtx))
hc = hclust(dij, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, inputdf.ids, labelinternalnodes=false)
open(joinpath(outdir, "hamming-supporttree.nw"), "w") do f
        println(f, nws)
end
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "hamming-supporttree.nw"),
    stdout = joinpath(pdir, "hamming_8x14MSA-supporttree.svg")
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

## Hamming Distances
outdir = joinpath(outputdir, "Hamming")
mkpath(outdir)
using SPI
using Distances
using StatsPlots
inputdf = SPI.readphylip(inputfile)
toymtx = SPI.onehotencode(inputdf.seqs)[:, 2:2:18]
dij = pairwise(Hamming(), permutedims(toymtx))
hc = hclust(dij, linkage=:average, branchorder=:optimal)
nws = nwstr(hc, inputdf.ids, labelinternalnodes=false)
open(joinpath(outdir, "hamming-supporttree.nw"), "w") do f
        println(f, nws)
end
heatmap(dij, ratio=1, c=:davos)
savefig(joinpath(pdir, "hamming_distance_heatmap.pdf"))
plot(hc)
run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(outdir, "hamming-supporttree.nw"),
    stdout = joinpath(pdir, "hamming_8x14MSA-supporttree.svg")
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
