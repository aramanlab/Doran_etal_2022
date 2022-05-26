using DrWatson
@quickactivate "Doran_etal_2022"
using CSV, DataFrames
using Gotree_jll, Goalign_jll, SeqGen_jll

treedir = joinpath(datadir(), "sims", "localsweep", "trees")
mkpath(treedir)

numleaves = [16, 32, 64, 128, 256, 512, 1024]
nfeatures = [16, 32, 64, 128, 256, 512, 1024]
baldepth = Int.(log2.(numleaves));
const seed = 42

## Make randomly generated trees ##
# branch lengths for each random tree are selected from an exponential distribution
# trees from GTDB are selected based on number of leaves an internal node posesses
for nl in numleaves
    run(pipeline(`$(gotree()) generate yuletree --seed $(seed) -l $nl`,
        `$(gotree()) brlen clear`, 
        `$(gotree()) brlen setmin -l 1.0 -o $(treedir)/yuletree-t$(nl)-s$(seed).nw`
        ))
end
for (d, nl) in zip(baldepth, numleaves)
    run(pipeline(`$(gotree()) generate balancedtree --seed $(seed) -d $d`,
    `$(gotree()) brlen clear`, 
    `$(gotree()) brlen setmin -l 1.0 -o $(treedir)/balancedtree-t$(nl)-s$(seed).nw`
    ))
end
# for nl in numleaves
#     run(`$(gotree()) generate uniformtree --seed $(seed) -l $nl -o $(treedir)/uniformtree-t$(nl)-s$(seed).nw`)
# end

## Make plots of simulated source trees ##
srctreeplotsdir = joinpath(projectdir(), "plots", "localtrees")
rm(srctreeplotsdir, force=true, recursive=true)
mkpath(srctreeplotsdir)
for tf in readdir(treedir)
    fn = first(split(tf, "."))
    run(pipeline(`$(gotree()) draw svg -w 400 -H 400`,
    stdin = joinpath(treedir, tf),
    stdout = joinpath(srctreeplotsdir, fn * ".svg")
    ))
end

## Make MSAs from generated trees ##
msadir = joinpath(datadir(), "sims", "localsweep", "MSAs")
rm(msadir, force=true, recursive=true)
mkpath(msadir)

# Nbits = 2 (binary)
for Nf in nfeatures
    mkpath(msadir)
    for f in readdir(treedir)
        fn = first(split(f, "."))
        run(pipeline(`$(seqgen()) -q -z$(seed) -s $(2/Nf) -or -l$Nf -mHKY -f0.5,0.0,0.0,0.5`,
            stdin=joinpath(treedir, f),
            stdout=joinpath(msadir, fn * "-l$Nf-b2.phy")))
    end
end
# # Nbits = 4 (nucleotide)
# # HKY with no other parameters, equivalent to JC69 model
# for Nf in nfeatures
#     mkpath(msadir)
#     for f in readdir(treedir)
#         fn = first(split(f, "."))
#         run(pipeline(`$(seqgen()) -q -z$(seed) -s $(gtdbmsalength/Nf) -or -l$Nf -mHKY`, 
#             stdin=joinpath(treedir, f),
#             stdout=joinpath(msadir, fn * "-l$Nf-b4.phy")))
#     end
# end
# # Nbits = 20 (aminoacids)
# for Nf in nfeatures
#     mkpath(msadir)
#     for f in readdir(treedir)
#         fn = first(split(f, "."))
#         run(pipeline(`$(seqgen()) -q -z$(seed) -s $(gtdbmsalength/Nf) -or -l$Nf -mWAG`,
#             stdin=joinpath(treedir, f),
#             stdout=joinpath(msadir, fn * "-l$Nf-b20.phy")))
#     end
# end

