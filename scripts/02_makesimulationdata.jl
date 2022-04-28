using DrWatson
@quickactivate :Doran_etal_2022
using Gotree_jll, Goalign_jll, SeqGen_jll

treedir = joinpath(datadir(), "sims", "trees")
mkpath(treedir)

numleaves = [16, 32, 64, 128]
baldepth = Int.(log2.(numleaves));
seed = 123456

## Make randomly generated trees ##
# branch lengths for each random tree are selected from an exponential distribution
# trees from GTDB are selected based on number of leaves an internal node posesses
for nl in numleaves
    run(`$(gotree()) generate yuletree --seed $(seed) -l $nl -o $(treedir)/yuletree-t$(nl)-s$(seed).nw`)
end
for (d, nl) in zip(baldepth, numleaves)
    run(`$(gotree()) generate balancedtree --seed $(seed) -d $d -o $(treedir)/balancedtree-t$(nl)-s$(seed).nw`)
end
for nl in numleaves
    run(`$(gotree()) generate caterpillartree --seed $(seed) -l $nl -o $(treedir)/caterpillartree-t$(nl)-s$(seed).nw`)
end
for nl in numleaves
    run(`$(gotree()) generate uniformtree --seed $(seed) -l $nl -o $(treedir)/uniformtree-t$(nl)-s$(seed).nw`)
end
for nl in numleaves
    run(`$(gotree()) generate startree --seed $(seed) -l $nl -o $(treedir)/startree-t$(nl)-s$(seed).nw`)
    run(`$(gotree()) resolve -i $(treedir)/startree-t$(nl)-s$(seed).nw -o $(treedir)/startree_resolved-t$(nl)-s$(seed).nw`)
    mv("$(treedir)/startree_resolved-t$(nl)-s$(seed).nw", "$(treedir)/startree-t$(nl)-s$(seed).nw", force=true)
end

## Make trees from GTDB subtrees ##
using NewickTree
using NewickTree: getleaves
using DataStructures: counter, inc!

# read metadata
df = CSV.read(joinpath(datadir(), "exp_raw", "GTDB", "bac120_taxonomy_r202.tsv"), DataFrame, header=false);
bacdf = df |>
    df->transform(df, 
        :Column1 => :ID, 
        :Column2 => ByRow((c)->split(c, ";")) => [:domain, :phylum, :class, :order, :family, :genus, :species]) |>
    df->select(df,[:ID, :domain, :phylum, :class, :order, :family, :genus, :species]);
# read full tree
gtdbtree = readnw(read(joinpath(datadir(), "exp_pro", "GTDB", "bac120_r202.nw"), String));

# extract 128 taxa subtree
gtdbtree128 = filter(prewalk(gtdbtree)) do node
    length(getleaves(node)) == 128
end |> first;
rootname128 = bacdf[indexin(name.(getleaves(gtdbtree128)), bacdf.ID), :family] |> unique |> first
open(joinpath(treedir,"gtdb-$rootname128-t128.nw"), "w") do io
    # strip all node names and branch supports because seqgen doesn't like them
    treestring = replace(nwstr(gtdbtree128), r"[a-zA-Z0-9._]*(?=:)" => s"") 
    # replace with unique tip names, tree topology is maintained
    c = counter(String)
    treestring = replace(treestring, r"(?<=[\(,])(?=:)" => m->"Tip$(lpad(inc!(c,m), 3, "0"))")
    write(io, treestring * "\n")
end

# extract 64 taxa subtree
gtdbtree64 = filter(prewalk(gtdbtree)) do node
    length(getleaves(node)) == 64
end |> first;
rootname64 = bacdf[indexin(name.(getleaves(gtdbtree64)), bacdf.ID), :genus] |> unique |> first
open(joinpath(treedir,"gtdb-$rootname64-t64.nw"),"w") do io
    # strip all node names and branch supports because seqgen doesn't like them
    treestring = replace(nwstr(gtdbtree64), r"[a-zA-Z0-9._]*(?=:)" => s"") 
    # replace with unique tip names, tree topology is maintained
    c = counter(String)
    treestring = replace(treestring, r"(?<=[\(,])(?=:)" => m->"Tip$(lpad(inc!(c,m), 2, "0"))")
    write(io, treestring * "\n")
end

# extract 32 taxa subtree
gtdbtree32 = filter(prewalk(gtdbtree)) do node
    length(getleaves(node)) == 32
end |> first;
rootname32 = bacdf[indexin(name.(getleaves(gtdbtree32)), bacdf.ID), :genus] |> unique |> first
open(joinpath(treedir,"gtdb-$rootname32-t32.nw"),"w") do io
    # strip all node names and branch supports because seqgen doesn't like them
    treestring = replace(nwstr(gtdbtree32), r"[a-zA-Z0-9._]*(?=:)" => s"") 
    # replace with unique tip names, tree topology is maintained
    c = counter(String)
    treestring = replace(treestring, r"(?<=[\(,])(?=:)" => m->"Tip$(lpad(inc!(c,m), 2, "0"))")
    write(io, treestring * "\n")
end

# extract 16 taxa subtree
gtdbtree16 = filter(prewalk(gtdbtree)) do node
    length(getleaves(node)) == 16
end |> x->getindex(x, 8); # choosing 8th tree to sample something other than psudomonas
rootname16 = bacdf[indexin(name.(getleaves(gtdbtree16)), bacdf.ID), :genus] |> unique |> first
open(joinpath(treedir,"gtdb-$rootname16-t16.nw"),"w") do io
    # strip all node names and branch supports because seqgen doesn't like them
    treestring = replace(nwstr(gtdbtree16), r"[a-zA-Z0-9._]*(?=:)" => s"") 
    # replace with unique tip names, tree topology is maintained
    c = counter(String)
    treestring = replace(treestring, r"(?<=[\(,])(?=:)" => m->"Tip$(lpad(inc!(c,m), 2, "0"))")
    write(io, treestring * "\n")
end

## Make MSAs from generated trees ##
nfeatures = [100, 1000, 10000, 50000]
gtdbmsalength = 34744
seed = 123456

msadir = joinpath(datadir(), "sims", "MSAs")
rm(msadir, force=true, recursive=true)
mkpath(msadir)

# Nbits = 2 (binary)
for Nf in nfeatures
    mkpath(msadir)
    for f in readdir(treedir)
        fn = first(split(f, "."))
        run(pipeline(`$(seqgen()) -q -z$(seed) -s $(gtdbmsalength/Nf) -or  -l$Nf -mHKY -f0.5,0.0,0.0,0.5`,
            stdin=joinpath(treedir, f),
            stdout=joinpath(msadir, fn * "-l$Nf-b2.phy")))
    end
end
# Nbits = 4 (nucleotide)
# HKY with no other parameters, equivalent to JC69 model
for Nf in nfeatures
    mkpath(msadir)
    for f in readdir(treedir)
        fn = first(split(f, "."))
        run(pipeline(`$(seqgen()) -q -z$(seed) -s $(gtdbmsalength/Nf) -or -l$Nf -mHKY`, 
            stdin=joinpath(treedir, f),
            stdout=joinpath(msadir, fn * "-l$Nf-b4.phy")))
    end
end
# Nbits = 20 (aminoacids)
for Nf in nfeatures
    mkpath(msadir)
    for f in readdir(treedir)
        fn = first(split(f, "."))
        run(pipeline(`$(seqgen()) -q -z$(seed) -s $(gtdbmsalength/Nf) -or -l$Nf -mWAG`,
            stdin=joinpath(treedir, f),
            stdout=joinpath(msadir, fn * "-l$Nf-b20.phy")))
    end
end

## Make plots of simulated source trees ##
srctreeplotsdir = joinpath(plotdir(), "simtrees")
rm(srctreeplotsdir, force=true, recursive=true)
mkpath(srctreeplotsdir)
for tf in readdir(treedir)
    fn = first(split(tf, "."))
    run(pipeline(`$(gotree()) draw svg -c -w 400 -H 400 --no-tip-labels`,
    stdin = joinpath(treedir, tf),
    stdout = joinpath(srctreeplotsdir, fn * ".svg")
    ))
end