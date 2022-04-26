using DrWatson
@quickactivate :Doran_etal_2022

## GTDB Data ##

# strip internal node names from GTDB tree because none of the newicktree
# parsers I've tried can handle the quotes and semicolons inside.
# also I am only using this tree for its topology, so the names
# will be meaningless in regards to the simulated alignments.
gtdbstring = read(joinpath(datadir(), "exp_raw", "GTDB", "bac120_r202.tree"), String)
stripped_gtdbstring = replace(gtdbstring, r"\'[^\']*\'" => s"")
mkpath(joinpath(datadir(), "exp_pro", "GTDB"))
open(joinpath(datadir(), "exp_pro", "GTDB", "bac120_r202.nw"), "w") do io
    write(io, stripped_gtdbstring)
end