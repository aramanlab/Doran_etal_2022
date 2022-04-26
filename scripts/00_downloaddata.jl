using DrWatson
@quickactivate :Doran_etal_2022


# GTDB data
outputdir = joinpath(datadir(), "exp_raw", "GTDB")
mkpath(outputdir)
download("https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_r202.tree",
    joinpath(outputdir, "bac120_r202.tree"))
download("https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_taxonomy_r202.tsv",
    joinpath(outputdir, "bac120_taxonomy_r202.tsv"))