using DrWatson
@quickactivate "Doran_etal_2022"

using Glob
using CSV, DataFrames
using Muon

@info "Searching for emapper.annotation files..."
cd(projectdir())
files = Glob.glob(joinpath("_research", "emapperUP8881", "*emapper.annotations"))
headernames = [
    "query", "seed_ortholog", "evalue", "score",
     "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
      "Preferred_name", "GOs", "EC", "KEGG_ko",
       "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
        "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
]

@info "Starting to read annotation files..."
alldfs = Vector()
for (i, fl) in enumerate(files)
    tsdf = CSV.read(fl, DataFrame; comment="#", delim="\t", header=headernames);
    tsdf[!, "UPID"] .= first(split.(basename(fl), "."))
    push!(alldfs, tsdf)
    print("\rRead: $(i)/$(length(files))")
end
println()

fulldf = vcat(alldfs...)

alldfs = nothing

rgs = Regex.("(.*)((?<=,)[^,@]*(?=@" .* tsdf.max_annot_lvl .* "))(.*)")
# rgs = Regex.("(.*)((?<=,)[^,@]*(?=@" .* "2|Bacteria" .* "))(.*)")
prs = rgs.=>s"\2"
tsdf[!, "OGG"] .= replace.(tsdf.eggNOG_OGs, prs)


@info "Starting transform to count matrix"
OGGmatrix = fulldf |> 
    df->groupby(df, ["UPID", "OGG"]) |>
    df->combine(df, nrow => "value") |>
    df->unstack(df, :OGG, :UPID, :value; fill=0.0)
    ;
oggm = Matrix(OGGmatrix)[:, 2:end];

nstrains = mapslices(c-> sum(c.>0.0), oggm, dims=1);