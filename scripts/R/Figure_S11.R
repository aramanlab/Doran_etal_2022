library(RColorBrewer)
library(tidyverse)
library(ape)
library(ggplot2)
library(ggtree)

projectdir = "/Users/ben/projects/Doran_etal_2022"
ddir = file.path(projectdir, "data", "exp_pro", "BB673")
pdir = file.path(projectdir, "plots", "fig_s11")

t = read.tree(file.path(ddir, "inferrednewicktree_UP7047.nw"))

metadata = read.table(file.path(ddir, "BB673_obsdata.tsv"), sep="\t", header=1)
metadata %>% head
metadata$Family = metadata$family
metadata$Family[grepl("vulgatus", metadata$species)] = "Bacteroidaceae"
metadata$Family["" == metadata$Family] = "<unknown family>"
mdata = metadata[match(t$tip.label, metadata$ID), ]
md = column_to_rownames(tibble(mdata), "ID")

md[, c("Family"), drop=FALSE]

palette = brewer.pal(11,'Set3');
brks = c(
  "Lachnospiraceae",
  "Bacteroidaceae",
  "Bifidobacteriaceae",
  "Tannerellaceae",
  "Prevotellaceae",
  "Rikenellaceae",
  "<unknown family>",
  "Oscillospiraceae",
  "Erysipelotrichaceae",
  "Odoribacteraceae",
  "Enterobacteriaceae"
)

p = ggtree(t, layout="fan", open.angle = 20)
gheatmap(p, md[, c("Family"), drop=FALSE], 
    offset=.8, 
    width=.2,
    colnames_angle=45, 
    colnames_offset_y=-.2,
) +
scale_fill_manual(breaks=brks, values=palette, name="Family")
ggsave(file.path(pdir, "fulltreepaintedbyfamily.pdf"))

