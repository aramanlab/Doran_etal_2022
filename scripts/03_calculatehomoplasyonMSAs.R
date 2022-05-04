## using R v4.1.1
## in case any of these packages are not present
# install.packages("tidyverse")
# install.packages("phangorn")
# install.packages("ape")
# install.packages("treeio")
# install.packages("ggplot2")
# install.packages("data.table")
projdir <- file.path("~", "projects", "Doran_etal_2022")
setwd(projdir)
suppressMessages({
    library(tidyverse)
    library(phangorn)
    library(ape)
    library(treeio)
    library(ggplot2)
})
print("Setting up filelists...")
msafiles <- Sys.glob("data/sims/MSAs/*")
treefiles <- msafiles %>%
    basename %>%
    gsub("-l[0-9]*-b[0-9]*.phy", "", .) %>%
    paste0("data/sims/trees/", ., ".nw")

df <- data.frame(
    treefile = treefiles,
    msafile = msafiles)

nbitstotype <- c("2" = "DNA", "4" = "DNA", "20" = "AA")
df <- df %>%
    transform(treetype = basename(df$treefile) %>% strsplit("-") %>% lapply(first) %>% unlist) %>%
    transform(nbits = gsub("(.*)((?<=-b)[0-9]*)(.*)", "\\2", df$msafile, perl=TRUE)) %>%
    transform(nfeatures = gsub("(.*)((?<=-l)[0-9]*)(.*)", "\\2", df$msafile, perl=TRUE)) %>%
    transform(ntaxa = gsub("(.*)((?<=-t)[0-9]*)(.*)", "\\2", df$msafile, perl=TRUE)) %>%
    transform(chartype = nbitstotype[nbits])
print("Calculating homoplasy...")
ridf <- apply(df, 1, function(x) {
    t <- read.tree(x[1])
    m <- read.phyDat(x[2], format = "phylip", type = x[7])
    ri <- RI(t, m)
    ci <- CI(t, m)
    hm <- 1 - ri
    al <- length(unique(m))
    return(data.frame(RI = ri, CI = ci, homoplasy = hm, ntaxa_unique = al))
})
ridf <- data.table::rbindlist(ridf)
df <- cbind(df, ridf)

if (!(sum(is.na(df$RI)) == 0)) {
    stop("some of RI = NA")
}

## Save df
print("Writing data...")
datadir <- file.path(projdir, "data", "exp_pro", "MSAs")
dir.create(datadir, recursive = TRUE)
write.csv(df, 
    file.path(datadir, "MSAs-metrics.csv"),
    row.names = FALSE)

## Plot homoplasy
print("Making plots...")
plotdir <- file.path("plots", "homoplasyplots")
dir.create(plotdir, recursive = TRUE)

# homoplasy x ntaxa
ggplot(df, aes(y = homoplasy, x = ntaxa, color = nfeatures)) +
    geom_point()
ggsave(file.path(plotdir, "homoplasy_x_ntaxa_allMSAs.pdf"), useDingbats = FALSE)

# homoplasy x nbits 
ggplot(df, aes(y = homoplasy, x = as.numeric(nbits), color = nfeatures)) +
    geom_point()
ggsave(file.path(plotdir, "homoplasy_x_nbits_allMSAs.pdf"), useDingbats = FALSE)

# homoplasy x nfeatures 
ggplot(df, aes(y = homoplasy, x = nfeatures, color = treetype)) +
    geom_point()
ggsave(file.path(plotdir, "homoplasy_x_nfeatures_allMSAs.pdf"), useDingbats = FALSE)

# homoplasy x treefile
ggplot(df, aes(y = homoplasy, x = basename(treefile), color = nfeatures)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(plotdir, "homoplasy_x_treefile_allMSAs.pdf"), useDingbats = FALSE)

# nunique x nfeatures
ggplot(df, aes(x = ntaxa_unique, fill = nfeatures)) +
    geom_histogram(binwidth = 1)
ggsave(file.path(plotdir, "nunique_x_nfeatures_allMSAs.pdf"), useDingbats = FALSE)

# nfeatures x is ntaxa_unique
ggplot(df, aes(x = nfeatures, fill = (ntaxa == ntaxa_unique))) +
    geom_bar()
ggsave(file.path(plotdir, "nfeatures_x_isntaxaunique_allMSAs.pdf"), useDingbats = FALSE)
