install.packages("VennDiagram") # Install VennDiagram package
library("VennDiagram") 

getwd()

# all CSB
grid.newpage()  # Move to new plotting page
vennplot = draw.triple.venn(area1 = 100,
                 area2 = 100,
                 area3 = 100,
                 n12 = 30,
                 n23 = 40,
                 n13 = 65,
                 n123 = 21,
                 c("early", "late", "middle"),
                 fill = c("yellow",  "purple", "green"),
                 lty = "blank",
                 overrideTriple=TRUE,
                 euler.d = TRUE,
                 scaled = FALSE)
pdf('./plots/donor_orthogonality/venn_csb.pdf');
grid.draw(vennplot);
dev.off()


# just gnavus & uniformis
grid.newpage()  # Move to new plotting page
vennplot = draw.triple.venn(area1 = 100,
                 area2 = 100,
                 area3 = 100,
                 n12 = 26,
                 n23 = 26,
                 n13 = 54,
                 n123 = 16,
                 c("early", "late", "middle"),
                 fill = c("yellow",  "purple", "green"),
                 lty = "blank",
                 overrideTriple=TRUE,
                 euler.d = TRUE,
                 scaled = FALSE)
pdf('./plots/donor_orthogonality/venn_br.pdf');
grid.draw(vennplot);
dev.off()
