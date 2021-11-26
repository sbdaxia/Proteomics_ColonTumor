library(tidyverse)
library(ggplot2)
library(ggrepel)

g12d_up <- read_csv("GSEA/Kinase enrichment analysis/Enriched in G12D.csv")
g12d_down <- read_csv("GSEA/Kinase enrichment analysis/Enriched in WT.csv")

kinase_act <- rbind(g12d_up, g12d_down)
kinase_act$ID <- ""
kinases <- c("EGFR", "MEK1", "CK2A1", "CK1D", "CK1A", "CK2B")
kinase_act$ID[kinase_act$NAME %in% kinases] <- kinase_act$NAME[kinase_act$NAME %in% kinases]
colnames(kinase_act) <- c("NAME", "NES", "p-value", "FDR", "ID")

pdf('PDF_figure/Kinase enrichment analysis.pdf',
    width = 3,
    height = 4)
ggplot(data = kinase_act, aes(x = NES, y = -log10(FDR), label = ID)) +
  geom_point(alpha = 0.6, colour = "grey", shape = 1) + 
  geom_text_repel(box.padding = 0.5, max.overlaps = Inf) +
  geom_point(data = subset(kinase_act, kinase_act$NAME %in% kinases), alpha = 2, colour = "#d32f76") + 
  xlab("NES (KRas-G12D / WT)") +
  ylab("-log10(FDR)") +
  geom_hline(yintercept = -log10(0.25)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=12, margin = margin(t = 10)),
        axis.title.y = element_text(size=12, margin = margin(r = 10)),
        axis.text = element_text(size=10),
        axis.line.y = element_line(size = 0.5),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0),
        axis.ticks.y = element_line(size = 0.5),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
  ylim(0,5)
dev.off()
