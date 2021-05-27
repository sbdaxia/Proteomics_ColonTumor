library(tidyverse)

crcMS_diff <- read.csv("crcMS_diff.csv")
crcMS <- read.csv("2017-03-19_Haigis_5v5_Pro.csv")

diff_protein <- crcMS_diff %>% filter(p_values < 0.05 & q_values < 0.1) %>% select(Protein.Id) %>% unlist()
diff_raw_quant <- crcMS %>% filter(Protein.Id %in% diff_protein)

diff_raw_quant_matrix <- as.matrix(diff_raw_quant[, -c(1:2)])
rownames(diff_raw_quant_matrix) <- diff_raw_quant$Protein.Id

library(gplots)
library(RColorBrewer)
library(mosaic)

for (i in 1:dim(diff_raw_quant_matrix)[1]) {
  diff_raw_quant_matrix[i,1:10] <- zscore(as.numeric(diff_raw_quant_matrix[i,1:10]))
}


# all DE proteins
my_palette <- colorRampPalette(c("blue", "white", "red"))(256)

pdf('Diff_crc_proteomics_heatmap.pdf',
    width = 6,
    height = 12)
heatmap.2(diff_raw_quant_matrix,
          main = "Differentially expressed\nin colon tumor\np < 0.05 and q < 0.1",
          density.info = "none",
          key = TRUE,
          lwid = c(3,7),
          lhei = c(1,7),
          col=my_palette,
          margins = c(10,2),
          symbreaks = TRUE,
          trace = "none",
          dendrogram = "row",
          labRow = FALSE,
          ylab = "Proteins",
          cexCol = 2,
          Colv = "NA")
dev.off()

# Top 10 up/down proteins
top10diff <- crcMS_diff %>% filter(p_values < 0.05 & q_values < 0.1)
order_list <- top10diff[order(top10diff$LFC, decreasing = TRUE), ]
top10up <-  order_list[1:10, ]
top10down <- order_list[2462:2471, ]
top10 <- rbind(top10up, top10down)

top10_quant <- crcMS %>% filter(Protein.Id %in% top10$Protein.Id)

top10_quant_matrix <- as.matrix(top10_quant[, -c(1:2)])

# ID conversion
library(org.Mm.eg.db)
library(ensembldb)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)

annotations_edb <- AnnotationDbi::select(org.Mm.eg.db,
                                         keys = top10_quant$Protein.Id,
                                         columns = c("ENSEMBL", "SYMBOL"),
                                         keytype = "UNIPROT")

# manual fill-in
annotations_edb$SYMBOL[c(14,16,19)] <- c("Ig.heavy.chain.V region.36-60", "Ig.heavy.chain.V.region.441", "Ly6g")

top10_quant <- inner_join(top10_quant, annotations_edb[,c(1,3)], by = c("Protein.Id" = "UNIPROT"))

rownames(top10_quant_matrix) <- top10_quant$SYMBOL

for (i in 1:dim(top10_quant_matrix)[1]) {
  top10_quant_matrix[i,1:10] <- zscore(as.numeric(top10_quant_matrix[i,1:10]))
}


# Top 10 DE proteins
pdf('Top10Diff_crc_proteomics_heatmap.pdf',
    width = 10,
    height = 10)
heatmap.2(top10_quant_matrix,
          main = "Top 10 Differentially expressed\nin colon tumor\np < 0.05 and q < 0.1",
          density.info = "none",
          key = TRUE,
          lwid = c(3,7),
          lhei = c(1,7),
          col=my_palette,
          margins = c(8,14),
          symbreaks = TRUE,
          trace = "none",
          dendrogram = "row",
          labRow = rownames(top10_quant_matrix),
          cexCol = 2,
          Colv = "NA")
dev.off()