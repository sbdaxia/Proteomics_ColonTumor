diff_id <- c()

for (i in 1:dim(crcMS_diff)[1]) {
  if (crcMS_diff$p_values[i] < 0.05 & crcMS_diff$q_values[i] < 0.1) {
    diff_id <- c(diff_id, crcMS_diff$X1[i])
  }
}

diff_index <- c()

for (i in 1: length(diff_id)) {
  exactname <- paste('^',diff_id[i], '$', sep = '')
  diff_index <- c(diff_index, grep(exactname, protein_quant$Protein.Id))
}

diff_raw_quant <- protein_quant[diff_index,]
library(mosaic)
for (i in 1:dim(diff_raw_quant)[1]) {
  diff_raw_quant[i, 3:12]<- zscore(as.numeric(as.character(diff_raw_quant[i, 3:12])))
}


diff_raw_quant_matrix <- as.matrix(diff_raw_quant[, -c(1:2)])
rownames(diff_raw_quant_matrix) <- diff_raw_quant$Protein.Id

library(gplots)
library(RColorBrewer)

my_palette <- colorRampPalette(c("red", "white", "blue"))(1000)

png('Diff_crc_proteomics_heatmap.png',
    width = 600,
    height = 1400,
    res = 200,
    pointsize = 8)
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
