###################
### PREPARACIÓN ###
###################
# Limpiamos el entorno
rm(list = ls())
gc()

# Librerías
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(corrplot)
library(ggfortify)
library(tidyverse)
library(pheatmap)
library(DEGreport)
library(ggrepel)
library(tibble)
library(ggplot2)
library(openxlsx)



###########################
### DISEÑO EXPERIMENTAL ###
###########################
# Ficheros de counts
# Debe tener mínimo dos caracteres por condición, y el siguiente formato en la cabecera: "geneid" ab_1 ab_2 ab_3 cd_1 cd_2 cd_3
file_counts <- "../Counts_tomato_ribosome.tsv" 

# Fichero de anotaciones
# Dejar "" si no se se utilizan anotaciones
file_annot <- "../tomato_annot.tsv" 

# Variables
comparison <- "test"
condition1 <- "SC_Ribo_Mock"
condition2 <- "SC_Ribo_TYLCV"
rep1 <- 3
rep2 <- 3
pvt <- 0.05

# Añadir etiquetas al volcano_plot y al scatter_plot
LABELS_VOLCANO <- TRUE # TRUE o FALSE
LABELS_SCATTER <- TRUE # TRUE o FALSE



###########################
####### CONDICIONES #######
###########################
group <- factor(c(rep(condition1, rep1), rep(condition2, rep2)))
df_counts <- read.table(file_counts, sep = "\t", header = TRUE, row.names = 1, fill = TRUE, quote = "")

if (file_annot == "") {
  message("No se han encontrado las anotaciones, generando los archivos sin ellas:")
} else {
  message("Se han encontrado las anotaciones, generando los archivos con ellas:")
  df_annot <- read.table(file_annot, sep = "\t", header = TRUE, row.names = 1, fill = TRUE, quote = "")
}

keep <- grep(condition1,colnames(df_counts))
keep <- append(keep, grep(condition2,colnames(df_counts)))
df_counts <- df_counts[,keep]



##############
### DESEQ2 ###
##############
# Directorio para los resultados
dir.create(comparison, showWarnings = FALSE)
dir.create(paste0(comparison, "/graphs"), showWarnings = FALSE)

# Función
sampleCondition <- group
sampleTable <- data.frame(row.names = colnames(df_counts), condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromMatrix(countData = df_counts, colData = sampleTable, design = ~condition)
ddsHTSeq <- DESeq(ddsHTSeq)
  
# Resultados
results_df <- results(ddsHTSeq, contrast=c("condition", condition1, condition2))
write.table(as.matrix(results_df), file = paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE, na = "")

# Counts
normalized_counts <- counts(ddsHTSeq, normalized = TRUE)
write.table(as.matrix(normalized_counts), file = paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

# FDR
ddsHTSeq.res <- results(ddsHTSeq, contrast=c("condition", condition1, condition2))
ddsHTSeq.res.fdr <- ddsHTSeq.res[!is.na(ddsHTSeq.res$padj), ]
ddsHTSeq.res.fdr <- ddsHTSeq.res.fdr[ddsHTSeq.res.fdr$padj < pvt, ]
write.table(as.matrix(ddsHTSeq.res.fdr), file = paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)



################
### GRÁFICOS ###
################
# Transformación rlog
ddsHTSeq.rld <- rlogTransformation(ddsHTSeq, blind = TRUE)

# Heatmap
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
distRL <- dist(t(assay(ddsHTSeq.rld)))
mat <- as.matrix(distRL)
rownames(mat) <- colnames(mat) <- with(colData(ddsHTSeq), paste0(condition))
png(file = paste0(comparison, "/graphs/heatmap.png"), 1000, 1000, pointsize = 20)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(13, 13))
dev.off()

# PCA
png(file = paste0(comparison, "/graphs/PCA.png"), width = 2100, height = 1500, res = 300)
print(plotPCA(ddsHTSeq.rld, intgroup = "condition"))
dev.off()

# Correlograma
png(file = paste0(comparison, "/graphs/corr_plot.png"), 1000, 1000, pointsize = 20)
corrplot(cor(normalized_counts), method = "square", addCoef.col = "white")
dev.off()

# Dispersion plot
png(file = paste0(comparison, "/graphs/dsp_plot.png"), 1000, 1000, pointsize = 20)
plotDispEsts(ddsHTSeq, main="Dispersion plot")
dev.off()

# Volcano plot
data <- read.table(paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", header = TRUE, row.names = 1)
data <- data %>%
  mutate(
    Expression = case_when(
      log2FoldChange > 1 & padj < pvt ~ "Up",
      log2FoldChange < -1 & padj < pvt ~ "Down",
      TRUE ~ "No significance"),
    Significance = -log10(padj)
  )

umbral_y <- -log10(pvt)

top_up <- data %>% filter(Expression == "Up") %>% arrange(-log2FoldChange) %>% slice_head(n = 3)
top_down <- data %>% filter(Expression == "Down") %>% arrange(log2FoldChange) %>% slice_head(n = 3)
top_sig <- data %>% arrange(-Significance) %>% slice_head(n = 3)
top_genes <- bind_rows(top_up, top_down, top_sig) %>%
  rownames_to_column(var = "gene_id") %>%
  distinct(gene_id, .keep_all = TRUE)

p1 <- ggplot(data, aes(x = log2FoldChange, y = Significance)) +
  geom_point(aes(color = Expression), alpha = 0.4, size = 1.6) +
  scale_color_manual(values = c("red3", "gray60", "green3")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  geom_vline(xintercept = -1, linetype = "dotdash", color = "gray25") +
  geom_vline(xintercept = 1, linetype = "dotdash", color = "gray25") +
  geom_hline(yintercept = umbral_y, linetype = "dotdash", color = "gray25") +
  coord_cartesian(xlim = c(-15, 15), clip = "off") +
  theme(
    panel.background = element_rect(fill = "gray96"),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")) +
  xlab(expression(log[2] * "FC")) +
  ylab(expression("-log"[10] * "Padj"))
if (LABELS_VOLCANO) {
  p1 <- p1 +
    geom_label_repel(
      data = top_genes,
      aes(label = gene_id),
      size = 3,
      max.overlaps = Inf,
      segment.color = "black",
      color = "black",
      fill = "white",
      min.segment.length = 0,
      arrow = arrow(length = unit(0.01, "npc")),
      box.padding = 1.75
    ) 
  }

suppressWarnings(ggsave(paste0(comparison, "/graphs/volcano_plot.png"), plot = p1, width = 8, height = 5, dpi = 600))

# Scatter plot
counts_log10 <- log10(assay(ddsHTSeq, "counts") + 1)
cols1 <- grep(condition1, colnames(counts_log10), value = TRUE)
cols2 <- grep(condition2, colnames(counts_log10), value = TRUE)
df_means <- data.frame(
  gene_id = rownames(counts_log10),
  cond1 = rowMeans(counts_log10[, cols1, drop = FALSE]),
  cond2 = rowMeans(counts_log10[, cols2, drop = FALSE])
)

res_df <- as.data.frame(results(ddsHTSeq, contrast = c("condition", condition1, condition2))) %>%
  rownames_to_column("gene_id")

df2 <- df_means %>%
  left_join(res_df, by = "gene_id") %>%
  mutate(
    status = case_when(
      cond1 - cond2 > 0.5  ~ "Up",
      cond1 - cond2 < -0.5 ~ "Down",
      TRUE                 ~ "No significance"
    )
  )

top_up   <- res_df %>%
  filter(!is.na(padj), log2FoldChange > 1, padj < pvt) %>%
  arrange(-log2FoldChange) %>% slice_head(n = 3) %>% pull(gene_id)
top_down <- res_df %>%
  filter(!is.na(padj), log2FoldChange < -1, padj < pvt) %>%
  arrange( log2FoldChange) %>% slice_head(n = 3) %>% pull(gene_id)
top_genes_scatter <- unique(c(top_up, top_down))

p2 <- ggplot(df2, aes(x = cond1, y = cond2, color = status)) +
  geom_point(alpha = 0.4, size = 1.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray25") +
  scale_color_manual(values = c("red3", "gray60", "green3")) +
  theme(panel.background = element_rect(fill = "gray96"),
        axis.line = element_line(color = "black")) +
  labs(
    x = bquote(.(condition1) ~ "(log"[10]*" counts)"),
    y = bquote(.(condition2) ~ "(log"[10]*" counts)"),
    color = "Status"
  )

if (LABELS_SCATTER) {
  p2 <- p2 +
    geom_label_repel(
      data = df2 %>% filter(gene_id %in% top_genes_scatter),
      aes(label = gene_id),
      size = 3,
      max.overlaps = Inf,
      segment.color = "black",
      color = "black",
      fill = "white",
      arrow = arrow(length = unit(0.01, "npc")),
      box.padding = 1.75
    )
}

ggsave(paste0(comparison, "/graphs/scatter_plot.png"), plot = p2, width = 8, height = 5, dpi = 600)



#############################
### UNIÓN DE LOS ARCHIVOS ###
#############################

############################
##### CON ANOTACIONES ######
############################
if (file.exists(file_annot)) {
  # Carga los df
  df_fdr <- read.table(paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", header = TRUE, row.names = 1)
  df_norm_counts <- read.table(paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", header = TRUE, row.names = 1)
  df_results <- read.table(paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", header = TRUE, row.names = 1)
  
  # Añadir gene_id al rowname
  df_fdr <- df_fdr %>% rownames_to_column(var = "gene_id")
  df_norm_counts <- df_norm_counts %>% rownames_to_column(var = "gene_id")
  df_annot <- df_annot %>% rownames_to_column(var = "gene_id")
  df_results <- df_results %>% rownames_to_column(var = "gene_id")
  
  # Combinar los df
  df_merged_temp <- df_fdr %>%
    left_join(df_norm_counts, by = "gene_id")
  df_merged <- df_merged_temp %>%
    left_join(df_annot, by = "gene_id")
  rownames(df_merged) <- df_merged$gene_id
  df_merged$gene_id <- NULL

  # Combinar all_counts_annot
  df_all_merged_temp <- df_results %>%
    left_join(df_norm_counts, by = "gene_id")
  df_all_merged <- df_all_merged_temp %>%
    left_join(df_annot, by = "gene_id")
  rownames(df_all_merged) <- df_all_merged$gene_id
  df_all_merged$gene_id <- NULL
  
  # FDR_counts_annot
  df_sorted_merged <- df_merged[order(-df_merged$log2FoldChange), ]
  write.table(df_sorted_merged, file = paste0(comparison, "/FDR_counts_annot.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

  # FDR_counts_annot log2FC>0
  df_positive_log2fc <- df_sorted_merged[df_sorted_merged$log2FoldChange > 0, ]
  write.table(df_positive_log2fc, file = paste0(comparison, "/FDR_counts_annot_positive.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

  # FDR_counts_annot log2FC<0
  df_sorted_merged_asc <- df_merged[order(df_merged$log2FoldChange), ]
  df_negative_log2fc <- df_sorted_merged_asc[df_sorted_merged_asc$log2FoldChange < 0, ]
  write.table(df_negative_log2fc, file = paste0(comparison, "/FDR_counts_annot_negative.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
  
  
  
  ########################
  ### GENERAR UN EXCEL ###
  ########################
  df_positive_excel <- rownames_to_column(df_positive_log2fc, "Geneid")
  df_negative_excel <- rownames_to_column(df_negative_log2fc, "Geneid")
  
  wb <- createWorkbook()
  addWorksheet(wb, paste0(condition1, ">", condition2))
  writeData(wb, paste0(condition1, ">", condition2), df_positive_excel)
  addWorksheet(wb, paste0(condition2, ">", condition1))
  writeData(wb, paste0(condition2, ">", condition1), df_negative_excel)
  
  # Estilos de la cabecera y las columnas
  numCols1 <- ncol(df_positive_excel)
  numCols2 <- ncol(df_negative_excel)
  headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "left", fontSize = 12, textDecoration = "bold")
  addStyle(wb, paste0(condition1, ">", condition2), style = headerStyle, rows = 1, cols = 1:numCols1, gridExpand = TRUE)
  addStyle(wb, paste0(condition2, ">", condition1), style = headerStyle, rows = 1, cols = 1:numCols2, gridExpand = TRUE)
  setColWidths(wb, paste0(condition1, ">", condition2), cols = 1:numCols1, widths = "auto")
  setColWidths(wb, paste0(condition2, ">", condition1), cols = 1:numCols2, widths = "auto")
  
  # Congelar las filas y las columnas
  freezePane(wb, paste0(condition1, ">", condition2), firstRow = TRUE, firstCol = TRUE)
  freezePane(wb, paste0(condition2, ">", condition1), firstRow = TRUE, firstCol = TRUE)
  
  # Guardar el excel
  saveWorkbook(wb, file = paste0(comparison, "/", comparison, ".xlsx"), overwrite = TRUE)
  
  
  
  ############################
  ### GENERAR UN EXCEL ALL ###
  ############################
  df_norm_counts_excel <- rownames_to_column(df_all_merged, "Geneid")

  wb <- createWorkbook()
  addWorksheet(wb, comparison)
  writeData(wb, comparison, df_norm_counts_excel)

  # Estilos de la cabecera y las columnas
  numCols1 <- ncol(df_norm_counts_excel)
  headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "left", fontSize = 12, textDecoration = "bold")
  addStyle(wb, comparison, style = headerStyle, rows = 1, cols = 1:numCols1, gridExpand = TRUE)
  setColWidths(wb, comparison, cols = 1:numCols1, widths = "auto")

  # Congelar las filas y las columnas
  freezePane(wb, comparison, firstRow = TRUE, firstCol = TRUE)

  # Guardar el excel
  saveWorkbook(wb, file = paste0(comparison, "/", comparison, "_all.xlsx"), overwrite = TRUE)
  
  
  
  ############################
  ##### SIN ANOTACIONES ######
  ############################
} else {
  message("No se han encontrado las anotaciones, generando los archivos sin ellas")
  
  # Carga los df
  df_fdr <- read.table(paste0(comparison, "/DESeq2_FDR.tsv"), sep = "\t", header = TRUE, row.names = 1)
  df_norm_counts <- read.table(paste0(comparison, "/DESeq2_counts.tsv"), sep = "\t", header = TRUE, row.names = 1)
  df_results <- read.table(paste0(comparison, "/DESeq2_results.tsv"), sep = "\t", header = TRUE, row.names = 1)
  
  # Añadir gene_id al rowname
  df_fdr <- df_fdr %>% rownames_to_column(var = "gene_id")
  df_norm_counts <- df_norm_counts %>% rownames_to_column(var = "gene_id")
  df_results <- df_results %>% rownames_to_column(var = "gene_id")
  
  # Combinar los df
  df_merged <- df_fdr %>%
    left_join(df_norm_counts, by = "gene_id")
  rownames(df_merged) <- df_merged$gene_id
  df_merged$gene_id <- NULL
  
  # Combinar all_counts_annot
  df_all_merged <- df_results %>%
    left_join(df_norm_counts, by = "gene_id")
  rownames(df_all_merged) <- df_all_merged$gene_id
  df_all_merged$gene_id <- NULL
  
  # FDR_counts
  df_sorted_merged <- df_merged[order(-df_merged$log2FoldChange), ]
  write.table(df_sorted_merged, file = paste0(comparison, "/FDR_counts.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
  
  # FDR_counts log2FC>0
  df_positive_log2fc <- df_sorted_merged[df_sorted_merged$log2FoldChange > 0, ]
  write.table(df_positive_log2fc, file = paste0(comparison, "/FDR_counts_positive.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
  
  # FDR_counts log2FC<0
  df_sorted_merged_asc <- df_merged[order(df_merged$log2FoldChange), ]
  df_negative_log2fc <- df_sorted_merged_asc[df_sorted_merged_asc$log2FoldChange < 0, ]
  write.table(df_negative_log2fc, file = paste0(comparison, "/FDR_counts_negative.tsv"), sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

  
  
  ########################
  ### GENERAR UN EXCEL ###
  ########################
  df_positive_excel <- rownames_to_column(df_positive_log2fc, "Geneid")
  df_negative_excel <- rownames_to_column(df_negative_log2fc, "Geneid")
  
  wb <- createWorkbook()
  addWorksheet(wb, paste0(condition1, ">", condition2))
  writeData(wb, paste0(condition1, ">", condition2), df_positive_excel)
  addWorksheet(wb, paste0(condition2, ">", condition1))
  writeData(wb, paste0(condition2, ">", condition1), df_negative_excel)
  
  # Estilos de la cabecera y las columnas
  numCols1 <- ncol(df_positive_excel)
  numCols2 <- ncol(df_negative_excel)
  headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "left", fontSize = 12, textDecoration = "bold")
  addStyle(wb, paste0(condition1, ">", condition2), style = headerStyle, rows = 1, cols = 1:numCols1, gridExpand = TRUE)
  addStyle(wb, paste0(condition2, ">", condition1), style = headerStyle, rows = 1, cols = 1:numCols2, gridExpand = TRUE)
  setColWidths(wb, paste0(condition1, ">", condition2), cols = 1:numCols1, widths = "auto")
  setColWidths(wb, paste0(condition2, ">", condition1), cols = 1:numCols2, widths = "auto")
  
  # Congelar las filas y las columnas
  freezePane(wb, paste0(condition1, ">", condition2), firstRow = TRUE, firstCol = TRUE)
  freezePane(wb, paste0(condition2, ">", condition1), firstRow = TRUE, firstCol = TRUE)
  
  # Guardar el excel
  saveWorkbook(wb, file = paste0(comparison, "/", comparison, ".xlsx"), overwrite = TRUE)


  
  ############################
  ### GENERAR UN EXCEL ALL ###
  ############################
  df_norm_counts_excel <- rownames_to_column(df_all_merged, "Geneid")

  wb <- createWorkbook()
  addWorksheet(wb, comparison)
  writeData(wb, comparison, df_norm_counts_excel)

  # Estilos de la cabecera y las columnas
  numCols1 <- ncol(df_norm_counts_excel)
  headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD", halign = "left", fontSize = 12, textDecoration = "bold")
  addStyle(wb, comparison, style = headerStyle, rows = 1, cols = 1:numCols1, gridExpand = TRUE)
  setColWidths(wb, comparison, cols = 1:numCols1, widths = "auto")

  # Congelar las filas y las columnas
  freezePane(wb, comparison, firstRow = TRUE, firstCol = TRUE)

  # Guardar el excel
  saveWorkbook(wb, file = paste0(comparison, "/", comparison, "_all.xlsx"), overwrite = TRUE)
}
