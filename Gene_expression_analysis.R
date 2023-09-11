getwd()
setwd("/home/alumno20/TFM_IAM/samples")

#librerías a utilizar
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library("VennDiagram", lib.loc = "/home/alumno20/TFM_IAM")
library("clValid", lib.loc = "/home/alumno20/TFM_IAM")


#Cargamos la matriz de conteo de datos genéticos
Count_Matrix <- `Galaxy251.[Multi.Join_on_data_248,_data_246,_and_others]`
#ajustamos la matriz para un entorno más facil de utilizar
colnames(Count_Matrix) <- Count_Matrix[1,]
colnames(Count_Matrix) <- c("Gene_ID","SRR5985090", "SRR5985103", "SRR5985112", "SRR5985116","SRR5985091","SRR5985092","SRR5985094","SRR5985095","SRR5985096","SRR5985104","SRR5985105","SRR5985106","SRR5985107","SRR5985108","SRR5985110", "SRR5985119")
Count_Matrix <- Count_Matrix[-1,]

rownames(Count_Matrix)<- Count_Matrix[,1]
Count_Matrix<- Count_Matrix[,-1]

head(Count_Matrix)
str(Count_Matrix)
#convertimos en variable snumericas
Count_Matrix <- Count_Matrix %>%
  mutate_if(is.character, as.numeric)

#cargamos y modificamos los datos que contienen información sobre las muestras
rownames(sample_info)<- sample_info[,1]
sample_info<- sample_info[,-1]
head(sample_info)
str

#Nos aseguramos que coincidan los encabezados de la matriz de conteo y el indice del guión de información de las muestras
all(colnames(Count_Matrix) %in% rownames(sample_info))

#Convertimos nuestras variables a factor
sample_info$IAM <- factor(sample_info$IAM)
sample_info$Paciente <- factor(sample_info$Paciente)

#Creamos el datatset para DESeq
dds <- DESeqDataSetFromMatrix(countData = Count_Matrix, 
                       colData = sample_info,
                       design = ~ IAM + Paciente)
dds

#ahora generamos un filtrado para eliminar aquellos conteos de genes con poca información y obtener un dataset más limpio
keep <- rowSums(counts(dds)) >=10
keep
dds <- dds[keep,]
dds

#Establecemos el nivel del factor
dds$IAM <- factor(dds$IAM, levels = c("STEMI", "NSTEMI"))
dds$IAM

#Ejecutamos el DESeq
dds2 <- dds
dds2 <- DESeq(dds2)

#Exploramos los resultados
res <- results(dds2)
res
summary(res)

resultsNames(dds2)
results(dds2, contrast = c("IAM", "STEMI", "NSTEMI"))

#Cambiamos los resultados a un dataframe
dds_res <- res
dds_res <- as.data.frame(dds_res)

#Ordenamos la tabla de resultados en orden creciente de p-valor
dds_res_sort <- dds_res[order(dds_res$pvalue),]
head(dds_res_sort)

#Extraemos los genes mayormente expresados en el dataset

##Paso 1: Filtrado en el padj
filtered <- dds_res %>% filter(dds_res$pvalue > 0.05)

##Paso 2: Filtrado basado en Foldchange, usamos límite de 1
filtered <- filtered %>% filter(abs(filtered$log2FoldChange)>1)

#Vemos que las dimensiones del dataset han variado
dim(filtered) 

#Sacamos los genes que pertenecen a cada tratamiento
Factor_IAM <- subset(filtered, padj > 0.05)

# Identify genes associated with a specific factor
factor_STEMI <- Factor_IAM[dds$IAM == "STEMI", ]
factor_NSTEMI <- Factor_IAM[dds$IAM == "NSTEMI", ]


#guardamos el resultado
write.csv(dds_res, "dds_res.all.csv")
write.csv(filtered, "filtered.dds_res.csv")
#guardamos los datos normalizados
normalized_counts <- counts(dds2,normalized=TRUE)
head(normalized_counts)
write.csv(normalized_counts, "normalized_counts.csv")

#analizamos la expresión de genes
log2dc <- results(dds2)$log2FoldChange
##treshold
ldc_threshold <- 1.0
##Identificamos sobreexpresión
overexpressed_genes <- rownames(results(dds2))[log2dc > ldc_threshold]
overexpressed_genes_table <- as.data.frame(overexpressed_genes)

overexpressed_dataset <- results(dds2)[log2dc > ldc_threshold]

balanced_genes <- rownames(results(dds2))[abs(log2dc) <= ldc_threshold]
balanced_genes_table <- as.data.frame(balanced_genes)

print(overexpressed_genes)
print(balanced_genes)

###genes overexpressed from STEMI
genes_OP_STEMI <- overexpressed_genes_table[dds2$IAM == "STEMI", ]
###genes overexpressed from NSTEMI
genes_OP_NSTEMI <- overexpressed_genes_table[dds2$IAM == "NSTEMI", ]


###genes balanceados STEMI
genes_balanced_STEMI <- balanced_genes_table[dds2$IAM == "STEMI", ]
###genes balanceados NSTEMI
genes_balanced_NSTEMI <- balanced_genes_table[dds2$IAM == "NSTEMI", ]

#clustering
k <- 3
kmeans_result <- kmeans(filtered, centers = k)


cluster_matrix <- matrix(kmeans_result$cluster, ncol = 1)

clust_valid <- clValid(cluster_matrix, 
                       measures = c("kl", "ch", "hartigan", "ccc", "scott"), 
                       validation = "internal")

#Plots

##Plot de dispersión
plotDispEsts(dds2)

##Plot MA
plotMA(res)


##PCA plot
###Utilizaremos un plot a apertir de un PCA para poder explicar la varianza en la expresión de genes 
vsd <- vst(dds2, blind = FALSE)
plotPCA(vsd, intgroup=c("IAM", "Paciente"))


#Heatmaps
library(pheatmap)


##Heatmap para la matriz de distancia (con clustering) basado en la normalización del conteo

###Generamos la matriz de distancia
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)

###Hacemos el esquema de colores
colorsDDS <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

###Construimos el mapa
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorsDDS)

##Heatmap para los datos de conteo normalizados transformados a logaritmo
log_counts <- dds_res[order(dds_res$padj),][1:10,]
log_counts <- row.names(log_counts)
log_counts

rld <- rlog(dds2,blind=FALSE)
pheatmap(assay(rld)[log_counts,])

#A este mapa podemos añadirle información sobre las muestras
annot_information <- as.data.frame(colData(dds2)[,c("IAM", "Paciente")])
pheatmap(assay(rld)[log_counts,],
         annotation_col = annot_information)

##Heatmap de los Z-scores para la desviación estándar 
Z_score <- function(x) {(x-mean(x)) / sd(x)}

zscore_all <- t(apply(normalized_counts, 1, Z_score))
zscore_subset <-  zscore_all[log_counts,]
pheatmap(zscore_subset)


####Enriquecimiento funcional#####
library(clusterProfiler)
library(org.Hs.eg.db)

gene_to_test <- rownames(filtered[filtered$log2FoldChange > 1,])

gene_list <- sub("\\..*", "", gene_to_test)


GO_result  <- enrichGO(gene = gene_list,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = FALSE,
                       pvalueCutoff = 0.05)

enrichment_results <- as.data.frame(GO_result)

dotplot(GO_result)


##Fligner DESeq2 results
tofligner <- normalized_counts[dds_res$log2FoldChange>1,]
tofligner <- as.data.frame(tofligner)
conditions <- tofligner[,dds2$IAM]
fligner.test(tofligner ~ conditions)

#Buscamos las proteinas relacionadas con estos ID de GO
GO_genes <- enrichment_results$ID
protein_ids <- mapIds(org.Hs.eg.db,
                      keys = GO_genes,
                      column = "ENSEMBLPROT",
                      keytype = "GO")
protein_ids <- as.data.frame(protein_ids)
protein_ids <- na.omit(protein_ids)
protein_ids <- protein_ids$protein_ids
protein_ids <- as.data.frame(protein_ids)


#Buscamos las proteinas relacionadas con estos ID de GO de muestras STEMI
GO_genes_STEMI <- enrichment_results_STEMI$ID
protein_ids_STEMI <- mapIds(org.Hs.eg.db,
                      keys = GO_genes_STEMI,
                      column = "ENSEMBLPROT",
                      keytype = "GO")
protein_ids_STEMI <- as.data.frame(protein_ids_STEMI)
protein_ids_STEMI <- na.omit(protein_ids_STEMI)
protein_ids_STEMI <- protein_ids_STEMI$protein_ids_STEMI
protein_ids_STEMI <- as.data.frame(protein_ids_STEMI)


#Buscamos las proteinas relacionadas con estos ID de GO de muestras STEMI
GO_genes_NSTEMI <- enrichment_results_NSTEMI$ID
protein_ids_NSTEMI <- mapIds(org.Hs.eg.db,
                            keys = GO_genes_NSTEMI,
                            column = "ENSEMBLPROT",
                            keytype = "GO")
protein_ids_NSTEMI <- as.data.frame(protein_ids_NSTEMI)
protein_ids_NSTEMI <- na.omit(protein_ids_NSTEMI)
protein_ids_NSTEMI <- protein_ids_NSTEMI$protein_ids_NSTEMI
protein_ids_NSTEMI <- as.data.frame(protein_ids_NSTEMI)


#GO genes para Cytoscape

gene_names_STEMI_CYTOS <- enrichment_results_STEMI$geneID
flat_genes_STEMI <- unlist(gene_names_STEMI_CYTOS)


#Guardamos los datos del enriquecimento
write.csv(enrichment_results, "enrichment_results.csv")
write.csv(enrichment_results_STEMI, "enrichment_results_STEMI.csv")
write.csv(enrichment_results_NSTEMI, "enrichment_results_NSTEMI.csv")
write.csv(gene_list, "gene_list.csv")
write.csv(protein_ids, "protein_ids.csv", sep = "/t")
write.csv(protein_ids_STEMI, "protein_ids_STEMI.csv", sep = "/t")
write.csv(protein_ids_NSTEMI, "protein_ids_NSTEMI.csv", sep = "/t")
 
#Enrichment analysis STEMI
gene_to_test_STEMI <- rownames(factor_STEMI[factor_STEMI$log2FoldChange > 1,])

gene_list_STEMI <- sub("\\..*", "", gene_to_test_STEMI)


GO_result_STEMI  <- enrichGO(gene = gene_list_STEMI,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP",
                       readable = FALSE,
                       pvalueCutoff = 0.2)

enrichment_results_STEMI <- as.data.frame(GO_result_STEMI)

dotplot(GO_result_STEMI)

#Enrichment analysis NSTEMI
gene_to_test_NSTEMI <- rownames(factor_NSTEMI[factor_NSTEMI$log2FoldChange > 1,])

gene_list_NSTEMI <- sub("\\..*", "", gene_to_test_NSTEMI)



GO_result_NSTEMI  <- enrichGO(gene = gene_list_NSTEMI,
                             OrgDb = "org.Hs.eg.db",
                             keyType = "ENSEMBL",
                             ont = "BP",
                             readable = FALSE,
                             pvalueCutoff = 0.2)

enrichment_results_NSTEMI <- as.data.frame(GO_result_NSTEMI)

dotplot(GO_result_NSTEMI)

#Representamos 

ggplot(enrichment_results, aes(x= ID, y = p.adjust)) +
  geom_bar(stat = "identity") +
  xlab("Gene Set") +
  ylab("Enrichment padjust")

ggplot(enrichment_results, aes(x= ID, y = GeneRatio)) +
  geom_bar(stat = "identity") +
  xlab("Gene Set") +
  ylab("GeneRatio")

##Buscamos las proteinas relacionadas con estos ID de GO para STEMI
GO_genes_STEMI <- enrichment_results_STEMI$ID
protein_ids_STEMI <- mapIds(org.Hs.eg.db,
                      keys = GO_genes_STEMI,
                      column = "ENSEMBLPROT",
                      keytype = "GO")
protein_ids_STEMI <- as.data.frame(protein_ids_STEMI)
protein_ids_STEMI <- na.omit(protein_ids_STEMI)
protein_ids_STEMI <- protein_ids_STEMI$protein_ids_STEMI
protein_ids_STEMI <- as.data.frame(protein_ids_STEMI)

##Buscamos las proteinas relacionadas con estos ID de GO para NSTEMI
GO_genes_NSTEMI <- enrichment_results_NSTEMI$ID
protein_ids_NSTEMI <- mapIds(org.Hs.eg.db,
                            keys = GO_genes_NSTEMI,
                            column = "ENSEMBLPROT",
                            keytype = "GO")
protein_ids_NSTEMI <- as.data.frame(protein_ids_NSTEMI)
protein_ids_NSTEMI <- na.omit(protein_ids_NSTEMI)
protein_ids_NSTEMI <- protein_ids_NSTEMI$protein_ids_NSTEMI
protein_ids_NSTEMI <- as.data.frame(protein_ids_NSTEMI)

###Guardamos los resultados
write.csv(protein_ids_STEMI, "protein_ids_STEMI.csv", sep = "/t")
write.csv(protein_ids_NSTEMI, "protein_ids_NSTEMI.csv", sep = "/t")

##volcano plot

### creamos el plot
volcano_plot <- ggplot(dds_res, aes(x = dds_res$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(dds_res$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()

### visualizar
print(volcano_plot)

##Hacemos un subset para analizarlos por partes

IAM <- subset(dds_res, padj > 0.05)

volcano_STEMI <- IAM[dds$IAM == "STEMI", ]
volcano_NSTEMI <- IAM[dds$IAM == "NSTEMI", ]

###volcano para STEMI
volcano_plot_STEMI <- ggplot(volcano_STEMI, aes(x = volcano_STEMI$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(volcano_STEMI$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()
print(volcano_plot_STEMI)

###volcano para NSTEMI
volcano_plot_NSTEMI <- ggplot(volcano_NSTEMI, aes(x = volcano_NSTEMI$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(volcano_NSTEMI$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()
print(volcano_plot_NSTEMI)

###volcano con datos filtrados a un lfc = 1 y pvalue > 0.05

volcano_plot_filtered <- ggplot(filtered, aes(x = filtered$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(filtered$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()
print(volcano_plot_filtered)

###plot volcano sobreexresión

#hacemos un sesgo para sacar muestras 
fold_change_threshold <- 1.0


overexpressed_genes_volcano <- dds_res[abs(dds_res$log2FoldChange) > fold_change_threshold , ]

###genes overexpressed from STEMI
genes_OP_STEMI_volcano <- overexpressed_genes_volcano[dds2$IAM == "STEMI", ]
###genes overexpressed from NSTEMI
genes_OP_NSTEMI_volcano <- overexpressed_genes_volcano[dds2$IAM == "NSTEMI", ]

###volcano sobreexpresado STEMI
volcano_plot_STEMI_OP <- ggplot(genes_OP_STEMI_volcano, aes(x = genes_OP_STEMI_volcano$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(genes_OP_STEMI_volcano$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()
print(volcano_plot_STEMI_OP)

##volcano sobreexpresado NSTEMI
volcano_plot_NSTEMI_OP <- ggplot(genes_OP_NSTEMI_volcano, aes(x = genes_OP_NSTEMI_volcano$log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 1, color = ifelse(genes_OP_NSTEMI_volcano$pvalue < ldc_threshold, "red", "black")) +
  geom_hline(yintercept = -log10(ldc_threshold), linetype = "dashed", color = "gray") +
  labs(x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw()
print(volcano_plot_NSTEMI_OP)

##bar plot

sorted_genes <- overexpressed_genes_volcano[order(overexpressed_genes_volcano$log2FoldChange, decreasing = TRUE), ]


top_genes <- head(sorted_genes, 20)
genes <- rownames(top_genes)

bar_plot <- ggplot(top_genes, aes(x = genes, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Gene", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(bar_plot)

###genes overexpressed from STEMI
genes_OP_STEMI_barplot <- sorted_genes[dds2$IAM == "STEMI", ]

top_genes_STEMI <- head(genes_OP_STEMI_barplot, 20)
genes_STEMI <- rownames(top_genes_STEMI)

bar_plot_STEMI <- ggplot(top_genes_STEMI, aes(x = genes_STEMI, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Gene", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(bar_plot_STEMI)
###genes overexpressed from NSTEMI
genes_OP_NSTEMI_barplot <- sorted_genes[dds2$IAM == "NSTEMI", ]

top_genes_NSTEMI <- head(genes_OP_NSTEMI_barplot, 20)
genes_NSTEMI <- rownames(top_genes_NSTEMI)

bar_plot_NSTEMI <- ggplot(top_genes_NSTEMI, aes(x = genes_NSTEMI, y = log2FoldChange)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(x = "Gene", y = "Log2 Fold Change") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(bar_plot_NSTEMI)

#Plot de enriquecimiento. 
##STEMI 

sorted_genes_STEMI <- enrichment_results_STEMI$ID[order(enrichment_results_STEMI$pvalue)]
top_genes_STEMI <- sorted_genes_STEMI[1:20]

###Dotplot
dotchart(enrichment_results_STEMI$pvalue[1:20], labels = top_genes_STEMI, xlim = c(0, max(enrichment_results_STEMI$pvalue)), xlab = "P-value", ylab = "Genes", main = "Top 20 Genes (Ordered by P-value) in STEMI")

###Barplot
barplot(enrichment_results_STEMI$pvalue[1:20], names.arg = top_genes_STEMI, xlab = "Genes", ylab = "P-value", main = "Top 20 Genes (Ordered by P-value) in STEMI")

##NSTEMI 

sorted_genes_NSTEMI <- enrichment_results_NSTEMI$ID[order(enrichment_results_NSTEMI$pvalue)]
top_genes_NSTEMI <- sorted_genes_NSTEMI[1:20]

###Dotplot
dotchart(enrichment_results_NSTEMI$pvalue[1:20], labels = top_genes_NSTEMI, xlim = c(0, max(enrichment_results_NSTEMI$pvalue)), xlab = "P-value", ylab = "Genes", main = "Top 20 Genes (Ordered by P-value) in NSTEMI")

###Barplot
barplot(enrichment_results_NSTEMI$pvalue[1:20], names.arg = top_genes_NSTEMI, xlab = "Genes", ylab = "P-value", main = "Top 20 Genes (Ordered by P-value) in NSTEMI")


