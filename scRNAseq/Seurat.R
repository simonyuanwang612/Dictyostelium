####basic analysis
library(Seurat)
library(dplyr)
library(cowplot)
library(harmony)
library(ggplot2)
library(pheatmap)
library(patchwork)
col = c("gray","yellow","orange","red")
setwd("/Users/plg/Dropbox (UMass Medical School)/Dicty_scRNAseq/Run3_unmapped/")
Dicty = readRDS("Dicty.rds")

Idents(Dicty) = Dicty$stage
Dicty3 = subset(Dicty, ident = "Streaming", invert = T)

df = Dicty3@meta.data
df[grep("Sl", df$orig.ident),]$stage = "Slug"
Dicty3@meta.data = df
Dicty3 = SCTransform(Dicty3)
Dicty3 <- RunPCA(Dicty3)
ElbowPlot(Dicty3, ndims = 40)
Dicty3 = RunHarmony(Dicty3, group.by.vars = "batch", plot_convergence = TRUE, assay.use="SCT", kmeans_init_nstart=20, kmeans_init_iter_max=500)
Dicty3 <- RunUMAP(Dicty3, reduction = "harmony", dims = 1:11)
Dicty3 <- FindNeighbors(Dicty3 ,reduction = "harmony", dims = 1:11) 
Dicty3 <- FindClusters(Dicty3 ,resolution = 0.5)


jpeg("3_stage_Sample_UMAP_harmony.jpeg", width = 10, height = 10, res = 300, units = "in")
DimPlot(Dicty3, group.by = "orig.ident")
dev.off()
jpeg("3_stage_Stage_UMAP_harmony.jpeg", width = 10, height = 10, res = 300, units = "in")
DimPlot(Dicty3, group.by = "stage")
dev.off()
jpeg("3_stage_Sample_UMAP_harmony_split.jpeg", width = 12, height = 15, res = 300, units = "in")
DimPlot(Dicty3, group.by = "orig.ident", split.by = "orig.ident", ncol = 3)
dev.off()
jpeg("3_stage_multicellular_genes.jpeg", width = 12, height = 12, res = 300, units = "in")
FeaturePlot(Dicty3, features = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013"), order = T, cols = col)
dev.off()

jpeg("3_stage_Cluster_UMAP.jpeg", width = 10, height = 10, res = 300, units = "in")
DimPlot(Dicty3, group.by = "seurat_clusters", label = T)+NoLegend()
dev.off()

markers <- FindAllMarkers(Dicty3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
saveRDS(markers,"3_stage_markers.rds")
write.table(markers, "3_stage_DEGs_of_24_clusters.txt", sep = "\t", col.names = NA)
top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

jpeg("3_stage_Top5_Markers_by_cluster.jpeg", width = 30, height = 5, res = 300, units = "in")
DotPlot(Dicty3, features = as.character(unique(top5$gene)), cols = "RdBu")+ theme(axis.text.x = element_text(angle = 90))
dev.off()

saveRDS(Dicty3, file = "Dicty3.rds")

####Characterize cell types
Germ = c("gerA","gerB","gerC")
Prespore = c("pspA","cotB","cotC","D7","dmtA")
Prestalk = c("ecmA","ecmB")
Cup = c("cupB","srfA","mybE","bzpH","bzpD")
Stalk = c("ecmB","bzpF","gtaG")
Spore = c("dr1","DDB-G0278179","DDB-G0274691","DDB-G0271886","gtaH","hbx10","gtaN","cdc5l","DDB-G0267638")

jpeg("3_stage_Cup.jpeg", width = 20, height = 4, units = "in", res = 300)
FeaturePlot(Dicty3, features = Cup, order = T, cols = col, ncol = 5)
dev.off()
jpeg("3_stage_Stalk.jpeg", width = 12, height = 4, units = "in", res = 300)
FeaturePlot(Dicty3, features = Stalk, order = T, cols = col, ncol = 3)
dev.off()
jpeg("3_stage_Spore.jpeg", width = 20, height = 8, units = "in", res = 300)
FeaturePlot(Dicty3, features = Spore, order = T, cols = col,ncol = 5)
dev.off()
jpeg("3_stage_Germ_cells.jpeg", width = 12, height = 4, units = "in", res = 300)
FeaturePlot(Dicty3, features = Germ, order = T, cols = col, ncol = 3)
dev.off()
jpeg("3_stage_Prespore.jpeg", width = 20, height = 4, units = "in", res = 300)
FeaturePlot(Dicty3, features = Prespore, order = T, cols = col, ncol = 5)
dev.off()
jpeg("3_stage_Prestalk.jpeg", width = 8, height = 4, units = "in", res = 300)
FeaturePlot(Dicty3, features = Prestalk, order = T, cols = col, ncol = 2)
dev.off()

####DEGs by state
Dicty3$cell_type = "Multi"
df = Dicty3@meta.data
df[df$seurat_clusters%in%c(22,12,8,18,14,4,11,6),]$cell_type = "Single"
Dicty3@meta.data = df
jpeg("3_stage_cell_type_UMAP.jpeg", width = 10, height = 10, res = 300, units = "in")
DimPlot(Dicty3, group.by = "cell_type")
dev.off()

saveRDS(Dicty3, file = "Dicty3.rds")

Idents(Dicty3) = Dicty3$cell_type
markers <- FindAllMarkers(Dicty3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(markers,"markers_by_cell_type.rds")
write.table(markers, "3_stage_DEGs_multi_single.txt", sep = "\t", col.names = NA)

#Prestalk genes
#Prestalk-A region tip orgnizer: cudA
Dicty3 = readRDS("Dicty3.rds")
FeaturePlot(Dicty3, features = "cudA", order = T, cols = col)
FeaturePlot(Dicty3, features = "ecmA", order = T, cols = col) #pstA
FeaturePlot(Dicty3, features = "ecmB", order = T, cols = col) #pstAB, pstB
FeaturePlot(Dicty3, features = "ecmF", order = T, cols = col) #pstA
FeaturePlot(Dicty3, features = "staA", order = T, cols = col) #pstA
FeaturePlot(Dicty3, features = "staB", order = T, cols = col) #prestalk basal disc
FeaturePlot(Dicty3, features = "abpC", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "DDB-G0283503", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "DDB-G0286193", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "DDB-G0283511", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "DDB-G0283501", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "tpsA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "DDB-G0283515", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "5NT", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "DDB-G0283465", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "rasC", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "DDB-G0293356", order = T, cols = col) #pstO
FeaturePlot(Dicty3, features = "DDB-G0283507", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "hssA", order = T, cols = col) #pstAO
FeaturePlot(Dicty3, features = "sigN2", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "cyrA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk

FeaturePlot(Dicty3, features = "DDB-G0269702", order = T, cols = col) #pstAB cells and in pstA cells, pstO cells, upper and lower cups, and stalk during culmination
FeaturePlot(Dicty3, features = "DDB-G0267898", order = T, cols = col) #pstAB
FeaturePlot(Dicty3, features = "expl7", order = T, cols = col) #expressed in pstAB cells and in pstAB and stalk cells during the Mexican hat stage and culmination
FeaturePlot(Dicty3, features = "DDB-G0293862", order = T, cols = col) #expressed in pstA and pstO cells and in pstAB cells during culmination
FeaturePlot(Dicty3, features = "dipA", order = T, cols = col) #expressed in pstAB cells and in upper cup during culmination
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk
FeaturePlot(Dicty3, features = "dtfA", order = T, cols = col) #prestalk

#Syntenic gene list from 
library(readxl)
Dicty3 = readRDS("Dicty3.rds")
Glist = as.data.frame(read_excel("Syntenic list update 0626.xlsx"))
gene_id_table = read.delim("/Users/plg/Dropbox (UMass Medical School)/Dicty_scRNAseq/gene_information.txt", sep = "\t")

Genelist = list()
for (i in 1:4) {
  expressed_id = as.character(gene_id_table$GENE.ID[gene_id_table$GENE.ID%in%Glist[,i]])
  expressed_id2 = as.character(gene_id_table$Gene.Name[gene_id_table$GENE.ID%in%Glist[,i]])
  expressed_id3 = as.character(gene_id_table$Gene.Name[gene_id_table$Gene.Name%in%Glist[,i]])
  expressed_id4 = as.character(gene_id_table$GENE.ID[gene_id_table$Gene.Name%in%Glist[,i]])
  tmp = unique(c(expressed_id2,expressed_id,expressed_id3,expressed_id4))
  tmp = gsub("DDB_","DDB-",tmp)
  tmp = tmp[tmp%in%rownames(Dicty3)]
  Genelist[[i]] = tmp
}

Idents(Dicty3) = Dicty3$cell_type
for (i in 1:4) {
  w = length(Genelist[[i]])*0.25+10
  jpeg(paste0("3_stage_genelist_",i,"_single_multi_dotplot.jpeg"), width = w, height = 5, res = 300, units = "in")
  print(DotPlot(Dicty3, features = Genelist[[i]], cols = "RdBu")+ theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(colnames(Glist)[i])))
  dev.off()
}

Idents(Dicty3) = Dicty3$stage
for (i in 1:4) {
  w = length(Genelist[[i]])*0.25+10
  jpeg(paste0("3_stage_genelist_",i,"_stage_dotplot.jpeg"), width = w, height = 5, res = 300, units = "in")
  print(DotPlot(Dicty3, features = Genelist[[i]], cols = "RdBu")+ theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(colnames(Glist)[i])))
  dev.off()
}

Idents(Dicty3) = Dicty3$seurat_clusters
for (i in 1:4) {
  w = length(Genelist[[i]])*0.25+10
  jpeg(paste0("3_stage_genelist_",i,"_cluster_dotplot.jpeg"), width = w, height = 5, res = 300, units = "in")
  print(DotPlot(Dicty3, features = Genelist[[i]], cols = "RdBu")+ theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(colnames(Glist)[i])))
  dev.off()
}

#Figure 3E
jpeg("Figure3E.jpeg", height = 5, width = 5, units = "in", res = 300)
DimPlot(Dicty3, group.by = "stage", cols = c("#F69324","#F4EB22","#498AC9"), pt.size = 0.05)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_line(colour = 'black', size = 1),
        axis.ticks=element_line(colour = 'black', size = 1))+NoLegend()+NoAxes()
dev.off()
jpeg("Figure3E_Original.jpeg", height = 10, width = 10, units = "in", res = 300)
DimPlot(Dicty3, group.by = "stage", cols = c("#F69324","#F4EB22","#498AC9"))
dev.off()
#Figure 3F
library(extrafont)
font_import()
i=2
pdf("Figure3F1.pdf", height = 5, width = 20)
DotPlot(Dicty3, group.by = "stage",features = Genelist[[i]], cols = "RdBu")+
  theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(colnames(Glist)[i]))+
  scale_size(range = c(2, 10))
dev.off()
jpeg("Figure3F1.jpeg", height = 2, width = 17, units = "in",res = 300)
DotPlot(Dicty3, group.by = "stage",features = Genelist[[i]], cols = "RdBu")+
  scale_size(range = c(2, 10))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_line(colour = 'black', size = 1),
        axis.ticks=element_line(colour = 'black', size = 1))+NoLegend()
dev.off()

jpeg("Figure3F1-2.jpeg", height = 3, width = 17, units = "in",res = 300)
DotPlot(Dicty3, group.by = "stage",features = Genelist[[i]], cols = "RdBu")+
  scale_size(range = c(2, 10))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_line(colour = 'black', size = 1),
        axis.ticks=element_line(colour = 'black', size = 1))+ theme(legend.position="top")
dev.off()

i=4
pdf("Figure3F2.pdf", height = 5, width = 20)
DotPlot(Dicty3, group.by = "stage",features = Genelist[[i]], cols = "RdBu")+
  theme(axis.text.x = element_text(angle = 90))+ggtitle(paste0(colnames(Glist)[i]))+
  scale_size(range = c(2, 10))
dev.off()
jpeg("Figure3F2.jpeg", height = 2, width = 17, units = "in",res = 300)
DotPlot(Dicty3, group.by = "stage",features = Genelist[[i]], cols = "RdBu")+
  scale_size(range = c(2, 10))+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_line(colour = 'black', size = 1),
        axis.ticks=element_line(colour = 'black', size = 1))+NoLegend()
dev.off()


#Figure 3H
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
plots = list()
for (i in 1:4) {
  plots[[i]] = FeaturePlot(Dicty3, features = genes[i], order = T, cols = col, pt.size = 1)+NoAxes()+NoLegend()+ggtitle(element_blank())

}
jpeg("Figure3H.jpeg", width = 20, height = 5, res = 300, units = "in")
wrap_plots(plots, ncol = 4)
dev.off()
jpeg("Figure3H_original.jpeg", width = 20, height = 5, res = 300, units = "in")
FeaturePlot(Dicty3, features = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013"), order = T, cols = col, ncol = 4)
dev.off()

###DDB_G0293674
jpeg("DDB_G0293674.jpeg", width = 5, height = 5, res = 300, units = "in")
FeaturePlot(Dicty3, features = "DDB-G0293674", order = T, cols = col)+NoAxes()+NoLegend()+ggtitle(element_blank())
dev.off()

#for GEO submission
write.table(Dicty@assays$RNA@counts, file = "raw_sc_Dicty_4_stage.tsv", sep = "\t", col.names = NA)
write.table(Dicty@assays$SCT@data, file = "normalized_sc_Dicty_4_stage.tsv", sep = "\t", col.names = NA)

df = Dicty@meta.data
df2 = df[,1:7]
UMAP = Dicty@reductions$umap@cell.embeddings
df2 = merge(df, UMAP, by = 0, all = T)  
harmony = Dicty@reductions$harmony@cell.embeddings
rownames(df2) = df2[,1]
df2 = df2[,-1]
df3 = merge(df2, harmony, by = 0, all = T)  
rownames(df3) = df3[,1]
df3 = df3[,-1]
write.table(df3, file = "Metadata_Dicty_4_stage.tsv", sep = "\t", col.names = NA)

write.table(Dicty3@assays$RNA@counts, file = "raw_sc_Dicty_3_stage.tsv", sep = "\t", col.names = NA)
write.table(Dicty3@assays$SCT@data, file = "normalized_sc_Dicty_3_stage.tsv", sep = "\t", col.names = NA)

df = Dicty3@meta.data
df2 = df[,c(1:7,11)]
UMAP = Dicty3@reductions$umap@cell.embeddings
df2 = merge(df, UMAP, by = 0, all = T)  
harmony = Dicty3@reductions$harmony@cell.embeddings
rownames(df2) = df2[,1]
df2 = df2[,-1]
df3 = merge(df2, harmony, by = 0, all = T)
rownames(df3) = df3[,1]
df3 = df3[,-1]
write.table(df3, file = "Metadata_Dicty_3_stage.tsv", sep = "\t", col.names = NA)

