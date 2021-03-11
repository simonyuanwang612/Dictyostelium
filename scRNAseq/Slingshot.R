library(destiny)
library(slingshot)
library(SingleCellExperiment)
library(ggplot2)
library(mclust)
library(Seurat)
library(plotly)
library(ggpubr)
library(RColorBrewer)
library(rgl)
library(emdbook)
require(gam)
library(pheatmap)
library(harmony)
Dicty3 = readRDS("Dicty3.rds")
metadf = Dicty3@meta.data
mat = GetAssayData(Dicty3,slot = "counts")

#take only subset of cells for slingshot 
cells = sample(1:ncol(Dicty3), round(ncol(Dicty3)/4, 0))
metadf = metadf[cells,]
mat = mat[,cells]
Dicty3_small = CreateSeuratObject(counts = mat,meta.data = metadf)
Dicty3_small = SCTransform(Dicty3_small)
Dicty3_small <- RunPCA(Dicty3_small)
Dicty3_small = RunHarmony(Dicty3_small, group.by.vars = "batch", plot_convergence = TRUE, assay.use="SCT")
ElbowPlot(Dicty3_small, reduction = "harmony")
Dicty3_small <- RunUMAP(Dicty3_small, reduction = "harmony", dims = 1:11)

Dicty3_small <- FindNeighbors(Dicty3_small ,reduction = "harmony", dims = 1:11) 
Dicty3_small <- FindClusters(Dicty3_small ,resolution = 0.5)
DimPlot(Dicty3_small, group.by = "seurat_clusters", label = T)
DimPlot(Dicty3_small, group.by = "stage")
FeaturePlot(Dicty3_small, features = "DDB-G0288579", order = T, cols = col)
saveRDS(Dicty3_small, file = "Dicty3_small.rds")

sim <- SingleCellExperiment(assays = List(counts = mat))
#filtering
geneFilter <- apply(assays(sim)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sim <- sim[geneFilter, ]

#normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)

#PCA
pca <- prcomp(t(log1p(assays(sim)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:3]
saveRDS(pca, "pca_for_slingshot.rds")
colData(sim)$stage = metadf[colnames(sim),]$stage
colData(sim)$pc1 = pca$x[,1]
colData(sim)$pc2 = pca$x[,2]
colData(sim)$pc3 = pca$x[,3]
coldatat = as.data.frame(colData(sim))
p1 = ggplot(coldatat, aes(x = pc1, y= pc2, color = stage))+ geom_point(size =.5, alpha = .5)
p2 = ggplot(coldatat, aes(x = pc2, y= pc3, color = stage))+ geom_point(size =.5, alpha = .5)
p3 = ggplot(coldatat, aes(x = pc1, y= pc3, color = stage))+ geom_point(size =.5, alpha = .5)
jpeg("PCA_stages.jpeg", width = 15, height = 5, res = 300, units = "in")
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
plot_ly(coldatat, x=~pc1, y=~pc2, z=~pc3, color = ~stage, marker=list(size=2))

#DM
dm <- DiffusionMap(t(log1p(assays(sim)$norm)), n_pcs = 50)
saveRDS(dm, "DM_for_slingshot.rds")
colData(sim)$DC1 = dm$DC1
colData(sim)$DC2 = dm$DC2
colData(sim)$DC3 = dm$DC3
coldatat = as.data.frame(colData(sim))
p1 = ggplot(coldatat, aes(x = DC1, y= DC2, color = stage))+ geom_point(size =.5, alpha = .5)
p2 = ggplot(coldatat, aes(x = DC2, y= DC3, color = stage))+ geom_point(size =.5, alpha = .5)
p3 = ggplot(coldatat, aes(x = DC1, y= DC3, color = stage))+ geom_point(size =.5, alpha = .5)
jpeg("DM_stages.jpeg", width = 15, height = 5, res = 300, units = "in")
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()
plot_ly(coldatat, x=~DC1, y=~DC2, z=~DC3, color = ~stage, marker=list(size=2))
rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2, DC3 = dm$DC3)

reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)

#cluster cells

cl1 <- Mclust(rd1,G = 1:15)$classification
colData(sim)$PCA_GMM <- as.factor(cl1)
ggplot(as.data.frame(colData(sim)), aes(x = pc1, y= pc2, color = PCA_GMM))+ geom_point(size =.5, alpha = .5)
plot_ly(as.data.frame(colData(sim)), x=~pc1, y=~pc2, z=~pc3, color = ~PCA_GMM, marker=list(size=2))

cl2 <- kmeans(rd1, centers = 15)$cluster
colData(sim)$PCA_kmeans <- as.factor(cl2)
ggplot(as.data.frame(colData(sim)), aes(x = pc1, y= pc2, color = PCA_kmeans))+ geom_point(size =.5, alpha = .5)
plot_ly(as.data.frame(colData(sim)), x=~pc1, y=~pc2, z=~pc3, color = ~PCA_kmeans, marker=list(size=2))

cl3 <- Mclust(rd2,G = 1:15)$classification
colData(sim)$DM_GMM <- as.factor(cl3)
ggplot(as.data.frame(colData(sim)), aes(x = DC1, y= DC2, color = DM_GMM))+ geom_point(size =.5, alpha = .5)
plot_ly(as.data.frame(colData(sim)), x=~DC1, y=~DC2, z=~DC3, color = ~DM_GMM, marker=list(size=2))

cl4 <- kmeans(as.data.frame(rd2), centers = 15)$cluster
colData(sim)$DM_kmeans <- as.factor(cl4)
ggplot(as.data.frame(colData(sim)), aes(x = DC1, y= DC2, color = DM_kmeans))+ geom_point(size =.5, alpha = .5)
plot_ly(as.data.frame(colData(sim)), x=~DC1, y=~DC2, z=~DC3, color = ~DM_kmeans, marker=list(size=2))


saveRDS(sim, "sim_before_slingshot.rds")


# check the root and tip for each setting
col = c("gray","yellow","orange","red")
coldata = as.data.frame(colData(sim))
rd1 = rd1[colnames(Dicty3_small),]
rd2 = rd2[colnames(Dicty3_small),]
Dicty3_small[["Slingshot_PCA"]] <- CreateDimReducObject(embeddings = rd1, key = "PCA_", assay = DefaultAssay(Dicty3_small))
Dicty3_small[["Slingshot_DM"]] <- CreateDimReducObject(embeddings = rd2, key = "DC_", assay = DefaultAssay(Dicty3_small))
coldata = coldata[,-(1:7)]
coldata = coldata[colnames(Dicty3_small),]
metadf = Dicty3_small@meta.data
metadf = cbind.data.frame(metadf, coldata)
Dicty3_small@meta.data = metadf
FeaturePlot(Dicty3_small, "DDB-G0288579", cols = col, order = T)
#root
DimPlot(Dicty3_small,group.by = "PCA_GMM", label = T) # 7
DimPlot(Dicty3_small,group.by = "PCA_kmeans", label = T) #6
DimPlot(Dicty3_small,group.by = "DM_GMM", label = T) #6
DimPlot(Dicty3_small,group.by = "DM_kmeans", label = T) #14
#tips
DimPlot(Dicty3_small,group.by = "PCA_GMM", label = T) #11,13,5,2,15
DimPlot(Dicty3_small,group.by = "PCA_kmeans", label = T) #4,14,9,11,7,12
DimPlot(Dicty3_small,group.by = "DM_GMM", label = T) #8,9,10,7,15
DimPlot(Dicty3_small,group.by = "DM_kmeans", label = T) #7,9,6,13,11



#######The following part are examination of the different dimension reduction and clustering for slingshot, could skip########
#Slingshot_Diffusion_map_with_GMM_clustering
sim_DG <- slingshot(sim, clusterLabels = 'DM_GMM', reducedDim = 'DiffMap', start.clus = 6, end.clus = c(8,9,10,7,15))
saveRDS(sim_DG, file = "DM_GMM_sim.rds")
coldata2 = as.data.frame(colData(sim_DG))

mycolors <- c("#F8766D", "#7CAE00","#00BFC4")
stages = c("Fruiting","Slug","Vegetative")
dirnames = "DM_GMM"####
Angle1 <- 5 
Angle <- rep(Angle1 * pi / 180, 360/Angle1) 
Yall <- log1p(assays(sim_DG)$norm)####

myplot = function(df, x_string, y_string) {
  ggplot(df, aes_string(x = x_string, y = y_string))
}
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
genes = genes[genes%in%rownames(Yall)]
var1000 <- names(sort(apply(Yall,1,var),decreasing = TRUE))[1:1000]
Y <- Yall[var1000,]
lineages = grep("^sling", colnames(coldata2), value = T)
curves = slingCurves(sim_DG)####
for (i in 1:length(lineages)) {
  curve = as.data.frame(curves[[i]]["s"])
  #colnames(curves) = gsub("^","pt_",colnames(curves))
  coldata3 = cbind.data.frame(coldata2, curve)
  #plot_ly() %>% 
  #  add_trace(data = coldata3, x=~DC1, y=~DC2, z=~DC3, color = ~stage, mode = "markers", marker=list(size=2)) %>% 
  #  add_trace(data = coldata3, x = ~s.DC1, y = ~s.DC2, z = ~s.DC3, name = 'Pseudotime_1', marker=list(size=2))
  coldata3$color = NA
  for (a in 1:length(stages)) {
    coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
  }
  par3d("windowRect"=4*(par3d("windowRect")))
  rgl.viewpoint(  zoom = .7 )
  plot3d(x= coldata3$s.DC1, y =coldata3$s.DC2, z = coldata3$s.DC3, size = 5, col="navy",
         xlab = "DC1", ylab = "DC3", zlab = "DC3", xlim = c(min(coldata3$DC1),max(coldata3$DC1)), ylim = c(min(coldata3$DC2),max(coldata3$DC2)), zlim = c(min(coldata3$DC3),max(coldata3$DC3)))
  points3d( coldata3$DC1, coldata3$DC2, coldata3$DC3, size=4, axes=T, radius = .1,
            colkey = F, col = coldata3$color, alpha =.5)
  legend3d("topright", legend = stages, pch = 16, col = mycolors, cex=1, inset=c(0))
  dir.create(paste0("~/Desktop/",dirnames,"_",lineages[i]))
  Animation.dir = paste0("~/Desktop/",dirnames,"_",lineages[i])
  for (b in seq(Angle)) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"),Angle[b], 0, 1, 0))
    rgl.snapshot(filename=paste(paste(Animation.dir, "/frame-", sep=""),sprintf("%03d", b), ".png", sep="")) }
  rgl.close() 
  move.to.animation.dir <- paste("cd", Animation.dir, "&&")
  ImageMagick.code <- paste("convert -delay ", 3, " -loop 0 frame*.png animated.gif", sep="") 
  system(paste(move.to.animation.dir, ImageMagick.code))
  ImageMagick.code2 = "convert animated.gif -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 patrol_cycle.gif"
  system(paste(move.to.animation.dir, ImageMagick.code2))
  ImageMagick.code3 = "convert -coalesce patrol_cycle.gif some%05d.png"
  system(paste(move.to.animation.dir, ImageMagick.code3))
  ImageMagick.code4 = paste0("ffmpeg -i some%05d.png ",lineages[i],".mov")
  system(paste(move.to.animation.dir, ImageMagick.code4))
  ImageMagick.code5 = paste0("ffmpeg -i ",lineages[i],".mov"," -acodec copy -vcodec copy ",lineages[i],".mp4")
  system(paste(move.to.animation.dir, ImageMagick.code5))
  ImageMagick.code6 = "rm some*.png"
  system(paste(move.to.animation.dir, ImageMagick.code6))
  t <- coldata3[,lineages[i]]
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
  })
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
  topgenes <- c(genes,topgenes)
  heatdata <- Yall[topgenes, order(t, na.last = NA)]
  heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
  
  #heatmap(log1p(heatdata), Colv = NA,
  #        ColSideColors = brewer.pal(9,"Set1")[heatclus])
  ann = coldata3[,c("stage",lineages[i],dirnames)]
  jpeg(paste(dirnames,lineages[i],"gene_expression.jpeg",sep = "_"), h = 20, w =12, unit = "in", res = 300)
  print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu"))))
  dev.off()
  
  
  
  jpeg(paste(dirnames,lineages[i],"cell_order.jpeg", sep = "_"), h = 5, w =12, unit = "in", res = 300)
  print(myplot(coldata3, paste0(lineages[i]),paste0("stage")) + 
          geom_point(aes(alpha=0.01, color = stage), size = 10) + 
          theme_minimal())
  dev.off()
}





#Slingshot_Diffusion_map_with_Kmeans_clustering
sim_DK <- slingshot(sim, clusterLabels = 'DM_kmeans', reducedDim = 'DiffMap', start.clus = 14, end.clus = c(7,9,6,13,11))
saveRDS(sim_DK, file = "DM_kmeans_sim.rds")

coldata2 = as.data.frame(colData(sim_DK))########

mycolors <- c("#F8766D", "#7CAE00","#C77CFF")
stages = c("Fruiting","Slug","Vegetative")
dirnames = "DM_kmeans"#######
Angle1 <- 5 
Angle <- rep(Angle1 * pi / 180, 360/Angle1) 
Yall <- log1p(assays(sim_DK)$norm)#######

myplot = function(df, x_string, y_string) {
  ggplot(df, aes_string(x = x_string, y = y_string))
}
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
genes = genes[genes%in%rownames(Yall)]
var1000 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Yall[var1000,]
lineages = grep("^sling", colnames(coldata2), value = T)
curves = slingCurves(sim_DK)###
for (i in 1:length(lineages)) {
  curve = as.data.frame(curves[[i]]["s"])
  #colnames(curves) = gsub("^","pt_",colnames(curves))
  coldata3 = cbind.data.frame(coldata2, curve)
  #plot_ly() %>% 
  #  add_trace(data = coldata3, x=~DC1, y=~DC2, z=~DC3, color = ~stage, mode = "markers", marker=list(size=2)) %>% 
  #  add_trace(data = coldata3, x = ~s.DC1, y = ~s.DC2, z = ~s.DC3, name = 'Pseudotime_1', marker=list(size=2))
  coldata3$color = NA
  for (a in 1:length(stages)) {
    coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
  }
  par3d("windowRect"=4*(par3d("windowRect")))
  rgl.viewpoint(  zoom = .7 )
  plot3d(x= coldata3$s.DC1, y =coldata3$s.DC2, z = coldata3$s.DC3, size = 5, col="navy",
         xlab = "DC1", ylab = "DC3", zlab = "DC3", xlim = c(min(coldata3$DC1),max(coldata3$DC1)), ylim = c(min(coldata3$DC2),max(coldata3$DC2)), zlim = c(min(coldata3$DC3),max(coldata3$DC3)))
  points3d( coldata3$DC1, coldata3$DC2, coldata3$DC3, size=4, axes=T, radius = .1,
            colkey = F, col = coldata3$color, alpha =.5)
  legend3d("topright", legend = stages, pch = 16, col = mycolors, cex=1, inset=c(0))
  dir.create(paste0("~/Desktop/",dirnames,"_",lineages[i]))
  Animation.dir = paste0("~/Desktop/",dirnames,"_",lineages[i])
  for (b in seq(Angle)) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"),Angle[b], 0, 1, 0))
    rgl.snapshot(filename=paste(paste(Animation.dir, "/frame-", sep=""),sprintf("%03d", b), ".png", sep="")) }
  rgl.close() 
  move.to.animation.dir <- paste("cd", Animation.dir, "&&")
  ImageMagick.code <- paste("convert -delay ", 3, " -loop 0 frame*.png animated.gif", sep="") 
  system(paste(move.to.animation.dir, ImageMagick.code))
  ImageMagick.code2 = "convert animated.gif -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 patrol_cycle.gif"
  system(paste(move.to.animation.dir, ImageMagick.code2))
  ImageMagick.code3 = "convert -coalesce patrol_cycle.gif some%05d.png"
  system(paste(move.to.animation.dir, ImageMagick.code3))
  ImageMagick.code4 = paste0("ffmpeg -i some%05d.png ",lineages[i],".mov")
  system(paste(move.to.animation.dir, ImageMagick.code4))
  ImageMagick.code5 = paste0("ffmpeg -i ",lineages[i],".mov"," -acodec copy -vcodec copy ",lineages[i],".mp4")
  system(paste(move.to.animation.dir, ImageMagick.code5))
  ImageMagick.code6 = "rm some*.png"
  system(paste(move.to.animation.dir, ImageMagick.code6))
  t <- coldata3[,lineages[i]]
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
  })
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
  topgenes <- c(genes,topgenes)
  heatdata <- Yall[topgenes, order(t, na.last = NA)]
  heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
  
  #heatmap(log1p(heatdata), Colv = NA,
  #        ColSideColors = brewer.pal(9,"Set1")[heatclus])
  ann = coldata3[,c("stage",lineages[i],dirnames)]
  jpeg(paste(dirnames,lineages[i],"gene_expression.jpeg",sep = "_"), h = 20, w =12, unit = "in", res = 300)
  print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu"))))
  dev.off()
  
  
  
  jpeg(paste(dirnames,lineages[i],"cell_order.jpeg", sep = "_"), h = 5, w =12, unit = "in", res = 300)
  print(myplot(coldata3, paste0(lineages[i]),paste0("stage")) + 
          geom_point(aes(alpha=0.01, color = stage), size = 10) + 
          theme_minimal())
  dev.off()
}








#Slingshot_PCA_with_GMM_clustering
sim_PG <- slingshot(sim, clusterLabels = 'PCA_GMM', reducedDim = 'PCA', start.clus = 7, end.clus = c(11,13,5,2,15))
saveRDS(sim_PG, file = "PCA_GMM_sim.rds")

coldata2 = as.data.frame(colData(sim_PG))########

mycolors <- c("#F8766D", "#7CAE00","#C77CFF")
stages = c("Fruiting","Slug","Vegetative")
dirnames = "PCA_GMM"#######
Angle1 <- 5 
Angle <- rep(Angle1 * pi / 180, 360/Angle1) 
Yall <- log1p(assays(sim_PG)$norm)#######
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
genes = genes[genes%in%rownames(Yall)]
myplot = function(df, x_string, y_string) {
  ggplot(df, aes_string(x = x_string, y = y_string))
}

var1000 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Yall[var1000,]
lineages = grep("^sling", colnames(coldata2), value = T)
curves = slingCurves(sim_PG)###

for (i in 1:length(lineages)) {
  curve = as.data.frame(curves[[i]]["s"])
  #colnames(curves) = gsub("^","pt_",colnames(curves))
  coldata3 = cbind.data.frame(coldata2, curve)
  #plot_ly() %>% 
  #  add_trace(data = coldata3, x=~DC1, y=~DC2, z=~DC3, color = ~stage, mode = "markers", marker=list(size=2)) %>% 
  #  add_trace(data = coldata3, x = ~s.DC1, y = ~s.DC2, z = ~s.DC3, name = 'Pseudotime_1', marker=list(size=2))
  coldata3$color = NA
  for (a in 1:length(stages)) {
    coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
  }
  par3d("windowRect"=4*(par3d("windowRect")))
  rgl.viewpoint(  zoom = .7 )
  plot3d(x= coldata3$s.PC1, y =coldata3$s.PC2, z = coldata3$s.PC3, size = 5, col="navy",
         xlab = "PC1", ylab = "PC2", zlab = "PC3", xlim = c(min(coldata3$pc1),max(coldata3$pc1)), ylim = c(min(coldata3$pc2),max(coldata3$pc2)), zlim = c(min(coldata3$pc3),max(coldata3$pc3)))
  points3d( coldata3$pc1, coldata3$pc2, coldata3$pc3, size=4, axes=T, radius = .1,
            colkey = F, col = coldata3$color, alpha =.5)
  legend3d("topright", legend = stages, pch = 16, col = mycolors, cex=1, inset=c(0))
  dir.create(paste0("~/Desktop/",dirnames,"_",lineages[i]))
  Animation.dir = paste0("~/Desktop/",dirnames,"_",lineages[i])
  for (b in seq(Angle)) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"),Angle[b], 0, 1, 0))
    rgl.snapshot(filename=paste(paste(Animation.dir, "/frame-", sep=""),sprintf("%03d", b), ".png", sep="")) }
  rgl.close() 
  move.to.animation.dir <- paste("cd", Animation.dir, "&&")
  ImageMagick.code <- paste("convert -delay ", 3, " -loop 0 frame*.png animated.gif", sep="") 
  system(paste(move.to.animation.dir, ImageMagick.code))
  ImageMagick.code2 = "convert animated.gif -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 patrol_cycle.gif"
  system(paste(move.to.animation.dir, ImageMagick.code2))
  ImageMagick.code3 = "convert -coalesce patrol_cycle.gif some%05d.png"
  system(paste(move.to.animation.dir, ImageMagick.code3))
  ImageMagick.code4 = paste0("ffmpeg -i some%05d.png ",lineages[i],".mov")
  system(paste(move.to.animation.dir, ImageMagick.code4))
  ImageMagick.code5 = paste0("ffmpeg -i ",lineages[i],".mov"," -acodec copy -vcodec copy ",lineages[i],".mp4")
  system(paste(move.to.animation.dir, ImageMagick.code5))
  ImageMagick.code6 = "rm some*.png"
  system(paste(move.to.animation.dir, ImageMagick.code6))
  t <- coldata3[,lineages[i]]
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
  })
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
  topgenes <- c(genes,topgenes)
  heatdata <- Yall[topgenes, order(t, na.last = NA)]
  heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
  
  #heatmap(log1p(heatdata), Colv = NA,
  #        ColSideColors = brewer.pal(9,"Set1")[heatclus])
  ann = coldata3[,c("stage",lineages[i],dirnames)]
  jpeg(paste(dirnames,lineages[i],"gene_expression.jpeg",sep = "_"), h = 20, w =12, unit = "in", res = 300)
  print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu"))))
  dev.off()
  
  
  
  jpeg(paste(dirnames,lineages[i],"cell_order.jpeg", sep = "_"), h = 5, w =12, unit = "in", res = 300)
  print(myplot(coldata3, paste0(lineages[i]),paste0("stage")) + 
          geom_point(aes(alpha=0.01, color = stage), size = 10) + 
          theme_minimal())
  dev.off()
}






#Slingshot_PCA_with_kmeans_clustering
sim_PK <- slingshot(sim, clusterLabels = 'PCA_kmeans', reducedDim = 'PCA', start.clus = 6, end.clus = c(4,14,9,11,7,12))
saveRDS(sim_PK, file = "PCA_kmeans_sim.rds")

coldata2 = as.data.frame(colData(sim_PK))########

mycolors <- c("#F8766D", "#7CAE00","#C77CFF")
stages = c("Fruiting","Slug","Vegetative")
dirnames = "PCA_kmeans"#######
Angle1 <- 5 
Angle <- rep(Angle1 * pi / 180, 360/Angle1) 
Yall <- log1p(assays(sim_PK)$norm)#######
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
genes = genes[genes%in%rownames(Yall)]
myplot = function(df, x_string, y_string) {
  ggplot(df, aes_string(x = x_string, y = y_string))
}

var1000 <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Yall[var1000,]
lineages = grep("^sling", colnames(coldata2), value = T)
curves = slingCurves(sim_PK)###

for (i in 1:length(lineages)) {
  curve = as.data.frame(curves[[i]]["s"])
  #colnames(curves) = gsub("^","pt_",colnames(curves))
  coldata3 = cbind.data.frame(coldata2, curve)
  #plot_ly() %>% 
  #  add_trace(data = coldata3, x=~DC1, y=~DC2, z=~DC3, color = ~stage, mode = "markers", marker=list(size=2)) %>% 
  #  add_trace(data = coldata3, x = ~s.DC1, y = ~s.DC2, z = ~s.DC3, name = 'Pseudotime_1', marker=list(size=2))
  coldata3$color = NA
  for (a in 1:length(stages)) {
    coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
  }
  par3d("windowRect"=4*(par3d("windowRect")))
  rgl.viewpoint(  zoom = .7 )
  plot3d(x= coldata3$s.PC1, y =coldata3$s.PC2, z = coldata3$s.PC3, size = 5, col="navy",
         xlab = "PC1", ylab = "PC2", zlab = "PC3", xlim = c(min(coldata3$pc1),max(coldata3$pc1)), ylim = c(min(coldata3$pc2),max(coldata3$pc2)), zlim = c(min(coldata3$pc3),max(coldata3$pc3)))
  points3d( coldata3$pc1, coldata3$pc2, coldata3$pc3, size=4, axes=T, radius = .1,
            colkey = F, col = coldata3$color, alpha =.5)
  legend3d("topright", legend = stages, pch = 16, col = mycolors, cex=1, inset=c(0))
  dir.create(paste0("~/Desktop/",dirnames,"_",lineages[i]))
  Animation.dir = paste0("~/Desktop/",dirnames,"_",lineages[i])
  for (b in seq(Angle)) {
    view3d(userMatrix = rotate3d(par3d("userMatrix"),Angle[b], 0, 1, 0))
    rgl.snapshot(filename=paste(paste(Animation.dir, "/frame-", sep=""),sprintf("%03d", b), ".png", sep="")) }
  rgl.close() 
  move.to.animation.dir <- paste("cd", Animation.dir, "&&")
  ImageMagick.code <- paste("convert -delay ", 3, " -loop 0 frame*.png animated.gif", sep="") 
  system(paste(move.to.animation.dir, ImageMagick.code))
  ImageMagick.code2 = "convert animated.gif -coalesce -duplicate 1,-2-1 -quiet -layers OptimizePlus -loop 0 patrol_cycle.gif"
  system(paste(move.to.animation.dir, ImageMagick.code2))
  ImageMagick.code3 = "convert -coalesce patrol_cycle.gif some%05d.png"
  system(paste(move.to.animation.dir, ImageMagick.code3))
  ImageMagick.code4 = paste0("ffmpeg -i some%05d.png ",lineages[i],".mov")
  system(paste(move.to.animation.dir, ImageMagick.code4))
  ImageMagick.code5 = paste0("ffmpeg -i ",lineages[i],".mov"," -acodec copy -vcodec copy ",lineages[i],".mp4")
  system(paste(move.to.animation.dir, ImageMagick.code5))
  ImageMagick.code6 = "rm some*.png"
  system(paste(move.to.animation.dir, ImageMagick.code6))
  t <- coldata3[,lineages[i]]
  gam.pval <- apply(Y,1,function(z){
    d <- data.frame(z=z, t=t)
    suppressWarnings({
      tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
    })
    p <- summary(tmp)[3][[1]][2,3]
    p
  })
  topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
  topgenes <- c(genes,topgenes)
  heatdata <- Yall[topgenes, order(t, na.last = NA)]
  heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
  
  #heatmap(log1p(heatdata), Colv = NA,
  #        ColSideColors = brewer.pal(9,"Set1")[heatclus])
  ann = coldata3[,c("stage",lineages[i],dirnames)]
  jpeg(paste(dirnames,lineages[i],"gene_expression.jpeg",sep = "_"), h = 20, w =12, unit = "in", res = 300)
  print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu"))))
  dev.off()
  
  
  
  jpeg(paste(dirnames,lineages[i],"cell_order.jpeg", sep = "_"), h = 5, w =12, unit = "in", res = 300)
  print(myplot(coldata3, paste0(lineages[i]),paste0("stage")) + 
          geom_point(aes(alpha=0.01, color = stage), size = 10) + 
          theme_minimal())
  dev.off()
}
####### end of evaluation #####


###PCA_Kmeans_lineage_6
#Slingshot_PCA_with_kmeans_clustering
library(readxl)
Glist = as.data.frame(read_excel("Syntenic list update 0623.xlsx"))
gene_id_table = read.delim("gene_information.txt", sep = "\t")
sim_PK = readRDS("PCA_kmeans_sim.rds")
Yall <- log1p(assays(sim_PK)$norm)#######
Genelist = list()
for (i in 1:4) {
  expressed_id = as.character(gene_id_table$GENE.ID[gene_id_table$GENE.ID%in%Glist[,i]])
  expressed_id2 = as.character(gene_id_table$Gene.Name[gene_id_table$GENE.ID%in%Glist[,i]])
  expressed_id3 = as.character(gene_id_table$Gene.Name[gene_id_table$Gene.Name%in%Glist[,i]])
  expressed_id4 = as.character(gene_id_table$GENE.ID[gene_id_table$Gene.Name%in%Glist[,i]])
  tmp = unique(c(expressed_id2,expressed_id,expressed_id3,expressed_id4))
  tmp = gsub("DDB_","DDB-",tmp)
  tmp = tmp[tmp%in%rownames(Yall)]
  Genelist[[i]] = tmp
}


coldata2 = as.data.frame(colData(sim_PK))

mycolors <- c("#F8766D", "#7CAE00","#C77CFF")
stages = c("Fruiting","Slug","Vegetative")
dirnames = "PCA_kmeans"


lineages = grep("^sling", colnames(coldata2), value = T)
curves = slingCurves(sim_PK)

for (i in 6) {
  curve = as.data.frame(curves[[i]]["s"])
  #colnames(curves) = gsub("^","pt_",colnames(curves))
  coldata3 = cbind.data.frame(coldata2, curve)
  #plot_ly() %>% 
  #  add_trace(data = coldata3, x=~DC1, y=~DC2, z=~DC3, color = ~stage, mode = "markers", marker=list(size=2)) %>% 
  #  add_trace(data = coldata3, x = ~s.DC1, y = ~s.DC2, z = ~s.DC3, name = 'Pseudotime_1', marker=list(size=2))
  coldata3$color = NA
  for (a in 1:length(stages)) {
    coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
  }
  t <- coldata3[,lineages[i]]
  for (x in 1:4) {
    topgenes <- Genelist[[x]]
    topgenes = topgenes[topgenes%in%rownames(Yall)]
    heatdata <- Yall[topgenes, order(t, na.last = NA)]
    heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
    
    #heatmap(log1p(heatdata), Colv = NA,
    #        ColSideColors = brewer.pal(9,"Set1")[heatclus])
    ann = coldata3[,c("stage",lineages[i],dirnames)]
    h = length(topgenes)*0.15+5
    jpeg(paste(dirnames,lineages[i],"Gene_list",x,"gene_expression.jpeg",sep = "_"), h = h, w =12, unit = "in", res = 300)
    print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = brewer.pal(9,"YlOrRd"),main = colnames(Glist)[x]))
    dev.off()
  }
}

#Figure 3G
#PCA_Kmeans_lineage_6
#Slingshot_PCA_with_kmeans_clustering
sim_PK = readRDS("PCA_kmeans_sim.rds")
Dicty3_small = readRDS("Dicty3_small.rds")
mat = Dicty3_small@assays$RNA@counts
sim <- SingleCellExperiment(assays = List(counts = mat))
#normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sim)$norm <- FQnorm(assays(sim)$counts)
mycolors <- c("#F8766D", "#7CAE00","#C77CFF")
stages = c("Fruiting","Slug","Vegetative")
Yall <- log1p(assays(sim)$norm)#######
var1000 <- names(sort(apply(Yall,1,var),decreasing = TRUE))[1:1000]
Y <- Yall[var1000,]
genes = c("srfA","DDB-G0279267","DDB-G0283437","DDB-G0288013")
genes%in%rownames(sim)
coldata2 = as.data.frame(colData(sim_PK))########
dirnames = "PCA_kmeans"#######
curves = slingCurves(sim_PK)###
curve = as.data.frame(curves$curve6$s)
coldata3 = cbind.data.frame(coldata2, curve)
coldata3$color = NA
for (a in 1:length(stages)) {
  coldata3[coldata3$stage==stages[a],]$color=mycolors[a]
}
t <- coldata3[,lineages[6]]
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  suppressWarnings({
    tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
  })
  p <- summary(tmp)[3][[1]][2,3]
  p
})
saveRDS(gam.pval, file = "genes_for_3G.rds")
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
topgenes <- c(genes,topgenes)
heatdata <- Yall[topgenes, order(t, na.last = NA)]
heatclus <- coldata3[,dirnames][order(t, na.last = NA)]
ann = coldata3[,c("stage",lineages[6],dirnames)]
aka3 = list(stage = c(Fruiting = "#F69324",Slug = "#F4EB22",Vegetative = "#498AC9"))
library(viridis)
jpeg("Figure3G.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3))
dev.off()
pdf("Figure3G.pdf", h = 20, w =12)
print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = scale_fill_viridis(), annotation_colors = aka3))
dev.off()

jpeg("Figure3G_scaled.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"Blues")), annotation_colors = aka3, scale = "row"))
dev.off()

jpeg("Figure3G_no_log.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3))
dev.off()
jpeg("Figure3G_no_log_scaled.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3, scale = "row"))
dev.off()

#clip_scale
heatdata2 = t(scale(t(heatdata)))
max(heatdata2)
min(heatdata2)
heatdata2[heatdata2>2]=2
heatdata2[heatdata2<(-2)]=(-2)
jpeg("Figure3G_test.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata2, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3))
dev.off()

heatdata3 = apply(heatdata, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
jpeg("Figure3G_test2.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(heatdata3, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3))
dev.off()

Dicty3_small2 = ScaleData(Dicty3_small, assay = "SCT", features = rownames(Dicty3_small))
mat2 = Dicty3_small2@assays$SCT@data

mat3 = mat2[rownames(heatdata),colnames(heatdata)]
mat3 = t(scale(t(mat3)))
max(mat3)
min(mat3)
mat3[mat3>2]=2
mat3[mat3<(-2)]=(-2)

jpeg("Figure3G_scale_SCT.jpeg", h = 20, w =12, unit = "in", res = 300)
print(pheatmap(mat3, show_colnames = F, cluster_cols = F, annotation_col = ann, cluster_rows = T, color = rev(brewer.pal(9,"RdBu")), annotation_colors = aka3))
dev.off()

Dicty3_small = readRDS("Dicty3_small.rds")
write.table(Dicty3_small@assays$RNA@counts, file = "raw_sc_Dicty_3_stage_ds.tsv", sep = "\t", col.names = NA)
write.table(Dicty3_small@assays$SCT@data, file = "normalized_sc_Dicty_3_stage_ds.tsv", sep = "\t", col.names = NA)


#For GEO submission
df = Dicty3_small@meta.data
df2 = df[,1:7]
UMAP = Dicty3_small@reductions$umap@cell.embeddings
df2 = merge(df, UMAP, by = 0, all = T)  
harmony = Dicty3_small@reductions$harmony@cell.embeddings
rownames(df2) = df2[,1]
df2 = df2[,-1]
df3 = merge(df2, harmony, by = 0, all = T)  
rownames(df3) = df3[,1]
df3 = df3[,-1]
write.table(df3, file = "Metadata_Dicty_3_stage_ds.tsv", sep = "\t", col.names = NA)


