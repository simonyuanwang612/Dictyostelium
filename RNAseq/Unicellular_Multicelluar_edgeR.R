library(edgeR)
counts_strand_UM <- read.delim("~/Desktop/Dicty_RNA_seq_kallisto/kallisto_dicty_matrix.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
x_UM <-  counts_strand_UM[c( "lib13", "lib14", "lib17", "lib18", "lib19", "lib20")]
head(x_UM)
dim(x_UM)
group_UM <- c(rep("Uni", 2), rep("Multi", 6))
group_UM
y_UM <- DGEList(counts=x_UM[,1:8], group=group_UM)
dim(y_UM)
y_UM
keep_UM <- rowSums(cpm(y_UM)>1) >= 2
y_UM <- y_UM[keep_UM,]
dim(y_UM)
y_UM$samples$lib.size <- colSums(y_UM$counts)
y_UM$samples
y_UM <- calcNormFactors(y_UM)
y_UM$samples
y_UM$samples$lib.size*y_UM$samples$norm.factors
plotMDS(y_UM, main ="MDS Plot for TMM Normalized Count Data uni_multi")
?plotMDS
pdf("MDS_Plot_uni_multi.pdf" , width = 11 , height = 7)
plotMDS(y_UM, main ="MDS Plot for TMM Normalized Count Data uni-multi")
dev.off() 
y_UM <- estimateCommonDisp(y_UM, verbose=TRUE)
y_UM <- estimateTagwiseDisp(y_UM)
summary(y_UM$tagwise.dispersion)
plotBCV(y_UM, main = "Plot of Estimated Dispersions_uni_multi")
pdf("Est_Dispersions_Plot_uni_multi.pdf" , width = 11 , height = 7)
plotBCV(y_UM, main = "Plot of Estimated Dispersions_uni_multi")
dev.off()
et_UM <- exactTest(y_UM, pair=2:1)
options(digits = 3)
top_UM <- topTags(et_UM)
top_UM
cpm(y_UM)[rownames(top_UM), ]
summary(de_UM <- decideTestsDGE(et_UM))
detags_UM <- rownames(y_UM)[as.logical(de_UM)]
plotSmear(et_UM, de.tags=detags_UM, main = "Smear Plot Unicellular VS Multicellular")
abline(h=c(-1, 1), col="blue")

pdf("Smear_Plot_uni_multi.pdf" , width = 7 , height = 7) 
plotSmear(et_UM, de.tags=detags_UM, main = "Smear Plot Uni-Multi")
abline(h=c(-1, 1), col="blue")
dev.off()
resultsTbl.tagwise_UM <- topTags(et_UM, n = nrow(et_UM$table))$table
head(resultsTbl.tagwise_UM)
dim(resultsTbl.tagwise_UM)
de.genes.tagwise_0.05_UM <- rownames(resultsTbl.tagwise_UM)[resultsTbl.tagwise_UM$FDR <= 0.05]
length(de.genes.tagwise_0.05_UM)
length(de.genes.tagwise_0.05_UM)/nrow(resultsTbl.tagwise_UM)*100
summary(de_UM <- decideTestsDGE(et_UM, p.value=0.05))
head(de.genes.tagwise_0.05_UM)
is.vector(de.genes.tagwise_0.05_UM)
wh.rows.tagwise_UM <- match(rownames(resultsTbl.tagwise_UM), rownames(y_UM$counts))
colnames(resultsTbl.tagwise_UM) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise_UM <- match(rownames(resultsTbl.tagwise_UM), rownames(y_UM$counts))
length(wh.rows.tagwise_UM)
combResults.tagwise_UM <- cbind(resultsTbl.tagwise_UM,
                               "Tgw.Disp" = y_UM$tagwise.dispersion[wh.rows.tagwise_UM],
                               "UpDown.Tgw" = decideTestsDGE(et_UM, p.value = 0.05)[wh.rows.tagwise_UM],
                               y_UM$counts[wh.rows.tagwise_UM,])

head(combResults.tagwise_UM)
dim(combResults.tagwise_UM)

de.tag_05_UM<-x_UM[which(row.names(x_UM)%in%(de.genes.tagwise_0.05_UM)),]
head(de.tag_05_UM)
dim(de.tag_05_UM)
head(combResults.tagwise_UM)
head(resultsTbl.tagwise_UM)
dim(combResults.tagwise_UM)
combResults.tagwise_05_UM<-merge(combResults.tagwise_UM [,1:6], de.tag_05_UM, by="row.names")
head(combResults.tagwise_05_UM)
colnames(combResults.tagwise_05_UM)[1] ="target_id"
head(combResults.tagwise_05_UM)
dim(combResults.tagwise_05_UM)
edgeR_05_UM=(combResults.tagwise_05_UM)
head(edgeR_05_UM)
dim(edgeR_05_UM)
edgeR_05_out_UM <-edgeR_05_UM[order(edgeR_05_UM[,7],decreasing=FALSE),]
head(edgeR_05_out_UM)
edgeR_05_out_geneid_UM <- merge(edgeR_05_out_UM, mart_export, by = "target_id")
head(edgeR_05_out_geneid_UM)
write.table(edgeR_05_out_geneid_UM, file="edgeR_tagwise_05_output_Uni_Multi_kali.txt", row.names=FALSE, sep="\t")
