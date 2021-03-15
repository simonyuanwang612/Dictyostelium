counts_strand_FV <- read.delim("~/Desktop/Dicty_RNA_seq_kallisto/kallisto_dicty_matrix.txt", comment.char="#",  row.names=1, stringsAsFactors=FALSE)
x_FV <-  counts_strand_FV[c( "lib13", "lib14", "lib19", "lib20")]
head(x_FV)
tail(x_FV)
dim(x_FV)

#create groups
group_FV <- c(rep("V", 2), rep("F", 2))
group_FV

#Building the edgeR Object
#create DGEList object to work on, where are the columns are the counts
#Put the counts and other information into a DGEList object
y_FV <- DGEList(counts=x_FV[,1:4], group=group_FV)
y_FV
dim(y_FV)

# Human introduced error
#Filter out low count reads since it would be impossible to detect differential expression. The method used in the
#edgeR vignette is to keep only those genes that have at least 1 read per million in at least 2 samples. Therefore, Filter out very lowly expressed tags, 
#keeping genes that are expressed at a reasonable level in at least one condition. Since the group sizes are 2, keep genes that achieve at least 
#one count per million (cpm) in 2 samples
#Once this filtering is done, calculate the normalization factors which correct for the different compositions of the samples. The effective
#library sizes are then the product of the actual library sizes and these factors.

keep_FV <- rowSums(cpm(y_FV)>1) >= 2
y_FV <- y_FV[keep_FV,]
dim(y_FV)
#Re-compute the library sizes
y_FV$samples$lib.size <- colSums(y_FV$counts)
y_FV$samples

##TMM Normalization
#Compute effective library sizes using Trimmed Mean of M Values (TMM) normalization
#This accounts for differences in transcriptome sizes
y_FV <- calcNormFactors(y_FV)
y_FV$samples

#effective library sizes
# This is relative normalization
y_FV$samples$lib.size*y_FV$samples$norm.factors

##Data exploration
#Generate a Multi-dimensional Scaling (MDS) plot
#An MDS plot measures the similarity of the samples and projects this measure into 2-dimensions.
#An MDS plot show distances, in terms of biological coefficient of variation (BCV), between samples
#Dimension 1 separates the Normal from TNBC samples. The MDS plot shows the replicates are consistent. 
#This is a good indication edgeR analysis will produce adequate DE genes.
#To view MDS plot
plotMDS(y_FV, main ="MDS Plot for TMM Normalized Count Data F_V")

#Write plot as a pdf
pdf("MDS_Plot_F_V.pdf" , width = 11 , height = 7) # in inches
plotMDS(y_FV, main ="MDS Plot for TMM Normalized Count Data_FV")
dev.off() # this tells [R] to close and stop writing to the pdf


##Estimating common and tagwise dispersion
#The common dispersion estimates the overall BCV of the dataset, averaged over all genes with the qCML method
y_FV <- estimateCommonDisp(y_FV, verbose=TRUE)
# Disp = 0.0138 , BCV = 0.117 

#Now estimate gene-specific (tagwise) dispersions with the qCML method
y_FV <- estimateTagwiseDisp(y_FV)
summary(y_FV$tagwise.dispersion)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0002  0.0002  0.0002  0.0049  0.0071  0.2938

#Plot the estimated dispersions 
#In general, BCV is higher with lower average log CPM
#To view estimated dispersions plot
plotBCV(y_FV, main = "Plot of Estimated Dispersions_F_V")


#Write plot as a pdf
pdf("Est_Dispersions_Plot_F_V.pdf" , width = 11 , height = 7) # in inches
plotBCV(y_FV, main = "Plot of Estimated Dispersions_F_M")
dev.off() # this tells [R] to close and stop writing to the pdf

##Differential expression
#FDR is less significant than the original P value, correcting for the P value. Set false discovery treshold.
#Compute NB exact genewise tests for differential expression between samples
et_FV <- exactTest(y_FV, pair=2:1)
options(digits = 3) #print only 3 digits
top_FV <- topTags(et_FV)
top_FV

#Check the individual cpm values for the top genes
#Check up/down expression to  *VERY IMPORTANT*
cpm(y_FV)[rownames(top_FV), ]

#print the total number of DE genes at 5% FDR
summary(de_FV <- decideTestsDGE(et_FV))

#Plot the log-fold-changes with a smear plot, highlighting the DE genes
#Blue lines indicate 2-fold changes
#to view Smear Plot
detags_FV <- rownames(y_FV)[as.logical(de_FV)]
plotSmear(et_FV, de.tags=detags_FV, main = "Smear Plot F-V")

abline(h=c(-1, 1), col="blue")

#Write plot as a pdf
pdf("Smear_Plot_F_V.pdf" , width = 7 , height = 7) # in inches
plotSmear(et_FV, de.tags=detags_FV, main = "Smear Plot F-V")
abline(h=c(-1, 1), col="blue")
dev.off() # this tells [R] to close and stop writing to the pdf

?plotSmear
#Create an object called "resultsTbl.tagwise" containing the four output parameters: logFC, logCPM, PValue,FDR
resultsTbl.tagwise_FV <- topTags(et_FV, n = nrow(et_FV$table))$table


#Check resultsTbl.tagwise
head(resultsTbl.tagwise_FV)
dim(resultsTbl.tagwise_FV)

#Generate objects according to three significance levels for tagwise dispersion estimates: 0.05, 0.001, 1e-06
#These will be used for GOseq analyses in classes 5 and 6.
#Names/IDs of DE genes
de.genes.tagwise_0.05_FV <- rownames(resultsTbl.tagwise_FV)[resultsTbl.tagwise_FV$FDR <= 0.05]
de.genes.tagwise_0.001_FV <- rownames(resultsTbl.tagwise_FV)[resultsTbl.tagwise_FV$FDR <= 0.001]
de.genes.tagwise_1e_06_FV <- rownames(resultsTbl.tagwise_FV)[resultsTbl.tagwise_FV$FDR <= 1e-06]

#Check number of genes at each significance level
length(de.genes.tagwise_0.05_FV) #2174 DE genes
length(de.genes.tagwise_0.001_FV) #1245 DE genes
length(de.genes.tagwise_1e_06_FV) #736 DE genes

#Check percentage of total genes
length(de.genes.tagwise_0.05_FV)/nrow(resultsTbl.tagwise_FV)*100 #22
length(de.genes.tagwise_0.001_FV)/nrow(resultsTbl.tagwise_FV)*100 #12.6
length(de.genes.tagwise_1e_06_FV)/nrow(resultsTbl.tagwise_FV)*100 #7.44

#Check Up/Down regulated summaries for tagwise results
summary(de_FV <- decideTestsDGE(et_FV, p.value=0.05))
#Down    854
#NotSig 7927
#Up      1320
summary(de_FV <- decideTestsDGE(et_FV, p.value=0.001))
#Down    440
#NotSig 8856
#Up      805
summary(de_FV <- decideTestsDGE(et_FV, p.value=1e-06))
#Down    201
#NotSig 9365
#Up      535
head(de.genes.tagwise_1e_06_FV)

#Double check
is.vector(de.genes.tagwise_1e_06_FV)

##Format Results and count data form DGEList object(y)
colnames(resultsTbl.tagwise_FV) <- c("logFC", "logCPM","PValue", "FDR")
wh.rows.tagwise_FV <- match(rownames(resultsTbl.tagwise_FV), rownames(y_FV$counts))


#Check for all 10101 genes
length(wh.rows.tagwise_FV)

#Combine exact test results, count data from DGEList object(y) and tagwise dispersion
combResults.tagwise_FV <- cbind(resultsTbl.tagwise_FV,
                                "Tgw.Disp" = y_FV$tagwise.dispersion[wh.rows.tagwise_FV],
                                "UpDown.Tgw" = decideTestsDGE(et_FV, p.value = 0.05)[wh.rows.tagwise_FV],
                                y_FV$counts[wh.rows.tagwise_FV,])

#Check combResults.tagwise 
head(combResults.tagwise_FV)
dim(combResults.tagwise_FV)

#Create significance level objects from original counts data (x) and FDR vectors
de.tag_05_FV<-x_FV[which(row.names(x_FV)%in%(de.genes.tagwise_0.05_FV)),]
de.tag_001_FV<-x_FV[which(row.names(x_FV)%in%(de.genes.tagwise_0.001_FV)),]
de.tag_1e_06_FV<-x_FV[which(row.names(x_FV)%in%(de.genes.tagwise_1e_06_FV)),]

#Check
head(de.tag_05_FV)
dim(de.tag_05_FV)
head(de.tag_001_FV)
dim(de.tag_001_FV)
head(de.tag_1e_06_FV)
dim(de.tag_1e_06_FV)

#Check again
head(combResults.tagwise_FV)
head(resultsTbl.tagwise_FV)

dim(combResults.tagwise_FV)

##Merge all output into FDR specific objects
combResults.tagwise_05_FV<-merge(combResults.tagwise_FV [,1:6], de.tag_05_FV, by="row.names")
head(combResults.tagwise_05_FV)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_05_FV)[1] ="target_id"
head(combResults.tagwise_05_FV)
dim(combResults.tagwise_05_FV)

#Create edgeR_05 output object and reorder columns 
edgeR_05_FV=(combResults.tagwise_05_FV)
head(edgeR_05_FV)
dim(edgeR_05_FV)

#Repeat for FDR of 0.001
combResults.tagwise_001_FV<-merge(combResults.tagwise_FV [,1:6], de.tag_001_FV, by="row.names")
head(combResults.tagwise_001_FV)

#Rename column 1 as "Gene"
colnames(combResults.tagwise_001_FV)[1] ="target_id"
head(combResults.tagwise_001_FV)
dim(combResults.tagwise_001_FV)

#Create edgeR_001 output object and reorder columns
edgeR_001_FV=combResults.tagwise_001_FV
head(edgeR_001_FV, 20)
dim(edgeR_001_FV)

#Repeat for FDR of 1e06
combResults.tagwise_1e_06_FV<-merge(combResults.tagwise_FV [,1:6], de.tag_1e_06_FV, by="row.names")
head(combResults.tagwise_1e_06_FV)

#Rename column 1 as "target_id"
colnames(combResults.tagwise_1e_06_FV)[1] ="target_id"
head(combResults.tagwise_1e_06_FV)
dim(combResults.tagwise_1e_06_FV)


#Create edgeR_1e_06 output object and reorder columns
edgeR_1e_06_FV = combResults.tagwise_1e_06_FV
head(edgeR_1e_06_FV)
dim(edgeR_1e_06_FV)

#Sort results by FDR as a check for row names.
edgeR_05_out_FV<-edgeR_05_FV[order(edgeR_05_FV[,7],decreasing=FALSE),]
head(edgeR_05_out_FV)

edgeR_001_out_FV<-edgeR_001_FV[order(edgeR_001_FV[,7],decreasing=FALSE),]
head(edgeR_001_out_FV)

edgeR_1e_06_out_FV<-edgeR_1e_06_FV[order(edgeR_1e_06_FV[,7],decreasing=FALSE),]
head(edgeR_1e_06_out_FV)



edgeR_05_out_geneid_FV <- merge(edgeR_05_out_FV, mart_export, by = "target_id")
edgeR_001_out_geneid_FV <- merge(edgeR_001_out_FV, mart_export, by = "target_id")
edgeR_1e_06_out_geneid_FV <- merge(edgeR_1e_06_out_FV, mart_export, by = "target_id")

write.table(edgeR_05_out_geneid_FV, file="edgeR_tagwise_05_output_F_V_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_001_out_geneid_FV, file="edgeR_tagwise_001_output_F_V_kali.txt", row.names=FALSE, sep="\t")
write.table(edgeR_1e_06_out_geneid_FV, file="edgeR_tagwise_1e_06_output_F_V_kali.txt", row.names=FALSE, sep="\t")

dim(edgeR_05_out_geneid_FV)
dim(edgeR_001_out_geneid_FV)
dim(edgeR_1e_06_out_geneid_FV)