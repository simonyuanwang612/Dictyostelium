library(clusterProfiler)
library(ggplot2)

#Uni Biological Process

david_uni <- enrichDAVID(gene = uni_enrich_gene$gene_id,
                           idType = "DICTYBASE_ID",
                           annotation = "GOTERM_BP_ALL",
                           species = "DICTYOSTELIUM_DISCOIDEUM",
                           david.user = "yuan.wang@childrens.harvard.edu")
#Uni KEGG
david_uni_KEGG <- enrichDAVID(gene = uni_enrich_gene$gene_id,
                         idType = "DICTYBASE_ID",
                         annotation = "KEGG_PATHWAY",
                         species = "DICTYOSTELIUM_DISCOIDEUM",
                         david.user = "yuan.wang@childrens.harvard.edu")
# Uni bar plot biological process
u <- barplot(david_uni, showCategory=50)
ggplot2::ggsave(u, filename = "barplot_Uni.png", width=20, height=20)
# Uni dot plot biological process
u <- dotplot(david_uni, showCategory=50)
ggplot2::ggsave(u, filename = "dotplot_Uni.png", width=20, height=20)
# Uni bar plot KEGG
u_kegg <- barplot(david_uni_KEGG, showCategory=50)
ggplot2::ggsave(u_kegg, filename = "barplot_Uni_KEGG.png", width=20, height=20)




#Multi biological process
david_multi <- enrichDAVID(gene = multi_enrich_genes$gene_id,
                     idType = "DICTYBASE_ID",
                     annotation = "GOTERM_BP_ALL",
                     species = "DICTYOSTELIUM_DISCOIDEUM",
                     david.user = "yuan.wang@childrens.harvard.edu")
#Multi barplot biological process
m <- barplot(david_multi, showCategory=50)
ggplot2::ggsave(m, filename = "barplot_Multi.png", width=20, height=20)
#Multi dotplot biological process
m <- dotplot(david_multi, showCategory=50)
ggplot2::ggsave(m, filename = "dotplot_Multi.png", width=20, height=20)


#Uni Biological Process
david_uni_CC <- enrichDAVID(gene = uni_enrich_gene$gene_id,
                         idType = "DICTYBASE_ID",
                         annotation = "GOTERM_CC_ALL",
                         species = "DICTYOSTELIUM_DISCOIDEUM",
                         david.user = "yuan.wang@childrens.harvard.edu")
#uni CC emapplot
ema <- emapplot(david_uni_CC, showCategory=100)
ggplot2::ggsave(ema, filename = "emapplot_uni_CC.png", width=20, height=20)
#uni CC cnetplot
cne <- cnetplot(david_uni_CC, showCategory=100)
ggplot2::ggsave(cne, filename = "cnetplot_uni_CC.png", width=20, height=20)

#Multi biological process
david_multi_CC <- enrichDAVID(gene = multi_enrich_genes$gene_id,
                           idType = "DICTYBASE_ID",
                           annotation = "GOTERM_CC_ALL",
                           species = "DICTYOSTELIUM_DISCOIDEUM",
                           david.user = "yuan.wang@childrens.harvard.edu")
#multi CC emapplot
multi_ema <- emapplot(david_multi_CC, showCategory=100)
ggplot2::ggsave(multi_ema, filename = "emapplot_multi_CC.png", width=20, height=20)
#uni CC cnetplot
multi_cne <- cnetplot(david_multi_CC , showCategory=100)
ggplot2::ggsave(multi_cne, filename = "cnetplot_multi_CC.png", width=20, height=20)




search_kegg_organism('ddi', by='kegg_code')
dicty <- search_kegg_organism('Dictyostelium discoideum', by='scientific_name')
dim(dicty)
head(dicty)
kk <- enrichKEGG( gene = multi_uniprot$Uniprot,
                  organism = 'hsa',
                  pvalueCutoff = 0.05)
library(DOSE)
library(GO.db)
library(GSEABase)
library(org.Hs.eg.db)

mkk <- enrichKEGG(edgeR_tagwise_05_output_Uni_Multi_CDS$Gene,
                   organism = 'ddi')
mkk
barplot(mkk)

F_V <- enrichKEGG(F_V_test$Gene, organism = 'ddi')
F_V
dotplot(F_V)

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("RDAVIDWebService", version = "devel")

library(clusterProfiler)
data(geneList, package='DOSE')                 
de <- names(geneList)[1:100]
yy <- enrichKEGG(de, organism='hsa', pvalueCutoff=0.01)
barplot(yy)
sessionInfo()
source("https://bioconductor.org/biocLite.R")
biocLite("goseq")

# S_V
david_S_V <- enrichDAVID(gene = S_V$gene_id,
                         idType = "DICTYBASE_ID",
                         annotation = "GOTERM_BP_ALL",
                         species = "DICTYOSTELIUM_DISCOIDEUM",
                         david.user = "yuan.wang@childrens.harvard.edu")

s_v <- barplot(david_S_V, showCategory=20)
ggplot2::ggsave(s_v, filename = "barplot_S_V.png", width=20, height=20)

s_v_dotplot <- dotplot(david_S_V, showCategory=5)
ggplot2::ggsave(s_v_dotplot, filename = "dotplot_s_v.pdf", width=20, height=20)

# M_S
david_M_S <- enrichDAVID(gene = M_S$gene_id,
                        idType = "DICTYBASE_ID",
                        annotation = "GOTERM_BP_ALL",
                        species = "DICTYOSTELIUM_DISCOIDEUM",
                        david.user = "yuan.wang@childrens.harvard.edu")

m_s <- barplot(david_M_S, showCategory=20)
ggplot2::ggsave(m_s, filename = "barplot_M_S.png", width=20, height=20)

# F_M

david_F_M <- enrichDAVID(gene = F_M$gene_id,
                         idType = "DICTYBASE_ID",
                         annotation = "GOTERM_BP_ALL",
                         species = "DICTYOSTELIUM_DISCOIDEUM",
                         david.user = "yuan.wang@childrens.harvard.edu")

f_m <- barplot(david_F_M, showCategory=20)
ggplot2::ggsave(f_m, filename = "barplot_F_M.png", width=20, height=20)

# V_F_Vhigh

david_V_F_vhigh <- enrichDAVID(gene = V_F_vhigh$gene_id,
                         idType = "DICTYBASE_ID",
                         annotation = "GOTERM_BP_ALL",
                         species = "DICTYOSTELIUM_DISCOIDEUM",
                         david.user = "yuan.wang@childrens.harvard.edu")

v_f_vhigh <- barplot(david_V_F_vhigh , showCategory=20)
ggplot2::ggsave(v_f_vhigh, filename = "barplot_V_F_vhigh.png", width=20, height=20)

# F_V_Fhigh

david_F_V_Fhigh <- enrichDAVID(gene = F_V_fhigh$gene_id,
                               idType = "DICTYBASE_ID",
                               annotation = "GOTERM_BP_ALL",
                               species = "DICTYOSTELIUM_DISCOIDEUM",
                               david.user = "yuan.wang@childrens.harvard.edu")

v_f_fhigh<- barplot(david_F_V_Fhigh, showCategory=20)
ggplot2::ggsave(v_f_fhigh, filename = "barplot_F_V_fhigh.png", width=20, height=20)

# Grouping different stage and do the GO analysis
Grouping <- read_excel("Grouping.xlsx")
lapply(Grouping, head)
david_combine <- compareCluster(geneCluster = Grouping, fun = "enrichDAVID", idType = "DICTYBASE_ID",
                                annotation = "GOTERM_BP_ALL", species = "DICTYOSTELIUM_DISCOIDEUM",
                                david.user = "yuan.wang@childrens.harvard.edu")
head(as.data.frame(david_combine))
v_f_david_combine <- dotplot(david_combine , showCategory=20)
ggplot2::ggsave(plot = david_combine, filename = "david_combine.png", width=20, height=20)
v_f_david_combine
