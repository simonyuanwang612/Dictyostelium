library(clusterProfiler)
library(ggplot2)
require(DOSE)
require(clusterProfiler)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(rJava)
library("RDAVIDWebService")
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