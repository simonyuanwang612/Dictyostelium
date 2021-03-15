merge_table <- merge(Gene_exp_PA, Gene_exp_SW, by = "gene_id")
head(merge_table)
dim(merge_table)
write.table(merge_table, "/Users/simon/Downloads/merge_table.txt", sep="\t")


# Cross check between SW and PA
library(ggplot2)
library(dplyr)
scatter_plot_V <- ggplot(data=merge_table, aes(merge_table$Vegetative_SW, merge_table$Vegetative_PA))
cor.test(log10(merge_table$Vegetative_SW+1), log10(merge_table$Vegetative_PA+1), method="pearson", conf.level = 0.05)
scatter_plot_V + geom_point() + labs(x= "Vegetative_SW", y ="Vegetative_PA", subtitle = "r = 0.72") + geom_smooth(method = "lm") + scale_x_log10() + scale_y_log10()

scatter_plot_S <- ggplot(data=merge_table, aes(merge_table$Streaming_SW, merge_table$Streaming_PA))
cor.test(log10(merge_table$Streaming_SW+1), log10(merge_table$Streaming_PA+1), method="pearson", conf.level = 0.05)
scatter_plot_S + geom_point() + labs(x= "Streaming_SW", y ="Streaming_PA", subtitle = "r = 0.74") + geom_smooth(method = "lm") +scale_x_log10() + scale_y_log10()

scatter_plot_M <- ggplot(data=merge_table, aes(merge_table$Mound_SW, merge_table$Mound_PA))
cor.test(log10(merge_table$Mound_SW+1), log10(merge_table$Mound_PA+1), method="pearson", conf.level = 0.05)
scatter_plot_M + geom_point() + labs(x= "Mound_SW", y ="Mound_PA", subtitle = "r = 0.81") + geom_smooth(method = "lm") + scale_x_log10() + scale_y_log10()

scatter_plot_F <- ggplot(data=merge_table, aes(merge_table$Fruiting_SW, merge_table$Fruiting_PA))
cor.test(log10(merge_table$Fruiting_SW+1), log10(merge_table$Fruiting_PA+1), method="pearson", conf.level = 0.05)
scatter_plot_F + geom_point() + labs(x= "Fruiting_SW", y ="Fruiting_PA", subtitle = "r = 0.83") + geom_smooth(method = "lm") + scale_x_log10() + scale_y_log10()

# Biological correlation in PA's result
scatter_PA_V1_V2 <- ggplot(data=PA_general, aes(PA_general$bio1.hr00, PA_general$bio2.hr00))
cor.test(PA_general$bio1.hr00, PA_general$bio2.hr00, method="pearson", conf.level = 0.05)
scatter_PA_V1_V2 + geom_point() + labs(x= "Vegetative_PA_rep1", y ="Vegetative_PA_rep2", subtitle = "r = 0.8") + geom_smooth(method = "lm")

scatter_PA_S1_S2 <- ggplot(data=PA_general, aes(PA_general$bio1.hr08, PA_general$bio2.hr08))
cor.test(PA_general$bio1.hr08, PA_general$bio2.hr08, method="pearson", conf.level = 0.05)
scatter_PA_S1_S2 + geom_point() + labs(x= "Streaming_PA_rep1", y ="Streaming_PA_rep2", subtitle = "r = 0.72") + geom_smooth(method = "lm")

scatter_PA_M1_M2 <- ggplot(data=PA_general, aes(PA_general$bio1.hr12, PA_general$bio2.hr12))
cor.test(PA_general$bio1.hr12, PA_general$bio2.hr12, method="pearson", conf.level = 0.05)
scatter_PA_M1_M2 + geom_point() + labs(x= "Mound_PA_rep1", y ="Mound_PA_rep2", subtitle = "r = 0.82") + geom_smooth(method = "lm")

scatter_PA_F1_F2 <- ggplot(data=PA_general, aes(PA_general$bio1.hr24, PA_general$bio2.hr24))
cor.test(PA_general$bio1.hr24, PA_general$bio2.hr24, method="pearson", conf.level = 0.05)
scatter_PA_F1_F2 + geom_point() + labs(x= "SFruiting_PA_rep1", y ="Fruiting_PA_rep2", subtitle = "r = 0.51") + geom_smooth(method = "lm")

# Biological correlation in SW's result
scatter_SW_V1_V2 <- ggplot(data=GSE137602_RNA_seq_normalized_count_matrix, aes(GSE137602_RNA_seq_normalized_count_matrix$lib13, GSE137602_RNA_seq_normalized_count_matrix$lib14))
cor.test(GSE137602_RNA_seq_normalized_count_matrix$lib13, GSE137602_RNA_seq_normalized_count_matrix$lib14, method="pearson", conf.level = 0.05)
scatter_SW_V1_V2 + geom_point() + labs(x= "Vegetative_SW_rep1", y ="Vegetative_SW_rep2", subtitle = "r = 0.98") + geom_smooth(method = "lm")

scatter_SW_S1_S2 <- ggplot(data=GSE137602_RNA_seq_normalized_count_matrix, aes(GSE137602_RNA_seq_normalized_count_matrix$lib15, GSE137602_RNA_seq_normalized_count_matrix$lib16))
cor.test(GSE137602_RNA_seq_normalized_count_matrix$lib15, GSE137602_RNA_seq_normalized_count_matrix$lib16, method="pearson", conf.level = 0.05)
scatter_SW_S1_S2 + geom_point() + labs(x= "Streaming_SW_rep1", y ="Streaming_SW_rep2", subtitle = "r = 0.95") + geom_smooth(method = "lm")

scatter_SW_M1_M2 <- ggplot(data=GSE137602_RNA_seq_normalized_count_matrix, aes(GSE137602_RNA_seq_normalized_count_matrix$lib17, GSE137602_RNA_seq_normalized_count_matrix$lib18))
cor.test(GSE137602_RNA_seq_normalized_count_matrix$lib17, GSE137602_RNA_seq_normalized_count_matrix$lib18, method="pearson", conf.level = 0.05)
scatter_SW_M1_M2 + geom_point() + labs(x= "Mound_SW_rep1", y ="Mound_SW_rep2", subtitle = "r = 0.98") + geom_smooth(method = "lm")

scatter_SW_F1_F2 <- ggplot(data=GSE137602_RNA_seq_normalized_count_matrix, aes(GSE137602_RNA_seq_normalized_count_matrix$lib19, GSE137602_RNA_seq_normalized_count_matrix$lib20))
cor.test(GSE137602_RNA_seq_normalized_count_matrix$lib19, GSE137602_RNA_seq_normalized_count_matrix$lib20, method="pearson", conf.level = 0.05)
scatter_SW_F1_F2 + geom_point() + labs(x= "SFruiting_SW_rep1", y ="Fruiting_SW_rep2", subtitle = "r = 0.98") + geom_smooth(method = "lm")


