V <- read.delim("~/Dicty_chip-seq_correlation/vegetative_bin_counts", header = T)
V
scatter_plot_VA_VB_input <- ggplot(data = V, aes(V$X.VA_input., V$X.VB_input.))
cor.test(V$X.VA_input., V$X.VB_input., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_input + geom_point() + labs(x = "VA_input", y = "VB_input", subtitle = "r=0.96") + geom_smooth(method="lm") 

scatter_plot_VA_VB_K4me1 <- ggplot(data = V, aes(V$X.VA_K4me1., V$X.VB_K4me1.))
cor.test(V$X.VA_K4me1., V$X.VB_K4me1., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_K4me1 + geom_point() + labs(x = "VA_K4me1", y = "VB_K4me1", subtitle = "r=0.87") + geom_smooth(method="lm") 

scatter_plot_VA_VB_K4me3 <- ggplot(data = V, aes(V$X.VA_K4me3., V$X.VB_K4me3.))
cor.test(V$X.VA_K4me3., V$X.VB_K4me3., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_K4me3 + geom_point() + labs(x = "VA_K4me3", y = "VB_K4me3", subtitle = "r=0.98") + geom_smooth(method="lm") 

scatter_plot_VA_VB_K9me3 <- ggplot(data = V, aes(V$X.VA_K9me3., V$X.VB_K9me3.))
cor.test(V$X.VA_K9me3., V$X.VB_K9me3., method = "pearson", conf.leFel = 0.95)
scatter_plot_VA_VB_K9me3 + geom_point() + labs(x = "VA_K9me3", y = "VB_K9me3", subtitle = "r=0.97") + geom_smooth(method="lm") 


scatter_plot_VA_VB_K27ac <- ggplot(data = V, aes(V$X.VA_K27ac., V$X.VB_K27ac.))
cor.test(V$X.VA_K27ac., V$X.VB_K27ac., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_K27ac + geom_point() + labs(x = "VA_K27ac", y = "VB_K27ac", subtitle = "r=0.98") + geom_smooth(method="lm") 

scatter_plot_VA_VB_K27me3 <- ggplot(data = V, aes(V$X.VA_K27me3., V$X.VB_K27me3.))
cor.test(V$X.VA_K27me3., V$X.VB_K27me3., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_K27me3 + geom_point() + labs(x = "VA_K27me3", y = "VB_K27me3", subtitle = "r=0.99") + geom_smooth(method="lm")

scatter_plot_VA_VB_H3 <- ggplot(data = V, aes(V$X.VA_H3., V$X.VB_H3.))
cor.test(V$X.VA_H3., V$X.VB_H3., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_H3 + geom_point() + labs(x = "VA_H3", y = "VB_H3", subtitle = "r=0.98") + geom_smooth(method="lm")

scatter_plot_VA_VB_K36me3 <- ggplot(data = V, aes(V$X.VA_K36me3., V$X.VB_K36me3.))
cor.test(V$X.VA_K36me3., V$X.VB_K36me3., method = "pearson", conf.level = 0.95)
scatter_plot_VA_VB_K36me3 + geom_point() + labs(x = "VA_K36me3", y = "VB_K36me3", subtitle = "r=0.66") + geom_smooth(method="lm") 


F <- read.delim("~/Dicty_chip-seq_correlation/fruiting_bin_counts", header = T)
F
scatter_plot_FA_FB_input <- ggplot(data = table, aes(table$V4, table$V35))
cor.test(table$V4, table$V35, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_input + geom_point() + labs(x = "FA_input", y = "FB_input", subtitle = "r=0.99") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3 <- ggplot(data = table, aes(table$V5, table$V20))
cor.test(table$V5, table$V20, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3 + geom_point() + labs(x = "FA_H3", y = "FB_H3", subtitle = "r=0.99") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K4me1 <- ggplot(data = table, aes(table$V9, table$V22))
cor.test(table$V9, table$V22, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K4me1 + geom_point() + labs(x = "FA_H3K4me1", y = "FB_H3K4me1", subtitle = "r=0.80") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K4me3 <- ggplot(data = table, aes(table$V6, table$V21))
cor.test(table$V6, table$V21, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K4me3 + geom_point() + labs(x = "FA_H3K4me3", y = "FB_H3K4me3", subtitle = "r=0.91") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K9me3 <- ggplot(data = table, aes(table$V10, table$V25))
cor.test(table$V10, table$V25, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K9me3 + geom_point() + labs(x = "FA_H3K9me3", y = "FB_H3K9me3", subtitle = "r=0.97") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K27ac <- ggplot(data = table, aes(table$V7, table$V22))
cor.test(table$V7, table$V22, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K27ac + geom_point() + labs(x = "FA_H3K27ac", y = "FB_H3K27ac", subtitle = "r=0.97") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K27me3 <- ggplot(data = table, aes(table$V11, table$V26))
cor.test(table$V11, table$V26, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K27me3+ geom_point() + labs(x = "FA_H3K27me3", y = "FB_H3K27me3", subtitle = "r=0.91") + geom_smooth(method="lm") 

scatter_plot_FA_FB_H3K36me3 <- ggplot(data = table, aes(table$V8, table$V23))
cor.test(table$V8, table$V23, method = "pearson", conf.level = 0.95)
scatter_plot_FA_FB_H3K36me3+ geom_point() + labs(x = "FA_H3K36me3", y = "FB_H3K36me3", subtitle = "r=0.91") + geom_smooth(method="lm") 


M <- read.delim("~/Dicty_chip-seq_correlation/mound_bin_counts", header = T)
M
scatter_plot_MA_MB_input <- ggplot(data = M, aes(M$X.MA_input., M$X.MB_input.))
cor.test(M$X.MA_input., M$X.MB_input., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_input + geom_point() + labs(x = "MA_input", y = "MB_input", subtitle = "r=0.98") + geom_smooth(method="lm") 

scatter_plot_MA_MB_K4me1 <- ggplot(data = M, aes(M$X.MA_K4me1., M$X.MB_K4me1.))
cor.test(M$X.MA_K4me1., M$X.MB_K4me1., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K4me1 + geom_point() + labs(x = "MA_K4me1", y = "MB_K4me1", subtitle = "r=0.73") + geom_smooth(method="lm") 

scatter_plot_MA_MB_K4me3 <- ggplot(data = M, aes(M$X.MA_K4me3., M$X.MB_K4me3.))
cor.test(M$X.MA_K4me3., M$X.MB_K4me3., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K4me3 + geom_point() + labs(x = "MA_K4me3", y = "MB_K4me3", subtitle = "r=0.95") + geom_smooth(method="lm") 

scatter_plot_MA_MB_K9me3 <- ggplot(data = M, aes(M$X.MA_K9me3., M$X.MB_K9me3.))
cor.test(M$X.MA_K9me3., M$X.MB_K9me3., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K9me3 + geom_point() + labs(x = "MA_K9me3", y = "MB_K9me3", subtitle = "r=0.95") + geom_smooth(method="lm") 

scatter_plot_MA_MB_K27ac <- ggplot(data = M, aes(M$X.MA_K27ac., M$X.MB_K27ac.))
cor.test(M$X.MA_K27ac., M$X.MB_K27ac., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K27ac + geom_point() + labs(x = "MA_K27ac", y = "MB_K27ac", subtitle = "r=0.98") + geom_smooth(method="lm") 

scatter_plot_MA_MB_K27me3 <- ggplot(data = M, aes(M$X.MA_K27me3., M$X.MB_K27me3.))
cor.test(M$X.MA_K27me3., M$X.MB_K27me3., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K27me3 + geom_point() + labs(x = "MA_K27me3", y = "MB_K27me3", subtitle = "r=0.98") + geom_smooth(method="lm")

scatter_plot_MA_MB_H3 <- ggplot(data = M, aes(M$X.MA_H3., M$X.MB_H3.))
cor.test(M$X.MA_H3., M$X.MB_H3., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_H3 + geom_point() + labs(x = "MA_H3", y = "MB_H3", subtitle = "r=0.93") + geom_smooth(method="lm")

scatter_plot_MA_MB_K36me3 <- ggplot(data = M, aes(M$X.MA_K36me3., M$X.MB_K36me3.))
cor.test(M$X.MA_K36me3., M$X.MB_K36me3., method = "pearson", conM.leMel = 0.95)
scatter_plot_MA_MB_K36me3 + geom_point() + labs(x = "MA_K36me3", y = "MB_K36me3", subtitle = "r=0.87") + geom_smooth(method="lm")

S <- read.delim("~/Dicty_chip-seq_correlation/streaming_bin_counts", header = T)
S
scatter_plot_SA_SB_input <- ggplot(data = table, aes(table$V12, table$V27))
cor.test(table$V12, table$V27, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_input + geom_point() + labs(x = "SA_input", y = "SB_input", subtitle = "r=0.90") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3 <- ggplot(data = table, aes(table$V13, table$V28))
cor.test(table$V13, table$V28, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3 + geom_point() + labs(x = "SA_H3", y = "SB_H3", subtitle = "r=0.74") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K4me1 <- ggplot(data = table, aes(table$V17, table$V32))
cor.test(table$V17, table$V32, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K4me1 + geom_point() + labs(x = "SA_H3K4me1", y = "SB_H3K4me1", subtitle = "r=0.88") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K4me3 <- ggplot(data = table, aes(table$V14, table$V29))
cor.test(table$V14, table$V29, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K4me3 + geom_point() + labs(x = "SA_H3K4me3", y = "SB_H3K4me3", subtitle = "r=0.98") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K9me3 <- ggplot(data = table, aes(table$V18, table$V33))
cor.test(table$V18, table$V33, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K9me3 + geom_point() + labs(x = "SA_H3K9me3", y = "SB_H3K9me3", subtitle = "r=0.92") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K27ac <- ggplot(data = table, aes(table$V15, table$V30))
cor.test(table$V15, table$V30, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K27ac + geom_point() + labs(x = "SA_H3K27ac", y = "SB_H3K27ac", subtitle = "r=0.71") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K27me3 <- ggplot(data = table, aes(table$V19, table$V34))
cor.test(table$V19, table$V34, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K27me3+ geom_point() + labs(x = "SA_H3K27me3", y = "SB_H3K27me3", subtitle = "r=0.99") + geom_smooth(method="lm") 

scatter_plot_SA_SB_H3K36me3 <- ggplot(data = table, aes(table$V16, table$V31))
cor.test(table$V16, table$V31, method = "pearson", conf.level = 0.95)
scatter_plot_SA_SB_H3K36me3+ geom_point() + labs(x = "SA_H3K36me3", y = "SB_H3K36me3", subtitle = "r=0.92") + geom_smooth(method="lm") 

