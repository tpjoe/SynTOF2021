library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyr)




#### Frequency plots ----------------------------------------------------------
df <- list()
# Importing data
for (region in c('BA9', 'DLCau', 'Hipp')){
    df[[region]] <- read.csv(paste0('R_py_exchange/df_freq_', region, '_noStd.csv'))[, -c(1)]
    group <- df[[region]][, 1]
    samples <- df[[region]][, 2]
    df[[region]] <- df[[region]][, -c(1, 2)]
    # colnames(df)[3:ncol(df)] <- sapply(colnames(df)[3:ncol(df)], function(x) strsplit(x, 'X')[[1]][2])
}

df <- do.call(cbind.data.frame, df)
df <- cbind.data.frame(group=group, samples=samples, df)

freq <- df[df$group!='ODC', ]





#### Importing data ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
    # df_ <- read.csv(paste0('R_py_exchange/df_medianAllMarkers_', region, '_noStd_exp0.csv'))[, -1]
    df_ <- read.csv(paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp2mc_5_13.csv'))[, -1]
    if (region != 'BA9') {
        #check that samples are in the same order
        print(all(df[, 2] == df_[, 2]))
        colnames(df_) <- paste0(region, '_', colnames(df_))
        df <- cbind.data.frame(df, df_[, 3:ncol(df_)])
    } else {
        colnames(df_)[3:ncol(df_)] <- paste0(region, '_', colnames(df_[, 3:ncol(df_)]))
        df <- df_
    }
    print(dim(df_))
}

exrp <- df[df$group!='ODC', ]


pro <- 'CD47'
region <- 'Hipp'



selected_freq_col <- grep(region, colnames(freq))
freq_sel <- freq[, selected_freq_col]

selected_expr_col <- sapply(colnames(exrp), function(x) (length(grep(pro, x))>0) & (length(grep(region, x))>0))
expr_sel <- exrp[, selected_expr_col]


weighted_expr <- apply(freq_sel*expr_sel, 1, sum)


wilcox.test(weighted_expr[freq$group=='LowNo'], weighted_expr[freq$group=='PHAD'])


boxplot(weighted_expr[freq$group=='LowNo'], weighted_expr[freq$group=='PHAD'])
dev.off()