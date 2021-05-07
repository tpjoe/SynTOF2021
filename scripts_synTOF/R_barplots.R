library(ggplot2)
library(tidyr)

#### Importing data ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
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


pair <- c('LowNo', 'PHAD')
X <- df[df$group %in% pair, !colnames(df) %in% c('group', 'sample')]
y <- df[df$group %in% pair, 'group']
y <- sapply(y, function(x) if (x==pair[1]) {0} else {1})

corrY <- cor(X, y, method='spearman')
corrY[is.na(corrY)] <- 0
zY <- sqrt((nrow(X)-3)/1.06)*atanh(corrY)
corrYP <- 2*pnorm(-abs(zY))

# corrYP <- data.frame()
# for (i in 1:ncol(X)){
#     corrYP <- rbind.data.frame(corrYP, wilcox.test(X[y==0, i], X[y==1, i])$p.value)
# }
# rownames(corrYP) <- colnames(X)
# colnames(corrYP) <- 'logP'

logP <- -log10(corrYP)


logP <- cbind.data.frame(feature=rownames(logP), logP)
logP <- logP[order(logP$logP, decreasing=FALSE), ]
logP$feature <- factor(logP$feature, levels=logP$feature)
region <- logP %>% separate(feature, c("region"), sep='_')
logP <- cbind.data.frame(logP, 'region'=factor(region[[1]], levels=c('BA9', 'DLCau', 'Hipp')))


line <- -log10(0.05)
bonLine <- -log10(0.05/length(corrY))

p <- ggplot(data=logP[(nrow(logP)-50):nrow(logP), ], aes(x=feature, y=logP, fill=region)) +
        geom_bar(stat="identity", alpha=1)  + 
        coord_flip() + labs(y='-log P value', x='')  +
        theme_bw() + scale_fill_manual(values=c("#605856", "#1C6E8C", "#FF4A1C")) +
        theme(legend.position="top", legend.title = element_blank()) + 
        geom_hline(yintercept=bonLine, linetype="dashed", color = "black", size=0.5) + 
        geom_hline(yintercept=line, color = "black", size=0.5) +
        annotate("text", x=10, y=bonLine+0.25, label="Bonferroni P value", angle=-90)

ggsave(paste0('figures/barplots/spearman_', pair[1], pair[2],'_top50_mean_noStd_exp2mc_5_13.pdf'), p)

