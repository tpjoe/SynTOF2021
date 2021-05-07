library(tidyr)
library(dplyr)
library(pvca)
library(data.table)
library(uwot)
library(sva)
library(Biobase)
library(ggplot2)
pvca_plt <- function(X, features, title='', threshold=0.7, getValue=FALSE){
    # Plot for pvca: showing the weighted average of variance from different
    # X: a p x n matrix where n columns are number of samples
    # features: factors of covariates, for example batch, sex, etc.
    # threshold: is how many % variance PCA has to explain.
    # getValue: if you want the returns to be plot values, otherwise return actual plots.
    
    pct_threshold <- threshold
    if (is.matrix(X) == FALSE){
        X <- as.matrix(X)
    }
    
    # converting inputs to an Expression Set object
    ExpressionSet()
    expressionSet <- ExpressionSet(assayData=X)
    for (i in seq_along(features)){
        pData(expressionSet)[i] <- as.data.frame(features[i])
    }
    batch.factors <- names(pData(expressionSet))

    # run PVCA (without progress report)
    invisible(capture.output(pvcaObj <- pvcaBatchAssess(expressionSet, batch.factors, 0.7)))

    # plot data
    values <- as.numeric(pvcaObj[[1]])
    label <- as.data.frame(pvcaObj[[2]])
    plot_data <- cbind(values, label)
    colnames(plot_data) <- c('Value', 'Effect')
    
    # plot
    if (getValue==FALSE) {
        p <- ggplot(data=plot_data, aes(x=Effect, y=Value)) +
        ylab('Weighted average proportion variance') +
        geom_bar(stat="identity") +
        geom_text(aes(label=sprintf("%0.2f", round(values, digits=2))), vjust=-0.3, size=5) +
        ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
        plot(p)
    }
    # values
    if (getValue == TRUE) {
        return(plot_data)
    }
}


### load batch key -------------------------------------------------------------------
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
mc_ <- mc[, c(3:6)]
batch_key <- mc_ %>% distinct()
demo <- as.data.frame(fread('../raw_data/demographics.csv'))


#### Importing data ------------------------------------------------------------------
df <- data.frame()
bk <- c()
sex <- c()
resv_colnames <- list()
regions <- rep(c('BA9', 'DLCau', 'Hipp'), each=(24-3)) #24 sampples, of which 3 were ODC
for (region in c('BA9', 'DLCau', 'Hipp')) {
    df_ <- read.csv(paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp2mc_5_13.csv'))[, -1]
    df_ <- df_[df_[, 1] != 'ODC', ]
    bk_ <- batch_key[batch_key$region==region, ]
    bk_ <- bk_[match(df_[, 2], bk_$sample), 'batch']
    sex_ <- demo[match(df_[, 2], paste0(gsub(' ', '', demo$sample), '.fcs')), 'sex']
    resv_colnames[[region]] <- paste(region, colnames(df_)[3:ncol(df_)], sep='_')
    colnames(df_) <- seq(dim(df_)[2])
    df <- rbind.data.frame(df, df_)
    bk <- c(bk, bk_)
    sex <- c(sex, sex_)
}



#### PVCA plots ------------------------------------------------------------------
group <- df[, c(1)]
batch <- factor(bk)
sex <- factor(sex)
X <- cbind.data.frame(batch, group, sex, df[, -c(1, 2)])
region <- 'DLCau'

pdf(paste0('figures/batch_effects/pvca_', region, '.pdf'))
plt <- pvca_plt(X=t(df[regions==region, -c(1, 2)]), features=data.frame(sex, group, batch)[regions==region, ], getValue=FALSE)
dev.off()

pdf(paste0('figures/batch_effects/pvca_', 'LowNo', '.pdf'))
plt <- pvca_plt(X=t(df[group=='LowNo', -c(1, 2)]), features=data.frame(sex, regions, batch)[group=='LowNo', ], getValue=FALSE)
dev.off()


#### UMAP plots ------------------------------------------------------------------
xy <- umap(as.matrix(df[, -c(1, 2)]), n_neighbors=15, min_dist=1)

pdf(paste0('figures/batch_effects/umap_col=batch_shape=regions.pdf'))
plot(xy[, 1], xy[, 2], pch=as.numeric(factor(regions)) , col=factor(batch), legned=TRUE)
legend("topright", inset=.02, title="Region", legend=unique(regions), col='black', 
        pch=unique(as.numeric(factor(regions))), horiz=TRUE, cex=0.8)
dev.off()



### Perform batch correction -----------------------------------------------------
BA9 <- ComBat(dat=t(df[regions=='BA9', -c(1, 2)]), batch=as.numeric(batch[regions=='BA9']))
DLCau <- ComBat(dat=t(df[regions=='DLCau', -c(1, 2)]), batch=as.numeric(batch[regions=='DLCau']))
Hipp <- ComBat(dat=t(df[regions=='Hipp', -c(1, 2)]), batch=as.numeric(batch[regions=='Hipp']))

merged_adjusted <- rbind.data.frame(t(BA9), t(DLCau), t(Hipp))
xy <- umap(aa, n_neighbors=15, min_dist=1)

pdf(paste0('figures/batch_effects/test.pdf'))
plot(xy[, 1], xy[, 2], pch=as.numeric(factor(regions)) , col=factor(batch), legned=TRUE)
dev.off()



### export
for (region in c('BA9', 'DLCau', 'Hipp')) {
    data_mat <- ComBat(dat=t(df[regions==region, -c(1, 2)]), batch=as.numeric(batch[regions==region]))
    export_mat <- cbind.data.frame(df[regions==region, c(1, 2)], t(data_mat))
    colnames(export_mat) <- c('group', 'sample', resv_colnames[[region]])
    write.csv(export_mat, paste0('R_py_exchange/df_meanAllMarkers_combat_', region, '_noStd_exp2_5,7,14.csv'))
}


# # limma with trying to preserve the region effect does not work
# df_limma <- t(limma::removeBatchEffect(t(df[, -c(1, 2)]), batch=batch, design=as.matrix(as.numeric(factor(regions)))))
# xy <- umap(df_limma, n_neighbors=15, min_dist=1)
# pdf(paste0('figures/batch_effects/test.pdf'))
# plot(xy[, 1], xy[, 2], pch=as.numeric(factor(regions)) , col=factor(batch), legned=TRUE)
# dev.off()