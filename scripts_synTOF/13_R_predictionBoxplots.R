"""
This script plots out prediction values for each disease group.
"""

library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggsignif)
library(yarrr)


#### boxplot AD ------------------------------------------------------------------
data <- read.csv('R_py_exchange_afterCluster/preds_PHAD_EN_noPerfect.csv')[, -1]
colnames(data) <- c('Predicted', 'Diagnosis')
data[, 2] <- sapply(data[, 2], function(x) if(x==0) {'Control'} else {'ADNC'})

p <- ggboxplot(data, x='Diagnosis', y='Predicted', color='#4c4b4d', outlier.shape=NA, fill='Diagnosis') + 
        geom_jitter(size=3.4, position=position_jitter(0.15)) + #aes(colour=Diagnosis),
        theme_classic() + labs(x="", y='EN Predicted Values') +
      #   scale_color_manual(values=c("black", "black")) +
        scale_fill_manual(values=c("#003049", "#F77F00")) +
        theme(legend.position="none", legend.title=element_blank(), text=element_text(size=16), 
        axis.text=element_text(size=16)) +
        stat_compare_means(cex=6,
            mapping=aes(label=paste0('p = ', format.pval(..p.., digits=2))),
            method="wilcox.test") + guides(colour=FALSE)+ ylim(c(0, 1))

ggsave(filename=paste0('figures/prediction_values/AD.pdf'), p, width=2.5, height=4)#, width=20, height=15)



#### boxplot LBD ------------------------------------------------------------------
data <- read.csv('R_py_exchange_afterCluster/preds_LBD_Ridge_noPerfect.csv')[, -1]
colnames(data) <- c('Predicted', 'Diagnosis')
data[, 2] <- sapply(data[, 2], function(x) if(x==0) {'Control'} else {'LB'})

p <- ggboxplot(data, x='Diagnosis', y='Predicted', color='#4c4b4d', outlier.shape=NA, fill='Diagnosis') + 
        geom_jitter(aes(colour=Diagnosis), size=3.4, position=position_jitter(0.15)) +
        theme_classic() + labs(x="", y='Ridge Predicted Values') +
        scale_color_manual(values=c("black", "black")) +
        scale_fill_manual(values=c("#003049", "#FCBF49")) +
        theme(legend.position="none", legend.title=element_blank(), text=element_text(size=16), 
        axis.text=element_text(size=16)) +
        stat_compare_means(cex=6,
            mapping=aes(label=paste0('p = ', format.pval(..p.., digits=2))),
            method="wilcox.test") + guides(colour=FALSE)+ ylim(c(0, 1))

ggsave(filename=paste0('figures/prediction_values/LBD.pdf'), p, width=2.5, height=4)#, width=20, height=15)