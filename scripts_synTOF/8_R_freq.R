"""
Plotting out frequencies of the resulting clusters
"""


library(ggplot2)
library(ggsignif)
library(ggpubr)
library(tidyr)

#### Frequency difference analysis ----------------------------------------------------------
# Importing data
pair <- c('LBD', 'PHAD')

pValue <- matrix(NA, 16, 3) #15-16 clusters in all 3 regions
colnames(pValue) <- c('BA9', 'DLCau', 'Hipp')
for (region in colnames(pValue)){
  df <- read.csv(paste0('R_py_exchange/df_freq_', region, '.csv'))[, -1]
  for (i in seq(ncol(df)-2)) {
    pValue[i, region] <- wilcox.test(df[df$group == pair[1], i+2], df[df$group == pair[2], i+2])$p.value
  }
}
write.csv(pValue, paste0('tables/freq_pValue/', pair[1], '_', pair[2], '.csv'))



#### Frequency plots ----------------------------------------------------------
# Importing data
for (region in c('BA9', 'DLCau', 'Hipp')){
    if (region=='BA9') {reg<-'BA9'} else if (region=='Hipp') {reg<-'Hippocampus'} else {reg<-'Caudate'}
    df <- read.csv(paste0('R_py_exchange/df_freq_', region, '_noStd.csv'))[, -1]
    # colnames(df)[3:ncol(df)] <- sapply(colnames(df)[3:ncol(df)], function(x) strsplit(x, 'X')[[1]][2])
    data <- df[, -2] %>% pivot_longer(-group, names_to="Cluster", values_to="Frequency")
    data <- as.data.frame(data)
    data <- data[data$group!='ODC', ]
    levels(data[, 'group']) <- c("LB", "Control", "ODC", "ADNC" )
    data$group <- droplevels(data$group)
    

    # grouped boxplot
    p <- ggplot(data, aes(x=ordered(factor(Cluster), levels=unique(data$Cluster)), y=Frequency,
                        fill=ordered(factor(group), levels=c('Control', 'ADNC', 'LB')))) + 
        geom_boxplot() + theme_classic() + labs(x="Cluster") + 
        scale_color_manual(values=c("black", "black", "black")) +
        scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e")) + 
        theme(legend.position="top", legend.title = element_blank(), legend.text=element_text(size=14), 
              axis.text=element_text(size=14), axis.title=element_text(size=14)) +
        ggtitle(region)
    ggsave(filename=paste0('figures/boxplots/frequency_', region, '.pdf'), p)

    # faceted boxplots
    data$group <- ordered(factor(data$group), levels=c('Control', 'ADNC', 'LB'))
    data$Cluster <- (factor(data$Cluster, levels=unique(data$Cluster), ordered=TRUE))
    my_comparisons <- list(c("Control", "LB"), c("Control", "ADNC"))
    p <- ggboxplot(data, x='group', y='Frequency', color='black', fill='group', facet.by='Cluster', 
                        outlier.shape=NA, short.panel.labs=FALSE, nrow=3, alpha=1) +
            geom_jitter(size=3.5, position=position_jitter(0.15)) + #aes(colour=group)
            theme_classic() + labs(x="", y=paste0('Frequency (', reg, ' Region)')) + ylim(0, 0.30) + 
            # scale_color_manual(values=c("black", "black", "black")) +
            scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e"))  + 
            theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=14), 
              axis.text.y=element_text(size=14), axis.text.x=element_text(size=14, angle=30, hjust=1),
              axis.title=element_text(size=14), strip.text=element_text(size=14)) + 
            
            stat_compare_means(ref.group='Control', label="p.signif", method="wilcox.test", 
            hide.ns=TRUE, label.y=0.23, cex=8)+ theme(legend.position="none") #comparisons=my_comparisons
    # if (reg!='BA9') {p <- p + theme(legend.position="none")}
    ggsave(filename=paste0('figures/boxplots/frequencyFaceted_', region, '_noStd.pdf'), p)
}
# dev.off()


#### Frequency plots of post synaptic ----------------------------------------------------------
# Importing data
for (region in c('BA9', 'DLCau', 'Hipp')){
    if (region=='BA9') {reg<-'BA9'} else if (region=='Hipp') {reg<-'Hippocampus'} else {reg<-'Caudate'}
    df <- read.csv(paste0('R_py_exchange/df_freqPost_', region, '_noStd.csv'))[, -1]
    df[is.na(df)] <- 0
    # colnames(df)[3:ncol(df)] <- sapply(colnames(df)[3:ncol(df)], function(x) strsplit(x, 'X')[[1]][2])
    
    # start pivoting
    data <- df[, -2] %>% pivot_longer(-group, names_to="Cluster", values_to="Frequency")
    data <- as.data.frame(data)
    data <- data[data$group!='ODC', ]
    levels(data[, 'group']) <- c("LB", "Control", "ODC", "ADNC" )
    data$group <- droplevels(data$group)
    
    # grouped boxplot
    p <- ggplot(data, aes(x=ordered(factor(Cluster), levels=unique(Cluster)), y=Frequency, 
                         fill=ordered(factor(group), levels=c('Control', 'ADNC', 'LB')),
                         color=ordered(factor(group), levels=c('Control', 'ADNC', 'LB')))) + 
        geom_boxplot(alpha=1) + theme_classic() + labs(x="Cluster", y=paste0('Frequency (', reg, ' Region)')) +
        # geom_jitter(aes(colour=group), size=2, position=position_jitter(0.15)) + 
        scale_color_manual(values=c("black", "black", "black")) +
        scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e"))  + 
        theme(legend.position="top", legend.title = element_blank(), legend.text=element_text(size=16), 
              axis.text=element_text(size=16), axis.title=element_text(size=16)) + theme(legend.position="none")
    ggsave(filename=paste0('figures/boxplots/frequencyPost_', region, '.pdf'), p)
}