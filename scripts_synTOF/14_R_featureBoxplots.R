"""
This scripts include main boxplots necessary for the manuscripts.
"""

# import libraries
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggsignif)
library(yarrr)



"""
For plotting all clusters of the selected marker in a brain region
"""
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

# reclassifying some LowNo to ODC
df <- df[df$group!='ODC', ]
levels(df[, 'group']) <- c('LBD', "Control", "ODC", "  ADNC" )
df$group <- droplevels(df$group)

# get marker names
splt_name <- sapply(colnames(df_), function(x) strsplit(x, '_')[[1]])
protein <- unique(sapply(splt_name, function(x) paste0(x[2:(length(x)-2)], collapse='_'))[-c(1, 2)])

# rename markers
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
          'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', 'p\\.Tau', 'X3NT')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', 'PHF-tau', '3-NT')
for (s in 1:length(to_sub)) {
    colnames(df) <- gsub(to_sub[s], sub_to[s], colnames(df))
}

# rename columns a bit
colnames(df) <- gsub('mean_', '', colnames(df))
colnames(df) <- gsub('_', ' ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippocampus', colnames(df))
colnames(df) <- gsub('DLCau', 'Caudate', colnames(df))

# selection of markers and brain regions
marker_list <- c('PHF-tau', 'CD47')
brain_region <- 'Hippocampus'

# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (Pro in marker_list) {
    # region <- c(rep('BA9', 15), rep('DLCau', 15), rep('Hippocampus', 15))
    region <- rep(brain_region, 15) #15 for 15 clusters
    selectedPro <- rep(Pro, length(region))
    # selectedFeature <- colnames(df)[grepl(Pro, colnames(df))]
    selectedFeature <- colnames(df)[grepl(Pro, colnames(df)) & grepl('Hippocampus', colnames(df))]

    df_sel <- df[, c('group', selectedFeature)]

    # long pivot the table for plotting
    data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('Control', '  ADNC', 'LBD', 'ODC'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    my_comparisons <- list(c('Control', '  ADNC'), c('Control', 'LBD'), c('Control', 'ODC'))
    p <- ggboxplot(data, x='group', y='Expression', color='group', facet.by='Marker', outlier.shape=NA, fill='group', alpha=1,
            short.panel.labs=TRUE, scales="free", strip.position="top", ncol=5) +
            geom_jitter(aes(colour=group), size=3, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression")  +
            scale_color_manual(values=c("black", "black", "black", "black"))  +
            scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e", "black"))  +
            # scale_color_manual(values=c("black", "#fa164f", "#2a18f2", 'orange')) + 
            # scale_fill_manual(values=c("black", "#fa164f", "#2a18f2", 'orange')) + 
            theme(legend.position="top", legend.title=element_blank(), #,
                  text=element_text(size=11), axis.text=element_text(size=11)) +
            stat_compare_means(
                mapping=aes(label=format.pval(..p.., digits=1)),
                method="wilcox.test", ref.group="Control", cex=4,
                label.y.npc=0.9) +
            theme(plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"))
    # ggsave(filename=paste0('figures/boxplots/All/expression_', Pro, '_noStd_exp2mc_5_13', '.pdf'), p, width=10.5, height=16)
    ggsave(filename=paste0('figures/boxplots/All/expression_', Pro, '_noStd_exp2mc_5_13', '.pdf'), p, width=10.5, height=6)
}



"""
For plotting individual clusters of the selected marker in a brain region
"""

#### Load data ------------------------------------------------------------------
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
df <- df[df$group!='ODC', ]
levels(df[, 'group']) <- c("LBD", "LowNo", "ODC", "AD" )
df$group <- droplevels(df$group)


splt_name <- sapply(colnames(df_), function(x) strsplit(x, '_')[[1]])
protein <- unique(sapply(splt_name, function(x) paste0(x[2:(length(x)-2)], collapse='_'))[-c(1, 2)])


# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
          'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', 'p\\.Tau', 'X3NT')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', 'PHF-tau', '3-NT')
for (s in 1:length(to_sub)) {
    colnames(df) <- gsub(to_sub[s], sub_to[s], colnames(df))
}
colnames(df) <- gsub('mean_', '', colnames(df))
colnames(df) <- gsub('_', ', ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippo', colnames(df))
colnames(df) <- gsub('DLCau', 'CN', colnames(df))
print(table(sapply(colnames(df), function(x) strsplit(x, ', ')[[1]][2])))

# change hippo to hippocampus and class name
colnames(df) <- gsub('Hippo', 'Hippocampus', colnames(df))
colnames(df) <- gsub('CN,', 'Caudate,', colnames(df))
levels(df$group) <- c('LBD', 'Control', '  ADNC')


##########################
# Fetures showing for AD #
##########################
desired_feature <- c(
                    #  'Hippocampus, PHF-tau, A1', 'Hippocampus, PHF-tau, B1', 'Hippocampus, PHF-tau, C1'
                    #  'Hippocampus, CD47, A2', 'Hippocampus, CD47, B2', 'Hippocampus, CD47, C1'
                    #  'BA9, DJ1, A1', 'BA9, DJ1, B2', 'BA9, DJ1, C7'
                    # 'Hippocampus, PHF-tau, C7', 'BA9, DJ1, C10'
                    'Hippocampus, CD47, C3'
                     )

for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- df[, c('group', selectedFeature)]

    data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('Control', '  ADNC', 'LBD'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    p <- ggboxplot(data, x='group', y='Expression', color='group', outlier.shape=NA, fill='group', alpha=1,
            short.panel.labs=TRUE, scales="free", strip.position="top") +
            geom_jitter(aes(colour=group), size=6, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression") +
            scale_color_manual(values=c("black", "black", "black"))  +
            scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e"))  +
            theme(legend.position="none", plot.title=element_text(hjust=0.5, size=22), text=element_text(size=22), axis.text=element_text(size=22)) +
            stat_compare_means(
                mapping=aes(label=format.pval(..p.., digits=1)),
                method="wilcox.test", ref.group="Control", cex=9, size=11,
                label.y.npc=0.95) + ylim(0, 3) +
            theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 
    p <- p + labs(title=selectedFeature)
    ggsave(filename=paste0('figures/boxplots/AD/ScaledsingleExpression_', selectedFeature, '.pdf'), p, width=4, height=6)
}

##########################
# Fetures showing for AD resilient #
##########################
desired_feature <- c(
                        'Hippocampus, ApoE, A1', 'Hippocampus, ApoE, C5'
                     )

df['group_resilience'] <- df$group
levels(df$group_resilience) <- c(levels(df$group_resilience), 'Permissive AD', 'Resilient AD') 
df[df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') &
   df$group=='  ADNC', 'group_resilience'] <- 'Resilient AD'
df[!df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') &
   df$group=='  ADNC', 'group_resilience'] <- 'Permissive AD'

# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- df[, c('group', 'group_resilience', selectedFeature)]
    df_sel <- df_sel[df_sel$group %in% c('Control', '  ADNC'), ]

    data <- df_sel %>% pivot_longer(-c(group, group_resilience), names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('Control', '  ADNC'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    my_comparisons <- list(c('Control', 'Permissive AD'))
    p <- ggboxplot(data, x='group_resilience', y='Expression', color='group_resilience', outlier.shape=NA, 
            fill='group_resilience', alpha=1,
            short.panel.labs=TRUE, scales="free", strip.position="top") +
            geom_jitter(aes(colour=group_resilience), size=7, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression") +
            scale_color_manual(values=c("black", "black", "black"))  +
            scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e"))  +
            theme(legend.position="none", plot.title=element_text(hjust=0.5, size=20), text=element_text(size=20.5), 
                  axis.text=element_text(size=25), axis.text.x=element_text(size=20)) +
            stat_compare_means(
                comparisons = my_comparisons,
                mapping=aes(label=format.pval(..p.., digits=1)),
                method="wilcox.test", cex=11, size=8, 
                label.y.npc=0.57) +
            theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 
    p <- p + labs(title=selectedFeature) + scale_x_discrete(labels=c('Permissive AD'='Permissive\nAD', 'Resilient AD'='Resilient\nAD'))
    ggsave(filename=paste0('figures/boxplots/AD_res/singleExpression_', selectedFeature, '.pdf'), p, width=5.5, height=6)
}


