library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggsignif)
library(yarrr)

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
df <- df[df$group!='ODC', ]
levels(df[, 'group']) <- c('LBD', "Control", "ODC", "  ADNC" )
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
colnames(df) <- gsub('_', ' ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippocampus', colnames(df))
colnames(df) <- gsub('DLCau', 'Caudate', colnames(df))


# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (Pro in c('PHF-tau', 'CD47')) {
    # region <- c(rep('BA9', 15), rep('DLCau', 15), rep('Hippocampus', 15))
    region <- rep('Hippocampus', 15)
    selectedPro <- rep(Pro, length(region))
    # selectedFeature <- colnames(df)[grepl(Pro, colnames(df))]
    selectedFeature <- colnames(df)[grepl(Pro, colnames(df)) & grepl('Hippocampus', colnames(df))]

    df_sel <- df[, c('group', selectedFeature)]

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



# for MIBI correlations 
hipp_df <- df[, grepl('group|Hippocampus', colnames(df))]
corr_df <- hipp_df[, grepl('group|CD47|PHF-tau', colnames(hipp_df))]
corr_df <- corr_df %>% pivot_longer(!group, names_to=c('region', 'marker', 'cl'), names_sep=' ', values_to="value")

cor(corr_df[corr_df$group=='Control' & corr_df$marker=='PHF-tau', 'value']/corr_df[corr_df$group=='Control' & corr_df$marker=='CD47', 'value'],
    corr_df[corr_df$group=='Control' & corr_df$marker=='CD47', 'value'], method='spearman')


cor(corr_df[corr_df$group=='AD' & corr_df$marker=='PHF-tau', 'value']/corr_df[corr_df$group=='AD' & corr_df$marker=='CD47', 'value'],
    corr_df[corr_df$group=='AD' & corr_df$marker=='CD47', 'value'], method='spearman')


# for ones with reselient

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
# df <- df[df$group!='ODC', ]
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
colnames(df) <- gsub('_', ' ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippo', colnames(df))
colnames(df) <- gsub('DLCau', 'CN', colnames(df))

df <- df[(df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') | df$group=='LowNo'), ]
df['resilient'] <- df$group
levels(df$resilient)[levels(df$resilient)=="AD"] <- "resilient"
levels(df$resilient)[levels(df$resilient)=="LBD"] <- "resilient"
df$resilient <- factor(df$resilient, levels = c("LowNo", "resilient", "ODC"))


# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (Pro in c('ApoE')) {
    region <- c(rep('BA9', 15), rep('DLCau', 15), rep('Hipp', 15))
    selectedPro <- rep(Pro, length(region))
    selectedFeature <- colnames(df)[grepl(Pro, colnames(df))]

    df_sel <- df[, c('group', 'resilient', selectedFeature)]

    data <- df_sel %>% pivot_longer(-c(group, resilient), names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('LowNo', 'AD', 'LBD'))#, 'ODC'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    my_comparisons <- list(c('LowNo', 'resilient'))#, c('LowNo', 'LBD'), c('LowNo', 'ODC'))
    p <- ggboxplot(data, x='resilient', y='Expression', color='black', facet.by='Marker', outlier.shape=NA,
            short.panel.labs=TRUE, scales="free", strip.position="top", fill="resilient", ncol=5) + #, add="jitter") +
            geom_jitter(aes(colour=group), size=5.5, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression")  +
            scale_color_manual(values=c("black", "black", "black"))  +
            scale_fill_manual(values=c("#341069", "#d6446d", "#fea16e"))  + #, 'orange')) + 
            theme(legend.position="top", legend.title = element_blank()) +
            stat_compare_means(
                mapping=aes(label=format.pval(..p.adj.., digits=1)),
                method="wilcox.test", ref.group="LowNo", cex=3,
                label.y.npc=0.9) +
            theme(plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"))
    ggsave(filename=paste0('figures/boxplots/expression_', Pro, '_noStd_exp2mc_5_13', '.pdf'), p, width=10, height=16)
}









#### Importing post data ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
    # df_ <- read.csv(paste0('R_py_exchange/df_medianAllMarkers_', region, '_noStd_exp0.csv'))[, -1]
    df_ <- read.csv(paste0('R_py_exchange/df_meanPostMarkers_', region, '_noStd.csv'))[, -1]
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

df$group <- factor(df$group, labels=c('LBD', 'LowNo', 'AD'))
splt_name <- sapply(colnames(df_), function(x) strsplit(x, '_')[[1]])
protein <- unique(sapply(splt_name, function(x) paste0(x[2:(length(x)-2)], collapse='_'))[-c(1, 2)])

Pro <- 'DAT'
region <- c(rep('BA9', 13), rep('DLCau', 15), rep('Hipp', 14))
selectedPro <- rep(Pro, length(region))
selectedFeature <- colnames(df)[grepl(Pro, colnames(df))]

df_sel <- df[, c('group', selectedFeature)]

data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
data <- as.data.frame(data)
data$group <- ordered(factor(data$group), levels=c('LowNo', 'LBD', 'AD'))
data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)


# faceted boxplots
my_comparisons <- list(c('LowNo', 'LBD'), c('LowNo', 'AD'))
p <- ggboxplot(data, x='group', y='Expression', color='group', facet.by='Marker', outlier.shape = NA,
        short.panel.labs=TRUE, scales="free", strip.position="top", add="jitter", ncol=8) +
        theme_classic() + labs(x="", y="Mean Expression") +
        scale_color_manual(values=c("#485696", "#E0A458", "#FC5130")) + 
        theme(legend.position="top", legend.title = element_blank()) +
        stat_compare_means(
            mapping=aes(label=format.pval(..p.adj.., digits=1)),
            method = "wilcox.test", ref.group="LowNo",
            label.y.npc=0.9) +
        theme(plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"))
ggsave(filename=paste0('figures/boxplots/expressionPost_', Pro, '_noStd', '.pdf'), p, width=15, height=15)











#### Plotting single feature ------------------------------------------------------------------
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


# For AD
desired_feature <- c(
                    #  'Hippocampus, PHF-tau, A1', 'Hippocampus, PHF-tau, B1', 'Hippocampus, PHF-tau, C1'
                    #  'Hippocampus, CD47, A2', 'Hippocampus, CD47, B2', 'Hippocampus, CD47, C1'
                    #  'BA9, DJ1, A1', 'BA9, DJ1, B2', 'BA9, DJ1, C7'
                    # 'Hippocampus, PHF-tau, C7', 'BA9, DJ1, C10'
                    'Hippocampus, CD47, C3'
                     )


# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
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


# For AD with resilient
desired_feature <- c(
                        'Hippocampus, ApoE, A1', 'Hippocampus, ApoE, C5'
                    #  'Caudate, 3-NT, B1', 'Caudate, 3-NT, C5', 'Caudate, 3-NT, C7',
                    #  'BA9, CD56, A1', 'BA9, CD56, A2', 'BA9, CD56, C4', 
                    #  'BA9, DJ1, A1', 'BA9, DJ1, B1', 'BA9, DJ1, C7', 'BA9, DJ1, C10',
                    #  'BA9, GAMT, B2', 'BA9, GAMT, C1', 'BA9, GAMT, C7',
                    #  'Hippocampus, GBA1, C1', 'Hippocampus, GBA1, C7',
                    #  'BA9, GFAP, A1', 'BA9, GFAP, C3', 'BA9, GFAP, C4'
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



# For LBD
desired_feature <- c('BA9, GFAP, C3', 'BA9, CD47, A2', 'CN, DAT, A2', 'BA9, GAMT, C11', 'BA9, Calreticulin, C1',
                     'CN, DAT, A1', 'BA9, GFAP, A1', 'BA9, EAAT1, A1')


# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- df[, c('group', selectedFeature)]

    data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('LowNo', 'AD', 'LBD'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    p <- ggboxplot(data, x='group', y='Expression', color='group', outlier.shape=NA,
            short.panel.labs=TRUE, scales="free", strip.position="top") +
            geom_jitter(aes(colour=group), size=4, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression") +
            scale_color_manual(values=c("black", "#fa164f", "#2a18f2"))  +
            theme(legend.position="none", plot.title=element_text(hjust=0.5, size=22), text=element_text(size=20.5), axis.text=element_text(size=19)) +
            stat_compare_means(
                mapping=aes(label=format.pval(..p.adj.., digits=1)),
                method="wilcox.test", ref.group="LowNo", cex=8,
                label.y.npc=0.9) +
            theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 
    p <- p + labs(title=selectedFeature)
    ggsave(filename=paste0('figures/boxplots/LBD/singleExpression_', selectedFeature, '.pdf'), p, width=4, height=6)
}








#### Plotting single feature for GFAP-EAAT- ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
    # df_ <- read.csv(paste0('R_py_exchange/df_medianAllMarkers_', region, '_noStd_exp0.csv'))[, -1]
    df_ <- read.csv(paste0('R_py_exchange/df_meanAllMarkersGFAPnegEAAT1neg_', region, '_noStd_exp0mc_5_13.csv'))[, -1]
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
levels(df$group) <- c('LBD', 'Control', 'AD')


# For AD with resilient
desired_feature <- c('Hippocampus, ApoE, A1', 'Hippocampus, ApoE, C5'
                    #  'Caudate, 3-NT, B1', 'Caudate, 3-NT, C5', 'Caudate, 3-NT, C7',
                    #  'BA9, CD56, A1', 'BA9, CD56, A2', 'BA9, CD56, C4', 
                    #  'BA9, DJ1, A1', 'BA9, DJ1, B1', 'BA9, DJ1, C7', 'BA9, DJ1, C10',
                    #  'BA9, GAMT, B2', 'BA9, GAMT, C1', 'BA9, GAMT, C7',
                    #  'Hippocampus, GBA1, C1', 'Hippocampus, GBA1, C7',
                    #  'BA9, GFAP, A1', 'BA9, GFAP, C3', 'BA9, GFAP, C4'
                     )

df['group_resilience'] <- df$group
levels(df$group_resilience) <- c(levels(df$group_resilience), 'Permissive AD', 'Resilient AD') 
df[df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') &
   df$group=='AD', 'group_resilience'] <- 'Resilient AD'
df[!df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') &
   df$group=='AD', 'group_resilience'] <- 'Permissive AD'

# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- df[, c('group', 'group_resilience', selectedFeature)]
    df_sel <- df_sel[df_sel$group %in% c('Control', 'AD'), ]

    data <- df_sel %>% pivot_longer(-c(group, group_resilience), names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('Control', 'AD'))
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
    ggsave(filename=paste0('figures/boxplots/AD_resGFAPnegEAAT1neg/singleExpression_', selectedFeature, '.pdf'), p, width=5.5, height=6)
}







#### Plotting multiple feature for GFAP-EAAT- ------------------------------------------------------------------

#### Importing data ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
    # df_ <- read.csv(paste0('R_py_exchange/df_medianAllMarkers_', region, '_noStd_exp0.csv'))[, -1]
    df_ <- read.csv(paste0('R_py_exchange/df_meanAllMarkersGFAPnegEAAT1neg_', region, '_noStd_exp0mc_5_13.csv'))[, -1]
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
levels(df[, 'group']) <- c("LBD", "Control", "ODC", "AD" )
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
colnames(df) <- gsub('_', ' ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippocampus', colnames(df))
colnames(df) <- gsub('DLCau', 'Caudate', colnames(df))


# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (Pro in c('ApoE', 'GFAP', 'EAAT1')) {
    region <- c(rep('BA9', 15), rep('DLCau', 15), rep('Hippocampus', 15))
    # region <- rep('Hippocampus', 15)
    selectedPro <- rep(Pro, length(region))
    selectedFeature <- colnames(df)[grepl(Pro, colnames(df))]
    # selectedFeature <- colnames(df)[grepl(Pro, colnames(df)) & grepl('Hippocampus', colnames(df))]

    df_sel <- df[, c('group', selectedFeature)]

    data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    data$group <- ordered(factor(data$group), levels=c('Control', 'AD', 'LBD', 'ODC'))
    data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    my_comparisons <- list(c('Control', 'AD'), c('Control', 'LBD'), c('Control', 'ODC'))
    p <- ggboxplot(data, x='group', y='Expression', color='group', facet.by='Marker', outlier.shape=NA, fill='group', alpha=0.2,
            short.panel.labs=TRUE, scales="free", strip.position="top", ncol=5) +
            geom_jitter(aes(colour=group), size=3, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression") +
            scale_color_manual(values=c("black", "#fa164f", "#2a18f2", 'orange')) + 
            scale_fill_manual(values=c("black", "#fa164f", "#2a18f2", 'orange')) + 
            theme(legend.position="top", legend.title = element_blank()) +
            stat_compare_means(
                mapping=aes(label=format.pval(..p.., digits=1)),
                method="wilcox.test", ref.group="Control", cex=4,
                label.y.npc=0.9) +
            theme(plot.margin = unit(c(0, 0.2, 0, 0.2), "cm"))
    ggsave(filename=paste0('figures/boxplots/AD_resGFAPnegEAAT1neg/expression_', Pro, '_noStd_exp0mc_5_13', '.pdf'), p, width=10.5, height=16)
    # ggsave(filename=paste0('figures/boxplots/All/expression_', Pro, '_noStd_exp2mc_5_13', '.pdf'), p, width=10.5, height=6)
}






#### Plotting single feature for other ground truths ------------------------------------------------------------------

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



gt <- read.csv('../raw_data/demographics.csv')
gt['sample_ID'] <- gsub(' ', '', gt$sample_ID)
gt['sample_ID'] <- paste0(gt$sample_ID, '.fcs')
gt <- gt[gt$sample_ID %in% as.character(df$sample), ]

gt <- gt[match(df$sample, gt$sample_ID), ]


# For Dementia
colnames(gt)
sel_y <- gt[, 'CognitiveStatus']

sel_y[sel_y=='MCI'] <- 'Dementia'
sel_y <- factor(sel_y, level=rev(c('Dementia', 'No dementia')), ordered=TRUE)

desired_feature <- c('BA9, GAMT, C7', 'Hippo, ApoE, A1', 'BA9, DJ1, C6', 'CN, AS, C6', 'Hippo, Parkin, A2', 'Hippo, CD47, C1')

# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- cbind.data.frame(df[, c(selectedFeature, 'group')], sel_y)
    colnames(df_sel) <- c(selectedFeature, 'diagnosis', 'group')

    data <- df_sel %>% pivot_longer(-c(group, diagnosis), names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    # data$group <- ordered(factor(data$group), levels=c('LowNo', 'AD', 'LBD'))
    # data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)

    # faceted boxplots
    p <- ggboxplot(data, x='group', y='Expression', color='#4c4b4d', outlier.shape=NA,
            short.panel.labs=TRUE, scales="free", strip.position="top") +
            geom_jitter(aes(colour=diagnosis), size=4, position=position_jitter(0.15)) +
            theme_classic() + labs(x="", y="Mean Expression") +
            scale_color_manual(values=c('#2a18f2', 'black', "#fac116", "#fa164f"))  +
            theme(legend.position="none", plot.title=element_text(hjust=0.5, size=22), text=element_text(size=20.5), 
            axis.text.x=element_text(angle=30, vjust=0.75), axis.text=element_text(size=19)) +
            stat_compare_means(method="wilcox.test", cex=6) + #, ref.group="No Dementia", label.y.npc=0.2,  mapping=aes(label=format.pval(..p.adj.., digits=1)),
            theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) 
    p <- p + labs(title=selectedFeature)
    ggsave(filename=paste0('figures/boxplots/Dementia/singleExpression_', selectedFeature, '.pdf'), p, width=4, height=6)
}









#### Plotting single feature for other regression ground truths ------------------------------------------------------------------

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



gt <- read.csv('../raw_data/demographics.csv')
gt['sample_ID'] <- gsub(' ', '', gt$sample_ID)
gt['sample_ID'] <- paste0(gt$sample_ID, '.fcs')
gt <- gt[gt$sample_ID %in% as.character(df$sample), ]

gt <- gt[match(df$sample, gt$sample_ID), ]


# For Dementia
colnames(gt)
sel_y <- gt[, 'calc_braakStage'] #'calc_NIA_AA', 'calc_thalPhase', 'calc_braakStage'

desired_feature <- c('Hippo, PHF-tau, B1')

# for (Pro in c('p.Tau', 'K48.Ubiquitin', 'CD47', 'DAT', 'a.Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')) {
for (selectedFeature in desired_feature) {
    print(selectedFeature)
    df_sel <- cbind.data.frame(df[, c(selectedFeature, 'group', 'sample')], sel_y)
    colnames(df_sel) <- c(selectedFeature, 'diagnosis', 'sample', 'group')

    data <- df_sel %>% pivot_longer(-c(group, diagnosis, sample), names_to="Marker", values_to="Expression")
    data <- as.data.frame(data)
    # faceted boxplots
    
    p <- ggscatter(data, x="group", y="Expression", fill='diagnosis', size=4.5, color=transparent('white', trans.val=1), 
        add="reg.line", add.params=list(color='#4c4b4d', fill="lightgray"), # Customize reg. line
        conf.int=TRUE, shape=21,
        cor.coef=TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method="spearman", label.sep="\n", size=5),
        palette=c("#2a18f2", "black", '#fac116', '#fa164f')) +
        theme_classic() + theme(text=element_text(size=14), axis.text=element_text(size=14), legend.position="none", plot.title = element_text(hjust = 0.5)) + 
        ylab('Mean Expression') + xlab('Braak Stage') + labs(title=selectedFeature)

    # ggsave(filename=paste0('figures/scatter/test.pdf'), p, width=4, height=4)
    # ggsave(filename=paste0('figures/scatter/Braak_singleExpression_', selectedFeature,'.pdf'), p, width=4, height=4)#, width=20, height=15)
}
