# libraries
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(segmented)

# import data
mibi <- read.csv('../raw_data/NormTifData.csv')
# mibi_norm <- mibi
# mibi_nonorm <- mibi
mibi <- mibi_nonorm
mibi <- mibi[, 2:ncol(mibi)]

# screen for rough typr of excitory/non-excitory synaposes
#| (mibi.CD56>0.004) (mibi.TotalTau>0.005) (mibi.VGLUT2<0.002)
# excite_condition_ad <- ((mibi$ApoE>0.004) | (mibi$VGLUT1>0.004) | (mibi$TotalTau>0.005)  | (mibi$CD56>0.004))
# excite_condition_hc <- ((mibi$ApoE>0.004) | (mibi$VGLUT1>0.004) | (mibi$TotalTau>0.015)  | (mibi$CD56>0.004))
# excite_condition_ad <- ((mibi$ApoE>0.004) | (mibi$VGLUT1>0.004)  | (mibi$VGLUT2>0.004) | (mibi$CD56>0.004) | (mibi$TotalTau>0.04))
# excite_condition_hc <- ((mibi$ApoE>0.004) | (mibi$VGLUT1>0.004)  | (mibi$VGLUT2>0.004) | (mibi$CD56>0.004) | (mibi$TotalTau>0.04))

# inspect individual by hist
hist(mibi[mibi$VGLUT1>0.004, 'VGLUT1'], breaks=100)
dev.off()

# inspect individual by zoom in cumulative
for (marker in c('PanAPOE')) { #c('CD47', 'CD56', 'TotalTau', 'ApoE', 'VGLUT1', 'PanGAD')){
    p <- ggplot(mibi[mibi$Point %in% seq(7, 12) & (mibi$Synaptophysin>0 & mibi$PSD95==0), ], aes_string(marker)) + stat_ecdf(geom="step", pad=FALSE, size=1) + 
    theme_bw(base_size = 20) + ylab('Accumulated Frequency') + xlab(paste0('Normalized ', marker, ' Expression Value'))
    # Piecewise regression
    cumu <- ggplot_build(p)$data[[1]]
    cumu <- cumu[order(cumu$x), ]
    lin.mod <- lm(y~x, data=cumu)
    set.seed(420)
    xx <- cumu$x
    # segmented.mod <- segmented(lin.mod, seg.Z=~xx, psi=list(xx=c(7.553e-05)))
    # my.fitted <- fitted(segmented.mod)
    # my.model <- data.frame(Marker=cumu$x, Elevation=my.fitted)
    # p <- p + geom_line(data=my.model, aes(x=Marker, y=Elevation), colour="tomato", linetype='dashed', size=1.25) + theme(legend.position="none") 
    # p <- p + annotate("text", x=segmented.mod$psi[2]+0.5, y=1, label=paste0('Value = ', round(segmented.mod$psi[2], 2)), color='tomato', size=8)  
    ggsave(paste0('figures/lineplots/', marker, '_mibiNoNorm_cumulativeFreq_HC_presyn.pdf'), p)
}

# excite_condition_ad <- ((mibi$VGLUT1>0.7) | (mibi$ApoE>0.49) | (mibi$CD56>1.78) | (mibi$TotalTau>0.51) | (mibi$CD47>0.53)) # all in one

# CDcondition_ad <- ((mibi$CD56>1.78) | (mibi$CD47>0.53))
# CDcondition_hc <- ((mibi$CD56>1.78) | (mibi$CD47>0.53))
# excite_condition_ad <- (mibi$TotalTau>0.51) | (mibi$VGLUT1>0.7 & CDcondition_ad) | CDcondition_ad
# excite_condition_hc <- (mibi$TotalTau>0.51) | (mibi$VGLUT1>0.7 & CDcondition_hc) | CDcondition_hc

CDcondition_ad <- ((mibi$CD56>1.72) | (mibi$CD47>0.47))
CDcondition_hc <- ((mibi$CD56>1.82) | (mibi$CD47>0.44))
excite_condition_ad <- mibi$TotalTau>0.34 | CDcondition_ad #| CDcondition_ad #(mibi$VGLUT1>0.7 & ((mibi$CD56>0) | (mibi$CD47>0))) | CDcondition_ad #| (mibi$TotalTau>0.34) 
excite_condition_hc <- mibi$TotalTau>0.9 | CDcondition_hc #| CDcondition_hc #(mibi$VGLUT1>2.58 & ((mibi$CD56>0) | (mibi$CD47>0))) | CDcondition_hc #| (mibi$TotalTau>0.9) 

# CDcondition_ad <- ((mibi$CD56>0.004) | (mibi$CD47>0.005))
# CDcondition_hc <- ((mibi$CD56>0.004) | (mibi$CD47>0.004))
# excite_condition_ad <- (mibi$VGLUT1>0.004) | CDcondition_ad# | (mibi$ApoE>0.004) #| (mibi$TotalTau>0.0175) 
# excite_condition_hc <- (mibi$VGLUT1>0.004) | CDcondition_hc# | (mibi$ApoE>0.004) #| (mibi$TotalTau>0.0175) 


excite_syn_ad <- mibi[excite_condition_ad & (mibi$Point %in% c(1, 2, 3, 4, 5, 6)) & (mibi$Synaptophysin>0 & mibi$PSD95==0), ]
excite_syn_hc <- mibi[excite_condition_hc & (mibi$Point %in% c(7, 8, 10, 9, 11, 12)) & (mibi$Synaptophysin>0 & mibi$PSD95==0), ]
nonExcite_syn_ad <- mibi[!excite_condition_ad & (mibi$Point %in% c(1, 2, 3, 4, 5, 6)) & (mibi$Synaptophysin>0 & mibi$PSD95==0), ]
nonExcite_syn_hc <- mibi[!excite_condition_hc & (mibi$Point %in% c(7, 8, 10, 9, 11, 12)) & (mibi$Synaptophysin>0 & mibi$PSD95==0), ]


# plot sc level exite comparison in one
plt_df <- excite_syn_hc %>% mutate(A='A*') %>% full_join(nonExcite_syn_hc %>% mutate(A='Non-A*')) %>% mutate(Group='Control')
plt_df <- plt_df %>% full_join(excite_syn_ad %>% mutate(A='A*') %>% full_join(nonExcite_syn_ad %>% mutate(A='Non-A*')) %>% mutate(Group='AD'))
plt_df <- plt_df[, c('Point', 'CD47', 'PHF1Tau', 'ApoE', 'Abeta', 'PanGAD', 'TotalTau', 'VGLUT1', 'CD56', 'Group', 'A')]
data.frame(plt_df %>% group_by(Point, A, Group) %>% summarise(n=n()))
plt_df <- plt_df[rev(order(plt_df$A)), ]
plt_df[, 'Group'] <- factor(plt_df$Group, levels=c('Control', 'AD'))


# max_y <- 0.006# for AD
# max_y <- 0.15# for HC



for (sel_marker in c('ApoE', 'PanGAD', 'CD56')){ #'Abeta', 'CD47', 'PHF1Tau', 'TotalTau', 'VGLUT1'
    max_y <- plt_df %>% group_by(Group, A) %>% summarise(means=mean(.data[[sel_marker]]), sd=sd(.data[[sel_marker]])) %>% summarise(maxM = max(means)) %>% summarise(max_y = max(maxM))
    max_y <- 1.5*max_y
    if (sel_marker=='PHF1Tau') {marker_label <-'PHF-Tau'
    } else if (sel_marker=='PanGAD') {marker_label <-'GAD65'
    } else if (sel_marker=='VGLUT1') {marker_label <-'vGLUT'
    } else if (sel_marker=='TotalTau') {marker_label <-'Tau'
    } else if (sel_marker=='Abeta') {marker_label <-'Ab40'
    } else {marker_label <- sel_marker}
    p <- ggbarplot(as.data.frame(plt_df), x='A', y=sel_marker, line.color="gray", line.size=0.4, point.size=3, facet.by='Group',
                color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y[[1]])) + 
                stat_compare_means(method='wilcox', label.y=max_y[[1]]-0.1*max_y[[1]], mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3))), size=6) +
                xlab('') + ylab(paste0('Pixel-wise ', marker_label, ' Expression in AD')) + theme_classic(base_size=23) + theme(legend.position='none') #+
                #    facet_wrap(~Group, scales="free", ncol=2)
    ggsave(paste0('figures/barplots/mibi_all_', marker_label, '_byVGLUT.pdf'), p)
}




# plot sc level exite comparison
plt_df <- excite_syn_ad %>% mutate(A='A*') %>% full_join(nonExcite_syn_ad %>% mutate(A='Non-A*'))
plt_df <- excite_syn_hc %>% mutate(A='A*') %>% full_join(nonExcite_syn_hc %>% mutate(A='Non-A*'))
plt_df <- plt_df[, c('Point', 'CD47', 'PHF1Tau', 'ApoE', 'Abeta', 'PanGAD', 'TotalTau', 'VGLUT1', 'CD56', 'A')]
plt_df %>% group_by(Point, A) %>% summarise(n=n())
plt_df <- plt_df[rev(order(plt_df$A)), ]

max_y <- 0.006# for AD
# max_y <- 0.15# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='PHF1Tau', line.color="gray", line.size=0.4, point.size=3, facet.by='Point',
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise PHF-tau Expression in AD') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_PHFtau_ad_noTau.pdf', p)

max_y <- 0.015# for AD
# max_y <- 0.15# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='CD47', line.color="gray", line.size=0.4, point.size=3, 
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise PHF-tau Expression in AD') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_PHFtau_ad_noTau.pdf', p)

max_y <- 0.3# for AD
# max_y <- 0.12# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='CD47', line.color="gray", line.size=0.4, point.size=3, facet.by='Point',
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise CD47 Expression in AD') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_CD47_ad_noTau.pdf', p)

max_y <- 0.0025# for AD
# max_y <- 0.3# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='ApoE', line.color="gray", line.size=0.4, point.size=3, facet.by='Point',
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise ApoE Expression in AD') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_ApoE_ad_noTau.pdf', p)

# max_y <- 0.8# for AD
max_y <- 0.3# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='Abeta', line.color="gray", line.size=0.4, point.size=3, facet.by='Point',
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise Abeta Expression in Controls') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_Abeta_hc.pdf', p)

# max_y <- 0.8# for AD
max_y <- 0.2# for HC
p <- ggbarplot(as.data.frame(plt_df), x='A', y='PanGAD', line.color="gray", line.size=0.4, point.size=3, facet.by='Point',
               color='black', palette="jco", fill="A", add="mean_se", ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('Pixel-wise GAD65 Expression in Controls') + theme_classic(base_size=15) + theme(legend.position='none')
ggsave('figures/barplots/mibi_PanGAD_hc.pdf', p)




# ggsave('test.pdf',
# ggsave('test.pdf', p)
# ggsave('figures/pairedplots/PHFtau_ad.pdf', p)
# sapply(seq(5), function(i) wilcox.test(plt_df[(plt_df$Point==i) & (plt_df$excited=='Excited'), 'CD47'], plt_df[(plt_df$Point==i) & (plt_df$excited=='NonExcited'), 'CD47']))



# # plot sc level AD vs HC comparison
# # plt_df <- excite_syn_hc %>% mutate(group='Control') %>% full_join(excite_syn_ad %>% mutate(group='AD'))
# plt_df <- nonExcite_syn_hc %>% mutate(group='Control') %>% full_join(nonExcite_syn_ad %>% mutate(group='AD'))
# plt_df <- plt_df[, c('Point', 'CD47', 'PHF1Tau', 'group', 'TotalTau', 'CD56', 'ApoE', 'PanGAD')]
# plt_df %>% group_by(Point, group) %>% summarise(n=n())
# plt_df <- plt_df[rev(order(plt_df$group)), ]
# # plt_df[plt_df$Point>6, 'Point'] <- plt_df[plt_df$Point>6, 'Point'] - 6
# plt_df['Point'] <- as.factor(plt_df$Point)
# plt_df['group'] <- as.factor(plt_df$group)

# # max_y <- 0.06 # for nonExcited
# max_y <- 0.6 # for Excited
# p <- ggbarplot(as.data.frame(plt_df), x='group', y='PHF1Tau', line.color="gray", line.size=0.4, point.size=3, position=position_dodge(0.8),
#                color='black', fill="Point", add="mean_se", ylim=c(0, max_y)) + #(, palette="jco")
#                scale_fill_manual(values=rep(c('#ff3d6a', 'black'), each=6)) +
#                stat_compare_means(method='kruskal.test()', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3))), label.x=0.75, size=6) +
#                xlab('') + ylab('Normalized PHF-tau Expression\nin Non-A* Pixels') + theme_classic(base_size=20) + theme(legend.position="none")
# ggsave('figures/barplots/mibi_PHFtau_nonExcited.pdf', p)

# # max_y <- 0.015 # for nonExcited
# max_y <- 0.3 # for Excited
# p <- ggbarplot(as.data.frame(plt_df), x='group', y='CD47', line.color="gray", line.size=0.4, point.size=3, position=position_dodge(0.8),
#                color='black', fill="Point", add="mean_se", ylim=c(0, max_y)) + #(, palette="jco")
#                scale_fill_manual(values=rep(c('#ff3d6a', 'black'), each=6)) +
#             #    scale_color_manual(values=rep('black', 12)) +
#                stat_compare_means(method='aov', size=6, label.y=max_y-0.1*max_y) +
#                xlab('') + ylab('Normalized CD47 Expression\nin A* Pixels') + theme_classic(base_size=20) + theme(legend.position="none")
# ggsave('test.pdf', p)
# ggsave('figures/barplots/mibi_CD47_excited.pdf', p)

# # summary(lmerTest::lmer(CD47~group + (1|Point), plt_df))
# # summary(lm(TotalTau~group, plt_df))
# # anova(nlme::lme(CD47 ~ group, random=~1|Point, data=plt_df))

# max_y <- 0.5 # for nonExcited
# # max_y <- 0.4 # for Excited
# p <- ggbarplot(as.data.frame(plt_df), x='group', y='TotalTau', line.color="gray", line.size=0.4, point.size=3, position=position_dodge(0.8),
#                color='black', fill="Point", add="mean_se") + #(, palette="jco")
#                scale_fill_manual(values=rep(c('#ff3d6a', 'black'), each=6)) +
#             #    scale_color_manual(values=rep('black', 12)) +
#                stat_compare_means(method='kruskal.test()', label.y=max_y-0.1*max_y, mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3))), size=6)+
#                xlab('') + ylab('Normalized Tau Expression\nin Non-A* Pixels') + theme_classic(base_size=20) + theme(legend.position="none")
# ggsave('figures/barplots/mibi_Tau_nonExcited.pdf', p)




# frequency
count_excite_ad <- excite_syn_ad %>% group_by(Point) %>% summarise(n=n())
count_excite_hc <- excite_syn_hc %>% group_by(Point) %>% summarise(n=n())
count_nonExcite_ad <- nonExcite_syn_ad %>% group_by(Point) %>% summarise(n=n())
count_nonExcite_hc <- nonExcite_syn_hc %>% group_by(Point) %>% summarise(n=n())

freq_ad <- count_excite_ad/(count_excite_ad + count_nonExcite_ad)*100
mean(freq_ad[, 2])
freq_hc <- count_excite_hc/(count_excite_hc + count_nonExcite_hc)*100
mean(freq_hc[, 2])


freq_all <- freq_ad %>% mutate(Group='AD') %>% full_join(freq_hc %>% mutate(Group='Control'))
freq_all[, 'Group'] = factor(freq_all[, 'Group'], levels=c('Control', 'AD'))

p <- ggbarplot(freq_all, "Point", "n", fill="grey", color = "black", add = c("mean_se", "dotplot"), jitter.size=6, facet.by='Group') +
    ylab('Frequency (%) of\nSubpopulation A*') + xlab('') + theme(
    text = element_text(size=20),
    axis.text.x = element_blank())
    # ggtitle('Subpopulations A*')
ggsave('figures/pairedplots/freq_byTau.pdf', p, width=5, height=5)


# marker <- 'CD47'
# marker <- 'PHF1Tau'

# # marker expressions
# mean_excite_ad <- excite_syn_ad %>% group_by(Point) %>% summarise(excite_ad_CD47=mean(CD47), excite_ad_PHFtau=mean(PHF1Tau))
# mean_excite_hc <- excite_syn_hc %>% group_by(Point) %>% summarise(excite_hc_CD47=mean(CD47), excite_hc_PHFtau=mean(PHF1Tau))
# mean_nonExcite_ad <- nonExcite_syn_ad %>% group_by(Point) %>% summarise(nonExcite_ad_CD47=mean(CD47), nonExcite_ad_PHFtau=mean(PHF1Tau))
# mean_nonExcite_hc <- nonExcite_syn_hc %>% group_by(Point) %>% summarise(nonExcite_hc_CD47=mean(CD47), nonExcite_hc_PHFtau=mean(PHF1Tau))

# mean_ad <- mean_nonExcite_ad %>% merge(mean_excite_ad, by='Point')
# mean_hc <- mean_nonExcite_hc %>% merge(mean_excite_hc, by='Point')

# # merge mean of all
# mean_all <- mean_ad %>% full_join(mean_hc)
# mean_all_long <- mean_all %>% pivot_longer(!Point, names_to=c('excite', 'group', 'marker'), names_sep='_', values_to="count") %>% drop_na()


# # plot expression
# plt_df <- data.frame(mean_all_long[(mean_all_long$marker=='PHFtau') & (mean_all_long$group=='ad'), ])
# plt_df <- plt_df[rev(order(plt_df$excite)), ]
# p <- ggpaired(plt_df, x='excite', y='count', line.color="gray", line.size=0.4, point.size=3,
#                color='black', palette="jco", fill="excite") + stat_compare_means(method='wilcox', pair=TRUE) 
#             #    ylim(0, 2.5e-5)
# # ggsave('test.pdf', p)
# ggsave('figures/pairedplots/PHFtau_ad.pdf', p)


# t.test(c(mean_all_long[(mean_all_long$marker=='PHFtau') & (mean_all_long$group=='ad') & (mean_all_long$excite=='excite'), 'count'])[[1]], 
#             c(mean_all_long[(mean_all_long$marker=='PHFtau') & (mean_all_long$group=='ad') & (mean_all_long$excite=='nonExcite'), 'count'])[[1]],
#             paired=TRUE
#             )







# syntof A comparison
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
levels(df[, 'group']) <- c("LBD", "Control", "ODC", "AD" )
df$group <- droplevels(df$group)


splt_name <- sapply(colnames(df_), function(x) strsplit(x, '_')[[1]])
protein <- unique(sapply(splt_name, function(x) paste0(x[2:(length(x)-2)], collapse='_'))[-c(1, 2)])

# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
          'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', 'p\\.Tau', 'X3NT', 'VGLUT')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', 'PHF-tau', '3-NT', 'vGLUT')
for (s in 1:length(to_sub)) {
    colnames(df) <- gsub(to_sub[s], sub_to[s], colnames(df))
}

colnames(df) <- gsub('mean_', '', colnames(df))
colnames(df) <- gsub('_', ' ', colnames(df))
colnames(df) <- gsub('Hipp', 'Hippocampus', colnames(df))
colnames(df) <- gsub('DLCau', 'Caudate', colnames(df))

Pro <- 'vGLUT'
region <- rep('Hippocampus', 15)
selectedPro <- rep(Pro, length(region))
selectedFeature <- colnames(df)[grepl(Pro, colnames(df)) & grepl('Hippocampus', colnames(df))]

df_sel <- df[, c('group', selectedFeature)]

data <- df_sel %>% pivot_longer(-group, names_to="Marker", values_to="Expression")
data <- as.data.frame(data)
data$group <- ordered(factor(data$group), levels=c('Control', 'AD', 'LBD', 'ODC'))
data$Marker <- factor(data$Marker, levels=unique(data$Marker), ordered=TRUE)
data <- data[data$group %in% c('Control', 'AD'), ]
data <- data %>% separate(Marker, c('Region', 'Marker', 'Cluster'), sep=' ') 
data$group <- droplevels(data$group)



#### Frequency plots ----------------------------------------------------------
# Importing data
region <- 'Hipp'
df <- read.csv(paste0('R_py_exchange/df_freq_', region, '_noStd.csv'))[, -1]
# colnames(df)[3:ncol(df)] <- sapply(colnames(df)[3:ncol(df)], function(x) strsplit(x, 'X')[[1]][2])
freq <- df[, -2] %>% pivot_longer(-group, names_to="Cluster", values_to="Frequency")
freq <- as.data.frame(freq)
freq <- freq[freq$group!='ODC', ]
levels(freq[, 'group']) <- c("LBD", "Control", "ODC", "AD" )
freq$group <- ordered(factor(freq$group), levels=c('Control', 'AD', 'LBD', 'ODC'))
freq <- freq[freq$group %in% c('Control', 'AD'), ]
freq <- freq %>% group_by(group, Cluster) %>% summarise(mean=mean(Frequency))
freq$group <- droplevels(freq$group)
freq[, 'A'] <- grepl('A', freq$Cluster)
sumFreq <- freq %>% group_by(group, A) %>% summarise(mean=sum(mean))
freq[freq$A & freq$group=='Control', 'mean'] <- freq[freq$A & freq$group=='Control', 'mean']/sumFreq[sumFreq$A & sumFreq$group=='Control', 'mean'][[1]]
freq[freq$A & freq$group=='AD', 'mean'] <- freq[freq$A & freq$group=='AD', 'mean']/sumFreq[sumFreq$A & sumFreq$group=='AD', 'mean'][[1]]


combined <- data %>% inner_join(freq, by=c('group', 'Cluster'))
combined[, 'wtedMean'] <- combined[, 'Expression'] * combined[, 'mean']
combined[, 'A'] <- grepl('A', combined$Cluster)
combined <- combined %>% mutate(A=ordered(factor(ifelse(A, 'A', 'Non-A')), levels=c('Non-A', 'A')))


combined %>% group_by(group, A) %>% summarise(mean=mean(wtedMean), sd=sd(wtedMean))

if(Pro=='PHF-tau') {max_y <- 0.17} else {max_y <- 0.75} # cd47
p <- ggbarplot(as.data.frame(combined), x='A', y='wtedMean', line.color="gray", line.size=0.4, point.size=3, facet.by='group',
               color='black', palette="jco", fill="A", add=c("mean_se"), ylim=c(0, max_y)) + 
               stat_compare_means(method='wilcox', label.y=max_y-0.1*max_y, size=6,
                                  mapping=aes(label=paste('Wilcoxon, p = ', format.pval(..p.adj.., digits=3)))) +
               xlab('') + ylab('SynTOF Weighted vGLUT Expression') + theme_classic(base_size=23) + theme(legend.position='none')
ggsave('figures/barplots/mibi_syntof_weighted_vGLUT_expression.pdf', p)


