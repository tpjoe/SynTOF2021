
"""
This script generates heatmap of the mean expression values from the LowNo group
"""
# Expression visualization and heatmaps -------------------------------------------------------------------------
# load data 
library(dplyr)
library(tidyr)
library(ggplot2)
library(flowCore)
library(data.table)
library(reshape2)
library(gplots)

# load cluster assignment
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
mc[mc$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), 'group'] <- 'ODC'
df_region <- readRDS(paste0('df_pre_save/exp', 2, 'mc_5_13.rds'))

# select group that you want
df_ <- df_region
group <- 'LowNo'
region <- '' #'' = all regions, otherwise specify among c('BA9', 'Hipp', 'DLCau')

# define marker types
functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
surfacePro <- sort(colnames(df_[['BA9']][[group]])[!colnames(df_[['BA9']][[group]]) %in% 
                   c('mc', 'sample', functionalPro)], decreasing=TRUE)

# start screening data
df__ <- data.frame()
if (region!='') {
    mc_sub <- mc[(mc$pp=='pre') & (mc$group==group) & (mc$region==region), 'mc'] 
    df__ <- df_[[region]][[group]][, surfacePro]
} else {
    mc_sub <- mc[(mc$pp=='pre') & (mc$group==group), 'mc'] # specify mc that corresponds to the df you loaded
    for (region in names(df_)) {
        df__ <- rbind.data.frame(df__, df_[[region]][[group]][, surfacePro])
    }
}

# rename columns a bit
colnames(df__) <- sapply(colnames(df__), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                   else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                   else if(x=='PrP_CD230') {'PrP'}
                                                   else if(x=='GATM_1') {'GATM'}
                                                   else if(x=='GAMT_2') {'GAMT'}
                                                   else if(x=='PARKIN') {'Parkin'}
                                                   else if (x=='a-Synuclein') {'AS'}
                                                   else if (x=='VGLUT') {'vGLUT'}
                                                   else x) 
mylevels <- colnames(df__)


# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc_sub <- sapply(mc_sub, function(x) new_label[x])
df_clus <- cbind.data.frame(mc_sub, df__)

# aggregate into table (by mean values)
aa <- aggregate(df_clus[, -1], by=list(df_clus$mc_sub), FUN=mean)
aaa <- aa[, sort(colnames(aa[, -1]))]
rownames(aaa) <- aa[, 1]

# specify color palette
redgreen <- c("black", "#20114bff", "#1c1043ff", '#56147dff', '#d6446dff', "#fea16eff", '#fcf7b9ff')
pal <- colorRampPalette(redgreen)(100)

# reorder cluster
aaa <- aaa[c("A1", "A2", "B1", "B2", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11"), ]

pdf(paste0('figures/heatmaps/all', '_', group, '_meanNoStd.pdf'))
# pdf(paste0('figures/heatmaps/all', '_', group, '_', region, '_meanNoStd.pdf'))
heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=as.dendrogram(hclust(dist((as.matrix(aaa))))), col=pal, #viridis_pal(begin=0, end=0.8, option='magma'), #magma(256, begin=0),
# heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=FALSE, col=pal, main="BA9",
          scale='row', trace='none', srtCol=90, adjCol = c(0.75, 0.5),
          lhei=c(2, 12), lwid=c(2, 5), key.par=list(cex=0.622), density.info="none",
          key.title='', keysize=0.1, margins=c(3,7)) #+ ggtitle('Caudate')
dev.off()