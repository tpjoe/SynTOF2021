"""
This script plots out the table of Q values in the supplementals for 9 AD vs. 6 Control and 
7 non-resilient AD vs. Control.
"""

#### libraries -----------------------------------------------------------------------
source('R_utils_plots.R')
library(Rtsne)
library(yarrr)
library(uwot)
library(randomcoloR)
library(plotrix)
library(ggpubr)
library(ltm)
library(yarrr)
library(DescTools)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)

get_fdrp <- function(r, n=24, log=TRUE) {
    # zY <- sqrt((n-2)/1.06)*atanh(r)
    # aaP <- 2*pnorm(-abs(zY))
    zY <- r*sqrt((n-2)/((r+1)*(1-r)))
    aaP <- 2*(1-pt(abs(zY), df=(n-2)))
    aaP <- p.adjust(aaP, method="fdr")
    if (log==TRUE) {
        aaP <- -log10(aaP)
    }
    aaP
}

#######################################################################
# For AD (non-resilient) vs Control #
#######################################################################

# Importing data 
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
# fill na and get rid of all-0 columns
for(i in 1:ncol(df)){
  df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
}
df <- df[df$group!='ODC', ]

# define X and y from diagnosis
colnames(df) <- gsub('_X3NT_', '_3NT_', colnames(df))
X <- df[df$group %in% pair, !colnames(df) %in% c('group', 'sample')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]
y <- df[df$group %in% pair, 'group']
y <- sapply(y, function(x) if (x==pair[1]) {0} else {1})

# ----------------------------------------------------------------------------------------------------------------------

# load weights
wts <- read.csv(paste0('R_py_exchange_afterCluster/wt_LowNo_', pair[2], '_EN.csv'))[, -1]
mean_median <- 'mean'

wts_mean <- apply(wts, 2, mean_median)
wts_mean <- cbind.data.frame(names(wts_mean), wts_mean)
rownames(wts_mean) <- NULL
wts_mean[, 1] <- make.names(wts_mean[, 1])
wts_mean <- wts_mean[match(colnames(X), wts_mean[, 1]), ]

# calculate Q values
corrP_fdr <- get_fdrp(apply(X, 2, function(x) cor(x, y, method='spearman')), log=FALSE, n=nrow(X))
logQ <- -log10(corrP_fdr)

# color by marker type
if (pair[2]=='PHAD') {rand<-5} else {rand<-99}
set.seed(rand) 
marker_names <- sapply(names(logQ), function(x) paste(strsplit(x, '_')[[1]][2], strsplit(x, '_')[[1]][3], sep='_'))
marker_names_numeric <- as.numeric(factor(marker_names))
names(marker_names_numeric) <- names(logQ)
cl_palette <- distinctColorPalette(k=length(unique(marker_names_numeric)))
cl_palette_ <- sapply(marker_names_numeric, function(x) cl_palette[x])
size_ <- rep(1.5, length(cl_palette_))
size_[corrP_fdr > 0.05] <- 1

# making df plot
df_plot <- cbind.data.frame(wt=wts_mean[, 2], logq=logQ, col=cl_palette_, names=names(cl_palette_))
df_plot <- df_plot[order(abs(df_plot$logq), decreasing=TRUE), ]
rownames(df_plot) <- NULL

# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
          'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', '3NT', 'p\\.Tau')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', '3-NT', 'PHF-tau')
for (s in 1:length(to_sub)) {
    df_plot$names <- gsub(to_sub[s], sub_to[s], df_plot$names)
}
df_plot$names <- gsub('_mean_', '_', df_plot$names)
df_plot$names <- gsub('_', ',', df_plot$names)
df_plot$names <- gsub('Hipp', 'Hippo', df_plot$names)
df_plot$names <- gsub('DLCau', 'CN', df_plot$names)


sel_names_ <- as.character(df_plot[order(abs(df_plot$logq), abs(df_plot$wt), decreasing=TRUE), 'names'][1:1000])
sel_names <- sel_names_[!duplicated(sapply(sel_names_, function(x) strsplit(x, ',')[[1]][2]))][1:5]
df_plot['shape'] <- factor(sapply(df_plot$names, function(x) strsplit(x, ',')[[1]][1]), labels=c('BA9', 'CN', 'Hippo'))


# keep only those with Q <0.05 and define marker names
df_plot2 <- df_plot[df_plot$logq > -log10(0.05), ]
df_plot2['marker'] <- sapply(df_plot2$names, function(x) strsplit(x, ',')[[1]][2])
df_plot2['sub'] <- sapply(df_plot2$names, function(x) paste(strsplit(x, ',')[[1]][1], strsplit(x, ',')[[1]][3], sep=','))


df_plot2_ <- df_plot2[order(df_plot2$marker, decreasing=TRUE), ]
allAD <- df_plot2_

# start plotting 
p2 <- ggscatter(df_plot2_, x="marker", y="logq", color='marker', size=5, palette=rep(c("#56147d", "#d6446d", "#fea16e"), 8))+ 
        geom_text_repel(data=df_plot2, aes(x = marker, y = logq, label=sub, color=marker), force=3, size=6) +
        theme(legend.position="none", legend.title=element_blank(), legend.text=element_blank(), 
          axis.text=element_text(size=16), axis.title=element_text(size=16), plot.margin=margin(2, 14, 2, 2)) +
      xlab('') + ylab('log Q (Control vs. AD)') + coord_flip()+
      geom_hline(yintercept=c(-log10(0.01)), linetype="dotted")

ggsave(filename=paste0('figures/Q_values/lowno_AD.pdf'), p2, width=7, height=12)#, width=20, height=15)




#######################################################################
# For AD (non-resilient) vs Control #
#######################################################################

# importing data
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

# fill na and get rid of all-0 columns
for(i in 1:ncol(df)){
  df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
}
gt <- read.csv('../raw_data/demographics.csv')
gt['sample_ID'] <- paste0(gsub(' ', '', gt[['sample_ID']]), '.fcs')
gt <- gt[gt[['sample_ID']] %in% df$sample, ]
gt <- gt[match(df$sample, gt[['sample_ID']]), ]
all(gt[['sample_ID']] == df$sample)
colnames(gt)
colnames(df) <- gsub('_X3NT_', '_3NT_', colnames(df))

# define X and Y from resilience status
gt <- gt[(!df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') | df$group=='LowNo'), ]
df <- df[(!df$sample %in% c('HF14-023.fcs', 'HF14-081.fcs', 'HF13-030.fcs', 'HF13-095.fcs') | df$group=='LowNo'), ]
df['resilient'] <- df$group
levels(df$resilient)[levels(df$resilient)=="AD"] <- "resilient"
levels(df$resilient)[levels(df$resilient)=="LBD"] <- "resilient"
df$resilient <- factor(df$resilient, levels = c("LowNo", "resilient", "ODC"))
X <- df[, !colnames(df) %in% c('group', 'sample', 'resilient')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]
y <- gt[['CognitiveStatus']]
y <- sapply(y, function(x) if (x=='No dementia') {0} else {1})
X <- X[df$group %in% c('LowNo', 'PHAD'), ]
y <- y[df$group %in% c('LowNo', 'PHAD')]

# calculate Q values
corr_ <- apply(X, 2, function(x) cor(x, y, method='spearman'))
corrP_fdr <- get_fdrp(corr_, n=nrow(X), log=FALSE)
logQ <- get_fdrp(corr_, n=nrow(X), log=TRUE)
names(logQ) <- colnames(X)

# color by marker type
set.seed(999)
marker_names <- sapply(names(logQ), function(x) paste(strsplit(x, '_')[[1]][2], strsplit(x, '_')[[1]][3], sep='_'))
marker_names_numeric <- as.numeric(factor(marker_names))
names(marker_names_numeric) <- names(logQ)
cl_palette <- distinctColorPalette(k=length(unique(marker_names_numeric)))
cl_palette_ <- sapply(marker_names_numeric, function(x) cl_palette[x])
cl_palette_[corrP_fdr > 0.05] <- "#e3e3e1" #transparent("grey", trans.val=0.85)
size_ <- rep(1.5, length(cl_palette_))
size_[corrP_fdr > 0.05] <- 1

# making df plot
df_plot <- cbind.data.frame(logq=logQ, col=cl_palette_, names=names(cl_palette_))
df_plot <- df_plot[order(abs(df_plot$logq), decreasing=TRUE), ]
rownames(df_plot) <- NULL

# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
          'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', '3NT', 'p\\.Tau')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', '3-NT', 'PHF-tau')
for (s in 1:length(to_sub)) {
    df_plot$names <- gsub(to_sub[s], sub_to[s], df_plot$names)
}
df_plot$names <- gsub('_mean_', '_', df_plot$names)
df_plot$names <- gsub('_', ',', df_plot$names)
df_plot$names <- gsub('Hipp', 'Hippo', df_plot$names)
df_plot$names <- gsub('DLCau', 'CN', df_plot$names)

# define shape
df_plot['shape'] <- factor(sapply(df_plot$names, function(x) strsplit(x, ',')[[1]][1]), labels=c('BA9', 'CN', 'Hippo'))

# keep only those with Q <0.05 and define marker names
df_plot2 <- df_plot[df_plot$logq > -log10(0.05), ]
df_plot2['marker'] <- sapply(df_plot2$names, function(x) strsplit(x, ',')[[1]][2])
df_plot2['sub'] <- sapply(df_plot2$names, function(x) paste(strsplit(x, ',')[[1]][1], strsplit(x, ',')[[1]][3], sep=','))

# reorder markers
df_plot2 <- df_plot2[order(df_plot2$marker, decreasing=TRUE), ]

# start plotting 
p2 <- ggscatter(df_plot2, x="marker", y="logq", color='marker', size=5, palette=rep(c("#56147d", "#d6446d", "#fea16e"), 8))+ 
        geom_text_repel(data=df_plot2, aes(x = marker, y = logq, label=sub, color=marker), force=3, size=6) +
        theme(legend.position="none", legend.title=element_blank(), legend.text=element_blank(), 
          axis.text=element_text(size=16), axis.title=element_text(size=16), plot.margin=margin(2, 14, 2, 2)) +
      xlab('') + ylab('log Q (Control vs. AD Dementia)') + coord_flip()  +
      geom_hline(yintercept=c(-log10(0.01)), linetype="dotted")

ggsave(filename=paste0('figures/Q_values/lowno_non-resAD.pdf'), p2, width=14, height=12)#, width=20, height=15)