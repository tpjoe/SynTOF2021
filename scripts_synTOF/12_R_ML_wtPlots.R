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

pair <- c('LowNo', 'LBD') #<----- This is where you choose between ['LowNo', LBD'] or ['LowNo', 'PHAD']

# fill na and get rid of all-0 columns
for(i in 1:ncol(df)){
  df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
}

# get rid of other disease controls
df <- df[df$group!='ODC', ]

# rename a bit
colnames(df) <- gsub('_X3NT_', '_3NT_', colnames(df))

# filter the dataframe to only those we want
X <- df[df$group %in% pair, !colnames(df) %in% c('group', 'sample')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]
y <- df[df$group %in% pair, 'group']
y <- sapply(y, function(x) if (x==pair[1]) {0} else {1})


# load weights
wts <- read.csv(paste0('R_py_exchange_afterCluster/wt_LowNo_', pair[2], '_Ridge.csv'))[, -1] #<-------- for LBD
# wts <- read.csv(paste0('R_py_exchange_afterCluster/wt_LowNo_', pair[2], '_EN.csv'))[, -1] #<-------- for AD
mean_median <- 'mean'

wts_mean <- apply(wts, 2, mean_median)
wts_mean <- cbind.data.frame(names(wts_mean), wts_mean)
rownames(wts_mean) <- NULL
wts_mean[, 1] <- make.names(wts_mean[, 1])


# add rownames that's not in there (those perfect features)
pTauName <- colnames(X)[!colnames(X) %in% wts_mean[, 1]]
pTau <- cbind.data.frame(pTauName, rep(mean((wts_mean[, 2])), length(pTauName)))
colnames(pTau) <- colnames(wts_mean)
wts_mean <- rbind.data.frame(wts_mean, pTau)
wts_mean <- wts_mean[match(colnames(X), wts_mean[, 1]), ]


# calculate q values
corrP_fdr <- get_fdrp(apply(X, 2, function(x) cor(x, y, method='spearman')), log=FALSE, n=nrow(X))
logQ <- -log10(corrP_fdr)

# color by marker type
if (pair[2]=='PHAD') {rand<-10301} else {rand<-10301} #99
set.seed(rand) 
marker_names <- sapply(names(logQ), function(x) paste(strsplit(x, '_')[[1]][2], strsplit(x, '_')[[1]][3], sep='_'))
marker_names_numeric <- as.numeric(factor(marker_names))
names(marker_names_numeric) <- names(logQ)
cl_palette <- distinctColorPalette(k=length(unique(marker_names_numeric)))
cl_palette_ <- sapply(marker_names_numeric, function(x) cl_palette[x])
# cl_palette_[corrP_fdr > 0.05] <- "#e3e3e1" #transparent("grey", trans.val=0.85)
size_ <- rep(1.5, length(cl_palette_))
size_[corrP_fdr > 0.05] <- 1


# making df plot
df_plot <- cbind.data.frame(wt=wts_mean[, 2], logq=logQ, col=cl_palette_, names=names(cl_palette_))
df_plot <- df_plot[order(abs(df_plot$logq), decreasing=TRUE), ]
rownames(df_plot) <- NULL
# if (pair[2] == 'PHAD') {
    df_plot <- df_plot[!df_plot[, 4] %in% pTauName, ]  #################### remove high correlaiton features for PHAD ################################
# }

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
df_plot$names <- gsub('Hipp', 'Hippocampus', df_plot$names)
df_plot$names <- gsub('DLCau', 'Caudate', df_plot$names)

# set number of label to be shown
if (pair[2]=='PHAD') {n_show<-11} else {n_show<-2} #LBD = 2, PHAD = 14
sel_names_ <- as.character(df_plot[order(abs(df_plot$logq), abs(df_plot$wt), decreasing=TRUE), 'names'][1:1000])
sel_names <- sel_names_[!duplicated(sapply(sel_names_, function(x) strsplit(x, ',')[[1]][2]))][1:n_show]


# define shape for different brain regions
df_plot['shape'] <- factor(sapply(df_plot$names, function(x) strsplit(x, ',')[[1]][1]), labels=c('BA9', 'Caudate', 'Hippocampus'))


# for main manuscript plot
p <- ggscatter(df_plot, x="wt", y="logq", color="col", size=2.5, label="names", label.select=sel_names,
               repel=TRUE, font.label=list(size=12), alpha=sapply(df_plot[, 'logq'], function(x) if(x>(-log10(0.05))) {0.8} else {0.2}), shape='shape') + 
     theme_classic() +
     scale_color_manual(breaks=unique(df_plot$col), values=as.character(unique(df_plot$col))) + 
     theme(legend.position="top", legend.title = element_blank(), legend.text=element_text(size=16), 
          axis.text=element_text(size=16), axis.title=element_text(size=16), plot.margin=margin(2, 14, 2, 2)) +
          xlab('Mean Ridge Coefficients') + ylab('log Q') +
     geom_hline(yintercept=-log10(0.05), linetype="dashed") + guides(color="none")

ggsave(filename=paste0('figures/coefficients/LBD_test.pdf'), p, width=5, height=5)#, width=20, height=15)