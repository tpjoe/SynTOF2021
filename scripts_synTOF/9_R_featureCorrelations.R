"""
This script plots out/calculate correlations between two markers
"""
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
# df <- df[df$group=='LowNo', ]
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



###########################################
# correlation plot between EAAT and GFAP #
###########################################

# calculations
eaat <- df[, (grepl('EAAT', colnames(df)))]
gfap <- df[,  (grepl('GFAP', colnames(df)))]
p <- sapply(1:ncol(eaat), function(x) cor.test(eaat[, x], gfap[, x], method='spearman')$p.value)
rho <- sapply(1:ncol(eaat), function(x) cor.test(eaat[, x], gfap[, x], method='spearman')$estimate)

cbind.data.frame(cl=colnames(eaat), p=p, rho=rho)

# start plotting

plot_list <- list()
for (i in 1:ncol(eaat)) {
    data <- cbind.data.frame(Diagnosis=df[, 'group'], EAAT1=(eaat[, i]), GFAP=(gfap[, i]))
    title <- gsub(' EAAT1,', '', colnames(eaat)[i])

    plot_list[[i]] <- ggscatter(data, x="GFAP", y="EAAT1", fill='Diagnosis', size=4.5, color=transparent('white', trans.val=1), 
        add="reg.line", add.params=list(color='#4c4b4d', fill="lightgray"), # Customize reg. line
        conf.int=TRUE, shape=21,
        cor.coef=TRUE, # Add correlation coefficient. see ?stat_cor
        cor.coeff.args = list(method="spearman", label.sep="\n", size=5),
        palette=c("#2a18f2", "black", '#fac116', '#fa164f')) + 
        theme_classic() + theme(text=element_text(size=14), axis.text=element_text(size=14), legend.position="none",
        plot.title=element_text(hjust=0.5)) + 
        ylab('EAAT1 Value') + xlab('GFAP Value') + ggtitle(title)
}

plot <- ggarrange(plotlist=plot_list, ncol=5, nrow=9)#, widths =c(2,1), labels = c("a", "b"))

ggsave(filename=paste0('figures/scatter/EAAT1vsGFAP.pdf'), plot, width=15, height=25)






###########################################
# correlation cal. between PHF-tau and CD47 #
###########################################

# for lowno and combined
df_lowno <- df[df$group=='LowNo', ]
region <- 'Hippo'

# get values
phftau_hipp <- df_lowno[, (grepl('PHF-tau', colnames(df_lowno))) & (grepl(paste0(region, ', '), colnames(df_lowno)))]
cd47_hipp <- df_lowno[,  (grepl('CD47', colnames(df_lowno))) & (grepl(paste0(region, ', '), colnames(df_lowno)))]
phftau_hipp <- as.vector(as.matrix(phftau_hipp))
cd47_hipp <- as.vector(as.matrix(cd47_hipp))

# calculate cor
cor.test(phftau_hipp, cd47_hipp, method='spearman')





###########################################
# correlation cal. between PHF-tau and CD47 #
###########################################

# for lowNo only
df_lowno <- df[df$group=='LowNo', ]
region <- 'BA9'
cd47 <- df_lowno[, (grepl('CD47', colnames(df_lowno))) & (grepl(paste0(region, ', '), colnames(df_lowno)))]
snap25 <- df_lowno[,  (grepl('SNAP25', colnames(df_lowno))) & (grepl(paste0(region, ', '), colnames(df_lowno)))]
cd47 <- as.vector(as.matrix(cd47))
snap25 <- as.vector(as.matrix(snap25))

cor.test(cd47, snap25, method='spearman')

# for all groups
cd47 <- df[, (grepl('CD47', colnames(df)))] <- df[, (grepl('CD47', colnames(df)))]
snap25 <- df[,  (grepl('SNAP25', colnames(df)))]
p <- sapply(1:ncol(cd47), function(x) cor.test(cd47[, x], snap25[, x], method='spearman')$p.value)
rho <- sapply(1:ncol(cd47), function(x) cor.test(cd47[, x], snap25[, x], method='spearman')$estimate)


