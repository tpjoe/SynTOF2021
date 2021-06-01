"""
This script generates the candybar plot that illustrates the changes in Q values between
1. 9 ADNC cases vs. 6 Control
2. 7 AD-non resilient vs. 6 Control
3. 7 Mixed cases (with 2 resilient cases included) vs. 6 Control
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

##########################################################################
# start with [5 non-resilient AD + 2 resilient] vs 6 Controls #
##########################################################################

# get data
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

# get ground truth (resilience)
gt <- read.csv('../raw_data/demographics.csv')
gt['sample_ID'] <- paste0(gsub(' ', '', gt[['sample_ID']]), '.fcs')
gt <- gt[gt[['sample_ID']] %in% df$sample, ]
gt <- gt[match(df$sample, gt[['sample_ID']]), ]
all(gt[['sample_ID']] == df$sample)
colnames(df) <- gsub('_X3NT_', '_3NT_', colnames(df))

# define resilient cases 
res_c <- c('HF13-030.fcs', 'HF13-095.fcs')
removal_pair <- t(combn(df$sample[!df$sample %in% res_c & df$group=='PHAD'], 2))

# store backup dataframes
df_backup <- df
gt_backup <- gt

# define blank q value results
qValue_table <- matrix(nrow=length(colnames(df)[3:ncol(df)]), ncol=nrow(removal_pair), 
                      dimnames=list(colnames(df)[3:ncol(df)], seq(nrow(removal_pair))))

# run all possible combinations of 2 non-resilient AD removals and calculate Q
for (n_pair in 1:nrow(removal_pair)) {
    pair <- removal_pair[n_pair, ]
    gt <- gt_backup
    df <- df_backup

    gt <- gt[!(df$sample %in% pair) & df$group %in% c('PHAD', 'LowNo'), ]
    df <- df[!(df$sample %in% pair) & df$group %in% c('PHAD', 'LowNo'), ]

    df['resilient'] <- df$group
    levels(df$resilient)[levels(df$resilient)=="AD"] <- "resilient"
    levels(df$resilient)[levels(df$resilient)=="LBD"] <- "resilient"
    df$resilient <- factor(df$resilient, levels = c("LowNo", "resilient", "ODC"))

    X <- df[, !colnames(df) %in% c('group', 'sample', 'resilient')]
    X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]

    y <- df$group
    y <- sapply(y, function(x) if (x=='LowNo') {0} else {1})
    # start plotting -------------------------------------------------------------------------------------------------------
    corr_ <- apply(X, 2, function(x) cor(x, y, method='spearman'))
    corrP_fdr <- get_fdrp(corr_, n=nrow(X), log=FALSE)
    qValue_table[, n_pair] <- get_fdrp(corr_, n=nrow(X), log=TRUE)
}

# only keep those with Q < 0.05
qValue_table[qValue_table < -log10(0.05)] <- NA 

# rename columns a bit
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
        'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', '3NT', 'p\\.Tau')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', '3-NT', 'PHF-tau')
for (s in 1:length(to_sub)) {
    rownames(qValue_table) <- gsub(to_sub[s], sub_to[s], rownames(qValue_table))
}
rownames(qValue_table) <- gsub('_mean_', '_', rownames(qValue_table))
rownames(qValue_table) <- gsub('_', ',', rownames(qValue_table))
rownames(qValue_table) <- gsub('Hipp', 'Hippo', rownames(qValue_table))
rownames(qValue_table) <- gsub('DLCau', 'CN', rownames(qValue_table))

# calculate mean values and median for dot sizes 
qValue_table_others <- apply(qValue_table, 1, function(x) median(x, na.rm=TRUE)) %>% as.data.frame %>% set_colnames('size')
qValue_mean <- apply(qValue_table, 1, function(x) mean(!is.na(x), na.rm=TRUE))
qValue_sd <- apply(qValue_table, 1, function(x) sd(!is.na(x), na.rm=TRUE))

# aggregate by brain region and values
qValue_all <- cbind.data.frame(
                               mean=qValue_mean, 
                               sd=qValue_sd,
                               region=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][1]),
                               cl=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][3]),
                               pro=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][2]),
                               size=10^-qValue_table_others
                               )
qValue_mean_reserved <- qValue_all

# calculate count of mean values that are not NAs (i.e. <0.05) and standard deviations
qValue_all %<>% group_by(region, pro) %>% summarise(mean_count=sum(mean), sd_count=sqrt(sum(sd^2)/15), size=median(size, na.rm=T)) #15=number of clusters




##########################################################################
# start with 7 non-resilient AD vs Control #
##########################################################################

gt <- gt_backup
df <- df_backup

# only one possible pair to remove, i.e. resilient AD
n_pair <- 1
removal_pair <- t(combn(df$sample[df$sample %in% res_c & df$group=='PHAD'], 2))

# define blank Q value result matrix
qValue_table <- matrix(nrow=length(colnames(df)[3:ncol(df)]), ncol=nrow(removal_pair), 
                      dimnames=list(colnames(df)[3:ncol(df)], seq(nrow(removal_pair))))
pair <- removal_pair[n_pair, ]

# define X and y essentially
gt <- gt[!(df$sample %in% pair) & df$group %in% c('PHAD', 'LowNo'), ]
df <- df[!(df$sample %in% pair) & df$group %in% c('PHAD', 'LowNo'), ]
df['resilient'] <- df$group
levels(df$resilient)[levels(df$resilient)=="AD"] <- "resilient"
levels(df$resilient)[levels(df$resilient)=="LBD"] <- "resilient"
df$resilient <- factor(df$resilient, levels = c("LowNo", "resilient", "ODC"))
X <- df[, !colnames(df) %in% c('group', 'sample', 'resilient')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]
y <- df$group
y <- sapply(y, function(x) if (x=='LowNo') {0} else {1})

# calculate Q values
corr_ <- apply(X, 2, function(x) cor(x, y, method='spearman'))
corrP_fdr <- get_fdrp(corr_, n=nrow(X), log=FALSE)
qValue_table[, n_pair] <- get_fdrp(corr_, n=nrow(X), log=TRUE)
qValue_table[qValue_table < -log10(0.05)] <- NA 

# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
        'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', '3NT', 'p\\.Tau')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', '3-NT', 'PHF-tau')
for (s in 1:length(to_sub)) {
    rownames(qValue_table) <- gsub(to_sub[s], sub_to[s], rownames(qValue_table))
}
rownames(qValue_table) <- gsub('_mean_', '_', rownames(qValue_table))
rownames(qValue_table) <- gsub('_', ',', rownames(qValue_table))
rownames(qValue_table) <- gsub('Hipp', 'Hippo', rownames(qValue_table))
rownames(qValue_table) <- gsub('DLCau', 'CN', rownames(qValue_table))

# calculate mean and dot sizes
qValue_table_res <- qValue_table %>% set_colnames(c("size"))
qValue_mean <- apply(qValue_table, 1, function(x) mean(!is.na(x), na.rm=TRUE))
qValue_mean <- cbind.data.frame(
                               mean=qValue_mean, 
                               region=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][1]),
                               cl=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][3]),
                               pro=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][2]),
                               size=10^-qValue_table_res
                               )

# count mean that is not NA
qValue_all <- qValue_mean %>% group_by(region, pro) %>% summarise(res_count=sum(mean), res_size=median(size, na.rm=TRUE)) %>% merge(qValue_all)




##########################################################################
# start with 7 non-resilient AD vs Control #
##########################################################################
# calculate q value of all cases
gt <- gt_backup
df <- df_backup

qValue_table <- matrix(nrow=length(colnames(df)[3:ncol(df)]), ncol=1, 
                      dimnames=list(colnames(df)[3:ncol(df)], seq(1)))

gt <- gt[df$group %in% c('PHAD', 'LowNo'), ]
df <- df[df$group %in% c('PHAD', 'LowNo'), ]

df['resilient'] <- df$group
levels(df$resilient)[levels(df$resilient)=="AD"] <- "resilient"
levels(df$resilient)[levels(df$resilient)=="LBD"] <- "resilient"
df$resilient <- factor(df$resilient, levels = c("LowNo", "resilient", "ODC"))

X <- df[, !colnames(df) %in% c('group', 'sample', 'resilient')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]

y <- df$group
y <- sapply(y, function(x) if (x=='LowNo') {0} else {1})

corr_ <- apply(X, 2, function(x) cor(x, y, method='spearman'))
corrP_fdr <- get_fdrp(corr_, n=nrow(X), log=FALSE)
qValue_table[, n_pair] <- get_fdrp(corr_, n=nrow(X), log=TRUE)

qValue_table[qValue_table < -log10(0.05)] <- NA 

# rename
to_sub <- c('DJ.1_PARK7', 'TMEM230_C20orf30', 'PrP_CD230', 'GATM_1', 'GAMT_2', 'PARKIN', 'a\\.Synuclein_pS129', 
        'b.Amyloid_X40', 'Casp3_Acti', 'b.Amyloid_X42', 'K48.Ubiquitin', 'a.Synuclein', '3NT', 'p\\.Tau')
sub_to <- c('DJ1', 'TEMEM230', 'PrP', 'GATM', 'GAMT', 'Parkin', 'p129-AS', 
            'Ab40', 'Ac-Casp3', 'Ab42', 'K48', 'AS', '3-NT', 'PHF-tau')
for (s in 1:length(to_sub)) {
    rownames(qValue_table) <- gsub(to_sub[s], sub_to[s], rownames(qValue_table))
}
rownames(qValue_table) <- gsub('_mean_', '_', rownames(qValue_table))
rownames(qValue_table) <- gsub('_', ',', rownames(qValue_table))
rownames(qValue_table) <- gsub('Hipp', 'Hippo', rownames(qValue_table))
rownames(qValue_table) <- gsub('DLCau', 'CN', rownames(qValue_table))

qValue_table_res <- qValue_table %>% set_colnames(c("size"))
qValue_mean <- apply(qValue_table, 1, function(x) mean(!is.na(x), na.rm=TRUE))
qValue_mean <- cbind.data.frame(
                               mean=qValue_mean, 
                               region=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][1]),
                               cl=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][3]),
                               pro=sapply(names(qValue_mean), function(x) strsplit(x, ',')[[1]][2]),
                               size=10^-qValue_table_res
                               )



# calculate mean q values from different combinations of 2 non- resilient AD removal
qValue_all <- qValue_mean %>% group_by(region, pro) %>% summarise(all_count=sum(mean), all_size=median(size, na.rm=TRUE)) %>% merge(qValue_all)
qValue_all <- qValue_all[!(qValue_all$all_count<0.5 & qValue_all$mean_count<0.5 & qValue_all$res_count<0.5), ]
qValue_all %<>% replace_na(list(size=0, res_size=0, all_size=0))
qValue_all['pro'] <- factor(qValue_all$pro, levels=rev(levels(qValue_all$pro)))
qValue_all['region'] <- factor(qValue_all$region, levels=(levels(qValue_all$region)))



##########################################################################
# Plotting #
##########################################################################
dodge <- 0.5
ggplot(qValue_all) +
geom_linerange(aes(y=res_count, ymin=mean_count, x=pro, ymax=res_count, group=desc(region) ), size=1.5,
                   color="#DCDCDC", position=position_dodge(width=dodge)) +
geom_point( aes(x=pro, y=all_count, shape=region, size=1/(all_size), group=desc(region) ), color=rgb(214/255, 68/255, 109/255), 
            fill=rgb(214/255, 68/255, 109/255, 0.75), position=position_dodge(width=dodge)) + # AD no res
geom_point( aes(x=pro, y=mean_count, shape=region, size=1/((size)), group=desc(region) ), color=rgb(254/255, 191/255, 132/255), 
            fill=rgb(254/255, 191/255, 132/255, 0.75), position=position_dodge(width=dodge)) + #AD dementia
geom_errorbar(aes(x=pro, ymin=mean_count-sd_count, ymax=mean_count+sd_count, width=0, group=desc(region) ), 
            position=position_dodge(width=dodge), size=0.5,  color='black') +
geom_point( aes(x=pro, y=res_count, shape=region, size=1/((res_size)), group=desc(region) ), color=rgb(9/255, 7/255, 34/255), 
            fill=rgb(9/255, 7/255, 34/255, 0.75), position=position_dodge(width=dodge)) + #AD All
scale_shape_manual(values=c(21, 22, 24)) +
coord_flip() +
theme_classic() +
theme(legend.position="bottom", panel.grid.major.x=element_line(colour="gray", size=0.3, linetype='dashed'),
      legend.title = element_text(size = 14), legend.text = element_text(size = 14),
      axis.text=element_text(size=14), axis.title=element_text(size=14), plot.margin=unit(c(0,1,0,0), "cm")) +
xlab("") +
ylab("Number of Subpopulations with Q Value < 0.05") +
scale_y_continuous(breaks=seq(0, 16, 1)) + 
guides(size=FALSE, shape=guide_legend(title="Region", override.aes=list(color='black', fill='black', size=c(7.5, 5, 3.75))))
# guides(fill=FALSE) # + guides(shape=)

ggsave(filename=paste0('figures/Q_values/candybar_compare.pdf'), height=11, width=6.5)