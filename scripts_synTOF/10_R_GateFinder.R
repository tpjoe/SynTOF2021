"""
This script plots out the automated gating scheme using GateFinder
"""

# import libraries ----------------------------------------------------
library(GateFinder)
library(data.table)
library(tidyr)
library(foreach)
library(doParallel)
library(viridis)
library(ggplot2)
library(ggcyto)
library(gridExtra)

cl <- makeCluster(50)
registerDoParallel(cl)

# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
df_region <- readRDS(paste0('df_pre_save/exp', 0, '.rds'))


# prop_markers <- c('CD56', 'CD47', 'GAD65', 'SNAP25', 'GATM_1', 'a-Synuclein', 'VGLUT', 'BIN1', 'TMEM230_C20orf30',
                #   'ApoE')

# select region and diagnosis group
region <- 'BA9'
group <- 'LowNo'


# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
colnames(df_region[[region]][[group]]) <- sapply(colnames(df_region[[region]][[group]]), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                   else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                   else if(x=='PrP_CD230') {'PrP'}
                                                   else if(x=='GATM_1') {'GATM'}
                                                   else if(x=='GAMT_2') {'GAMT'}
                                                   else if(x=='PARKIN') {'Parkin'}
                                                   else if (x=='a-Synuclein') {'AS'}
                                                   else x) 


# define surface and functional markers
functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))
surfacePro <- sort(colnames(df_region[[region]][[group]])[!colnames(df_region[[region]][[group]]) %in% functionalPro])


# start GateFinder for each cluster
for (cl in new_label) {
    X <- df_region[[region]][[group]][, surfacePro]
    mc_region <- mc[(mc$region==region) & (mc$pp=='pre') & (mc$group==group), ]
    mc_ <- mc_region

    B <- as.matrix(mc_$mc==which(new_label==cl))
    gf <- GateFinder(X, B, max.iter=2)#, nstart=20, subsample=10000) #apply GF

    # density plot bins
    bin_no <- 128

    # step 1 gating plot
    gt <-as.data.frame(attr(gf, 'gates')[1])
    p1 <- ggplot(X, aes_string(x=(colnames(gt)[1]), y=(colnames(gt)[2]))) + geom_hex(bins=bin_no) + scale_fill_viridis(option='magma') + # + scale_fill_distiller(palette ="YlGnBu")
        geom_hex(data=X[mc_$mc==which(new_label==cl), ], bins=bin_no, fill='#348aa7', alpha=0.5) +
        geom_polygon(data=gt, lwd=7, fill=NA, color='#348aa7') +
        theme_bw() +
        theme(legend.position='none', axis.text=element_text(size=17*6), axis.title=element_text(size=20*6),
        panel.border=element_rect(colour="black", fill=NA, size=1*6),
        panel.grid.minor=element_line(size=2), panel.grid.major=element_line(size=2.5))

    # step 2 gating plot
    gt <-as.data.frame(attr(gf, 'gates')[2])
    p2 <- ggplot(X[!attr(gf, 'pops')[[1]], ], aes_string(x=(colnames(gt)[1]), y=(colnames(gt)[2]))) + geom_hex(bins=bin_no, fill='gray') +
        geom_hex(data=X[(attr(gf, 'pops')[[1]]) | (attr(gf, 'pops')[[2]]), ], bins=bin_no)  + scale_fill_viridis(option='magma') + # + scale_fill_distiller(palette ="YlGnBu")
        geom_hex(data=X[mc_$mc==which(new_label==cl), ], bins=bin_no, fill='#348aa7', alpha=0.5) +

        geom_polygon(data=gt, lwd=7, fill=NA, color='#348aa7') +
        theme_bw() +
        theme(legend.position='none', axis.text=element_text(size=17*6), axis.title=element_text(size=20*6),
        panel.border=element_rect(colour="black", fill=NA, size=1*6),
        panel.grid.minor=element_line(size=2), panel.grid.major=element_line(size=2.5))
    
    # save plots together into 1 gifure
    png(paste0('figures/gate/exp00', region, '_', cl,'_', group, '.png'), width=4000, height=2000, pointsize=50)
    grid.arrange(p1, p2, ncol=2, nrow=1)
    dev.off()
}













library(flowCore)
library(GateFinder)
library(ggpubr)
library(dplyr)
library(tidyr)
source('R_utils_preprocessing.R')

# load data -----------------------------------------------------------------------------------
region <- 'Hipp'

# load data
df <- list()
for (group in c('LowNo', 'LBD', 'PHAD')) {
    # get files
    fcs_path <- '../raw_data/max_events/fcs/'
    file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
    df[[group]] <- data.frame()
    # sample <- c()
    # cytoftrans <- arcsinhTransform(transformationId='cytofTransform', a=0, b=(1/5), c=0)
    for (i in seq(length(file_list))) {
        frame <- read.FCS(paste0(fcs_path, file_list[i]))
        ####for applying transformation
        # translist <- transformList(colnames(exprs(frame)), cytoftrans)
        # frame <- transform(frame, translist)
        df[[group]] <- rbind.data.frame(df[[group]], exprs(frame))
        file_name_splt <- strsplit(file_list[i], '_')[[1]]
        # sample <- c(sample, rep(file_name_splt[length(file_name_splt)], nrow(exprs(frame))))
    }
    cutoff <- read.csv('../raw_data/cutoff_gate.csv')
    for (col in colnames(df[[group]])){
        min_value <- min(df[[group]][col])
        cutoff_ <- cutoff[cutoff[, 'SynTOF.Antibody.Panel']==col, 'arcsinh0.2.cutoff']
        df[[group]][df[[group]][col] < cutoff_, col] <- min_value
    }
    df[[group]]['NET'] <- NULL
}




# scaling data by LowNo only
# df_scaled <- list()
# for (group in names(df)) {
#     if (group == 'LowNo') {
#         df_scaled[[group]] <- scale(df[[group]])
#     } else {
#         df_scaled[[group]] <- scale(df[[group]], attr(df_scaled[['LowNo']], "scaled:center"), 
#                                     attr(df_scaled[['LowNo']], "scaled:scale"))
#     }
# }
df_scaled <- df

targetpop <- as.matrix(df_scaled[[3]][1:5000, 'GFAP'] > 1)
##Subset the markers that should be considered for gating.
prop_markers <- c(surfacePro, functionalPro) #c('CD56', 'GBA1', 'ApoE', 'APP', 'PrP_CD230', 'DAT', 'Tau', 'SERT') c('CD56', 'GBA1', 'CD47')
x <- df_scaled[[3]][1:5000, prop_markers]
colnames(x) <- prop_markers
##Run GateFinder.
ans <- GateFinder(x, targetpop, max.iter=12)
##Make the plots.
png('test.png')
plot.GateFinder(x=x, y=ans, c(4, 3), targetpop=targetpop)#, xlim=c(-3, 3), ylim=c(-3, 3))
dev.off()


sum(df_scaled[[3]][1:10000, 'GFAP'] > 1)

summary(df_scaled[[1]][1:5000, 'GFAP'])



