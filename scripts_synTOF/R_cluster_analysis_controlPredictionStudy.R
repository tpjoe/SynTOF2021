library(flowCore)
library(FlowSOM)
library(ggpubr)
library(clue)
library(aricode)
library(data.table)

# metaclustering----------------------------------------------------------------------------
file_list <- c(
                # 'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                # 'presynTOF_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                # 'presynTOF_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_H,F,1,3,-,1,1,7_sess_1_for_H,F,1,3,-,1,1,7.csv'
                # 'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1_for_H,F,1,3,-,1,1,7.csv'
               )
filname <- paste(strsplit(file_list[1], '_')[[1]][5:15], collapse='_')


cl_mat <- data.frame()
for (i in seq(length(file_list))) {
    print(file_list[i])
    cl_mat_ <- as.data.frame(fread(paste0('R_py_exchange/', file_list[i]), stringsAsFactors=FALSE, header=TRUE)[, -1])
    cl_mat <- rbind.data.frame(cl_mat, cl_mat_)
}

AllClusters <- list()
for (ii in 1:(ncol(cl_mat) -1)){
    AllClusters[[ii]] <- as.cl_partition(as.vector(cl_mat[, ii]))
}
##use the CLUE package to create the consensus partitions
a <- cl_ensemble(list=AllClusters)
b <- cl_consensus(a, method='DWH')$.Data
mc <- apply(b, 1, which.max)
mc <- cbind.data.frame(mc, sample=cl_mat$sample)
# save meta results
write.csv(mc, paste0('R_py_exchange/', filname), row.names=FALSE)
# write.csv(mc, paste0('R_py_exchange/test', 'DLCau', '.csv'), row.names=FALSE)


x <- seq(10)
pairs <- t(combn(x, 2))
Fs <- c()
nmi <- c()
for (row in 1:nrow(pairs)) {
    print(paste(pairs[row, 1], pairs[row, 2]))
    Fs <- c(Fs, FMeasure(cl_mat[, pairs[row, 1]], cl_mat[, pairs[row, 2]]))
    nmi <- c(nmi, NMI(cl_mat[, pairs[row, 1]], cl_mat[, pairs[row, 2]]))
}
mean(Fs)
mean(nmi)


# load corrected pre-synaptic data ----------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)

# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/', filname)))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')


# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc$mc <- sapply(mc$mc, function(x) new_label[x])


df_region <- list()
adjust_markers <- c('CD47', 'DAT', 'a-Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')

for (region in c('BA9', 'DLCau', 'Hipp')) {
    # load data to get post-synaptic reference
    # df_post <- list()
    # for (group in c('LowNo', 'LBD', 'PHAD')) {
    #     # load postsynaptic data
    #     fcs_path <- '../raw_data/max_events/fcs_post_synap/'
    #     file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
    #     df_post[[group]] <- data.frame()
    #     # check that sample order align with mc
    #     df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
    #     if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
    #         print('shitttttttttttttttttttttttttttttttt')
    #     }
    #     for (i in seq(length(file_list))) {
    #         frame <- read.FCS(paste0(fcs_path, file_list[i]))
    #         df_post[[group]] <- rbind.data.frame(df_post[[group]], exprs(frame))
    #         file_name_splt <- strsplit(file_list[i], '_')[[1]]
    #     }
    #     df_post[[group]]['NET'] <- NULL
    # }

    functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
                    #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1') #possible
    surfacePro <- sort(colnames(df_post[[group]])[!colnames(df_post[[group]]) %in% functionalPro], decreasing=TRUE)

    # attaching cluster_id and sample_id
    # for (group in names(df_post)) {
    #     df_post[[group]] <- cbind.data.frame(mc[(mc$pp=='post') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_post[[group]])
    # }

    # df_post_ref <- list()
    # # get mean of the post
    # for (group in names(df_post)) {
    #     df_post_ref[[group]] <- df_post[[group]] %>% group_by(sample, mc) %>% summarise_at(vars(c(functionalPro, surfacePro)), list(mean=mean))
    #     df_post_ref[[group]] <- df_post_ref[[group]] %>% pivot_wider(names_from=mc, 
    #                             values_from=one_of(paste0(c(functionalPro, surfacePro), '_mean')))
    #     df_post_ref[[group]] <- cbind.data.frame(group=rep(group, nrow(df_post_ref[[group]])), df_post_ref[[group]])
    # }

    # loading the presynaptic data with subtraction from the postsynaptic mean
    df_pre <- list()
    group <- 'LowNo'

    # load presynaptic data
    fcs_path <- '../raw_data/max_events/fcs/'
    file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
    # file_list <- file_list[!grepl('HF14-017', file_list) & !grepl('HF14-083', file_list) & !grepl('HF14-025', file_list) & !grepl('HF13-117', file_list) & !grepl('HF14-008', file_list)]
    file_list <- file_list[grepl('HF13-117', file_list) | grepl('HF14-008', file_list)]
    df_pre[[group]] <- data.frame()
    # check that sample order align with mc
    df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
    if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
        print('shitttttttttttttttttttttttttttttttt')
    }

    for (i in seq(length(file_list))) {
        frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
    #         # subtraction begins here
    #         sample <- strsplit(file_list[i], '_')[[1]][4]
    #         mc_sample <- mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group) & (mc$sample==sample), ]
    #         post_strsplit <- strsplit(colnames(df_post_ref[[group]][df_post_ref[[group]][, 'sample']==sample, ]), '_')
    #         post_cl <- sapply(post_strsplit, function(x) x[length(x)])
    #         post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
    #         for (marker in adjust_markers) {
    #             uni_cl <- sort(unique(post_cl))
    #             # for (cl in uni_cl[!uni_cl %in% c('group', 'sample')]) {
    #             for (cl in new_label[c(5, 13)]) {
    #                 # ref <- df_post_ref[[group]][df_post_ref[[group]][, 'sample']==sample, (post_marker==marker) & (post_cl==cl)]
    #                 post_strsplit <- strsplit(colnames(df_post_ref[['LowNo']][df_post_ref[['LowNo']][, 'sample']==sample, ]), '_')
    #                 post_cl <- sapply(post_strsplit, function(x) x[length(x)])
    #                 post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
    #                 ref_LowNo <- tryCatch(df_post_ref[['LowNo']][, (post_marker==marker) & (post_cl==cl)],
    #                                       error=function(x) {0})
    #                 if (length(ref_LowNo)==0) {ref_LowNo<-0}
    #                 post_strsplit <- strsplit(colnames(df_post_ref[['LBD']][df_post_ref[['LBD']][, 'sample']==sample, ]), '_')
    #                 post_cl <- sapply(post_strsplit, function(x) x[length(x)])
    #                 post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
    #                 ref_LBD <- tryCatch(df_post_ref[['LBD']][, (post_marker==marker) & (post_cl==cl)],
    #                                     error=function(x) {0})
    #                 if (length(ref_LBD)==0) {ref_LBD<-0}
    #                 post_strsplit <- strsplit(colnames(df_post_ref[['PHAD']][df_post_ref[['PHAD']][, 'sample']==sample, ]), '_')
    #                 post_cl <- sapply(post_strsplit, function(x) x[length(x)])
    #                 post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
    #                 ref_PHAD <- tryCatch(df_post_ref[['PHAD']][, (post_marker==marker) & (post_cl==cl)],
    #                                       error=function(x) {0})
    #                 if (length(ref_PHAD)==0) {ref_PHAD<-0}
    #                 mean_ref <- c()
    #                 for (r in c(list(ref_LowNo), list(ref_LBD), list(ref_PHAD))) {
    #                     # if (sum(!is.na(r))/length(r) > 0.5) {
    #                         mean_ref <- c(mean_ref, mean(r, na.rm=TRUE))
    #                     # }
    #                 }
    #                 # if (length(mean_ref) > 0){
    #                     ref <- mean(mean_ref, na.rm=TRUE)
    #                     if (!is.na(ref)) {
    #                         frame[mc_sample$mc==cl, marker] <- frame[mc_sample$mc==cl, marker] - ref
    #                         # frame[mc_sample$mc==cl, marker][frame[mc_sample$mc==cl, marker]<0] <- rep(0, length(
    #                             # frame[mc_sample$mc==cl, marker][frame[mc_sample$mc==cl, marker]<0])) # flooring
    #                     }
    #                 # }
    #             }
    #         }
        df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
    }
    # attaching cluster_id and sample_id
    for (group in names(df_pre)) {
        df_pre[[group]] <- cbind.data.frame(mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_pre[[group]])
    }
    df_region[[region]] <- df_pre
    df_pre[[group]]['NET'] <- NULL
}


for (region in c('BA9', 'DLCau', 'Hipp')){
    print(dim(df_region[[region]][['LowNo']]))
    df_LowNo <- df_region[[region]][['LowNo']]
    df_region[[region]][['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    print(dim(df_region[[region]][['LowNo']]))
}

# saveRDS(df_region, 'df_pre_save/exp2mc_5_13.rds')



# Create freq. and exprs matrix ------------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)


# for (exp_no in c(2)){ #c(0, seq(6))) {
# df_region <- readRDS(paste0('df_pre_save/exp', exp_no, 'mc_5_13.rds'))
for (region in names(df_region)) {
    df_mean <- list()
    df_freq <- list()
    for (group in names(df_region[[region]])) {
    # reorder the column
    for_level <- sort(unique(df_region[[region]][[group]][, 'mc']))
    for_level <- c(for_level[1:5], for_level[8:15], for_level[6:7])
    df_region[[region]][[group]][, 'mc'] <- factor(df_region[[region]][[group]][, 'mc'], levels=for_level, ordered=TRUE)
    # start grouping the mean
    functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
    surfacePro <- sort(colnames(df_region[['BA9']][[group]])[!colnames(df_region[['BA9']][[group]]) %in% 
                        c('mc', 'sample', functionalPro)], decreasing=TRUE)
    df_mean[[group]] <- df_region[[region]][[group]] %>% group_by(sample, mc) %>% summarise_at(vars(c(functionalPro, surfacePro)), list(mean=mean))
    df_mean[[group]] <- df_mean[[group]] %>% pivot_wider(names_from=mc, 
                        values_from=one_of(paste0(c(functionalPro, surfacePro), '_mean')))
    df_mean[[group]] <- cbind.data.frame(group=rep(group, nrow(df_mean[[group]])), df_mean[[group]])
    df_freq[[group]] <- df_region[[region]][[group]] %>% group_by(sample, mc) %>% tally() %>% 
                        pivot_wider(names_from=mc, values_from=n)
    nonNA <- which(sapply(1:ncol(df_freq[[group]]), function(x) any(is.na(df_freq[[group]][, x]))))
    df_freq[[group]][, c(-1, -nonNA)] <- df_freq[[group]][, c(-1, -nonNA)]/rowSums(df_freq[[group]][, c(-1, -nonNA)])
    df_freq[[group]] <- cbind.data.frame(group=rep(group, nrow(df_freq[[group]])), df_freq[[group]])
    # colnames(df_freq[[group]]) <- c(colnames(df_freq[[group]])[1:2], label_map[[region]])
    }
    # export the final tables
    df_freq_region <- do.call(rbind.data.frame, df_freq)
    df_mean_region <- do.call(rbind.data.frame, df_mean)
    # write.csv(df_freq_region, paste0('R_py_exchange/df_freq_', region, '_noStd.csv'))
    # write.csv(df_mean_region, paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp', exp_no, 'mc_5_13.csv'))
}
# }



# load the one trained by 6 samples
df_mean <- df_mean_region[[1]]
df_mean_region_all <- read.csv(paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp', 0, 'mc_5_13.csv'))[, -1]
colnames(df_mean_region_all) <- gsub('DJ.1_PARK7', 'DJ\\-1_PARK7', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('a.Synuclein', 'a\\-Synuclein', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('b.Amyloid_X40', 'b\\-Amyloid_X40', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('b.Amyloid_X42', 'b\\-Amyloid_X42', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('p.Tau', 'p\\-Tau', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('K48.Ubiquitin', 'K48\\-Ubiquitin', colnames(df_mean_region_all))
colnames(df_mean_region_all) <- gsub('X3NT', '3NT', colnames(df_mean_region_all))
df_mean_region_all <- df_mean_region_all[, Reduce("&", lapply(c(functionalPro, 'NET'), function(x) !(grepl(x, colnames(df_mean_region_all)))))]
df_mean <- df_mean[, Reduce("&", lapply(c(functionalPro, 'NET'), function(x) !(grepl(x, colnames(df_mean)))))]

allsix <- df_mean_region_all[df_mean_region_all$sample %in% df_mean_region$sample, 3:ncol(df_mean)] 


(allsix[, 1] - df_mean[, 3:ncol(df_mean)][, 1])/df_mean[, 3:ncol(df_mean)][, 1]


allsix[, 1]
df_mean[, 3:ncol(df_mean)][, 1]


# Expression visualization and heatmaps -------------------------------------------------------------------------
# load data 
library(dplyr)
library(tidyr)
library(ggplot2)
library(flowCore)
library(data.table)
library(reshape2)
library(gplots)


mc <- as.data.frame(fread(paste0('R_py_exchange/', filname)))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
df_region <- readRDS(paste0('df_pre_save/exp', 2, 'mc_5_13.rds'))




df_ <- df_region
group <- 'LowNo'
# region <- 'DLCau'

functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
surfacePro <- sort(colnames(df_[['BA9']][[group]])[!colnames(df_[['BA9']][[group]]) %in% 
                   c('mc', 'sample', functionalPro)], decreasing=TRUE)

df__ <- data.frame()
for (region in names(df_)) {
    df__ <- rbind.data.frame(df__, df_[[region]][[group]][, surfacePro])
}
# df__ <- df_[[region]][[group]][, surfacePro]
# df2 <- apply(df__, 2, function(x) (x-min(x))/(max(x)-min(x)))
colnames(df__) <- sapply(colnames(df__), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                   else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                   else if(x=='PrP_CD230') {'PrP'}
                                                   else if(x=='GATM_1') {'GATM'}
                                                   else if(x=='GAMT_2') {'GAMT'}
                                                   else if(x=='PARKIN') {'Parkin'}
                                                   else if (x=='a-Synuclein') {'AS'}
                                                   else x) 
mylevels <- colnames(df__)


mc_sub <- mc[,'mc']# mc[(mc$pp=='pre') & (mc$group==group), 'mc'] # specify mc that corresponds to the df you loaded

# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc_sub <- sapply(mc_sub, function(x) new_label[x])
df_clus <- cbind.data.frame(mc_sub, df__)

aa <- aggregate(df_clus[, -1], by=list(df_clus$mc_sub), FUN=mean)
aaa <- aa[, sort(colnames(aa[, -1]))]
rownames(aaa) <- aa[, 1]

redgreen <- c("#13192b", "#1c2541", "#3a506b", '#5bc0be', '#e4fde1', 'white')
pal <- colorRampPalette(redgreen)(100)
pdf(paste0('figures/heatmaps/all', '_', group, '_meanNoStd_no_HF13-117,HF14-008.pdf'))
heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=as.dendrogram(hclust(dist((as.matrix(aaa))))), col=pal,
# heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=FALSE, col=pal,
          scale='row', trace='none', srtCol=90, adjCol = c(0.75, 0.5),
          lhei=c(2, 12), lwid=c(2, 5), key.par=list(cex=0.622), density.info="none", 
          key.title='', keysize=0.1, margins=c(3,7))
dev.off()


a <- t(scale(t(aaa)))
b <- t(scale(t(for2)))


abs(a['B1', ] - b['B1', ]) / (abs(a['B1', ])/2 + abs(b['B1', ])/2)





