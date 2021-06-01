"""
This script repeats the analysis process done earlier, but for GFAP-EAAT1- data instead. In this script, this
includes loading the clustered data, assign final cluster assignment through metaclustering.
"""

library(flowCore)
library(FlowSOM)
library(ggpubr)
library(clue)
library(aricode)
library(data.table)

# metaclustering----------------------------------------------------------------------------
file_list <- c(
                'presynTOFGFAPnegEAAT1neg_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOFGFAPnegEAAT1neg_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOFGFAPnegEAAT1neg_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv'
               )
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
write.csv(mc, paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynapticGFAPnegEAAT1neg_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv'), row.names=FALSE)


x <- seq(10)
pairs <- t(combn(x, 2))
Fs <- c()
nmi <- c()
for (row in 1:nrow(pairs)) {
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
mc <- as.data.frame(fread(paste0('R_py_exchange/', 'mcResultsDWH_allGroups_maxK40_allLowNoPresynapticGFAPnegEAAT1neg_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "neglect", "group", "sample"), sep='_')
mc <- mc[, !colnames(mc) %in% c('neglect')]

# for adjustiing column names
new_label <-c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc$mc <- sapply(mc$mc, function(x) new_label[x])


df_region <- list()
adjust_markers <- c('CD47', 'DAT', 'a-Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')

for (region in c('BA9', 'DLCau', 'Hipp')) {
    # loading the presynaptic data with subtraction from the postsynaptic mean
    df_pre <- list()
    for (group in c('LowNo', 'PHAD', 'LBD')) {
            # load presynaptic data
        fcs_path <- '../raw_data/max_events/fcs_GFAPnegEAAT1neg/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        df_pre[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }

        for (i in seq(length(file_list))) {
            frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
            functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                        'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
            surfacePro <- sort(colnames(frame)[!colnames(frame) %in% functionalPro], decreasing=TRUE)
            df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
        }
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

saveRDS(df_region, 'df_pre_save/exp0mc_GFAPnegEAAT1neg_5_13.rds')