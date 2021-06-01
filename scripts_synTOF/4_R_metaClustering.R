library(flowCore)
library(FlowSOM)
library(ggpubr)
library(clue)
library(aricode)
library(data.table)
library(viridis)

# metaclustering----------------------------------------------------------------------------
# load cluster results for each group in pre- and post-synaptic data
file_list <- c(
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOF_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOF_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'postsynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'postsynTOF_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'postsynTOF_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv'
               )

# concatenate results row-wise
cl_mat <- data.frame()
for (i in seq(length(file_list))) {
    print(file_list[i])
    cl_mat_ <- as.data.frame(fread(paste0('R_py_exchange/', file_list[i]), stringsAsFactors=FALSE, header=TRUE)[, -1])
    cl_mat <- rbind.data.frame(cl_mat, cl_mat_)
}

# assign as cluster partition
AllClusters <- list()
for (ii in 1:(ncol(cl_mat) -1)){
    AllClusters[[ii]] <- as.cl_partition(as.vector(cl_mat[, ii]))
}

# use the CLUE package to create the consensus partitions
a <- cl_ensemble(list=AllClusters)
b <- cl_consensus(a, method='DWH')$.Data
mc <- apply(b, 1, which.max)
mc <- cbind.data.frame(mc, sample=cl_mat$sample)

# save meta results
write.csv(mc, paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv'), row.names=FALSE)

# check cluster similarity
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