"""
This script performs metaclustering on the results of the leave-one-out AE clustering. It also
analyses the similarity of the resulting metaclustering.
"""

library(flowCore)
library(FlowSOM)
library(clue)
library(aricode)
library(data.table)
library(dplyr)
library(tidyr)

########################################################################################
# perform meta clustering #
########################################################################################

# load AE clustering results for (6 of them where each leave one LowNo sample out) 
file_list <- c(
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF13-117_sess_1_for_HF13-117.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF14-076_sess_1_for_HF14-076.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF14-057_sess_1_for_HF14-057.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF14-053_sess_1_for_HF14-053.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF14-051_sess_1_for_HF14-051.csv', 
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_HF14-008_sess_1_for_HF14-008.csv'
               )

# concatenate and perform meta clustering
for (i in seq(length(file_list))) {
    filname <- paste(strsplit(file_list[i], '_')[[1]][5:15], collapse='_')
    cl_mat <- data.frame()
    print(file_list[i])
    cl_mat <- as.data.frame(fread(paste0('R_py_exchange/', file_list[i]), stringsAsFactors=FALSE, header=TRUE)[, -1])

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
}

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



########################################################################################
# perform evaluation #
########################################################################################

# load file after metaclustering
filname <- 'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv'
filname <- paste(strsplit(filname, '_')[[1]][5:15], collapse='_')
mc_all <- as.data.frame(fread(paste0('R_py_exchange/', filname)))
mc_all <- mc_all %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')

# for adjustiing column names
new_label <-c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc_all$mc <- sapply(mc_all$mc, function(x) new_label[x])

# list sample name of control samples
all_control_samples <- c('HF13-117', 'HF14-008', 'HF14-051', 'HF14-053', 'HF14-057', 'HF14-076')

# blank matrix for recording results
r_cl <- matrix(NA, length(all_control_samples), 15, dimnames=list(all_control_samples, seq(15)))
r_m <-  matrix(NA, length(all_control_samples), 25, dimnames=list(all_control_samples, seq(25)))
new_list <- list()

# start evaluation of leave-one-out results
for (sel_sample in all_control_samples) {
    ########
    # get the sample expressions obtained from clustering using ALL samples #
    ########
    df_region_all <- list()
    for (region in c('BA9', 'DLCau', 'Hipp')) {
        df_pre <- list()
        group <- 'LowNo'
        # load presynaptic data
        fcs_path <- '../raw_data/max_events/fcs/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        file_list <- file_list[grepl(sel_sample, file_list)]
        df_pre[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc_all[(mc_all$region==region) & (mc_all$group==group) & grepl(sel_sample, mc_all$sample), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
            functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                        'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
            surfacePro <- sort(colnames(frame)[!colnames(frame) %in% functionalPro], decreasing=TRUE)
            df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
        }
        # attaching cluster_id and sample_id
        for (group in names(df_pre)) {
            df_pre[[group]] <- cbind.data.frame(mc_all[(mc_all$pp=='pre') & (mc_all$region==region) & (mc_all$group==group) & grepl(sel_sample, mc_all$sample), c('mc', 'sample')], df_pre[[group]])
        }
        df_region_all[[region]] <- df_pre
        # remove low quality marker
        df_pre[[group]]['NET'] <- NULL
    }
    for (region in c('BA9', 'DLCau', 'Hipp')){
        # print(dim(df_region_all[[region]][['LowNo']]))
        df_LowNo <- df_region_all[[region]][['LowNo']]
        df_region_all[[region]][['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
        # print(dim(df_region_all[[region]][['LowNo']]))
    }

    # concatenate all regions to 1 dataframe
    df__all <- data.frame()
    for (region in names(df_region_all)) {
        df__all <- rbind.data.frame(df__all, df_region_all[[region]][[group]][, surfacePro])
    }
    # rename columns a bit
    colnames(df__all) <- sapply(colnames(df__all), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                    else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                    else if(x=='PrP_CD230') {'PrP'}
                                                    else if(x=='GATM_1') {'GATM'}
                                                    else if(x=='GAMT_2') {'GAMT'}
                                                    else if(x=='PARKIN') {'Parkin'}
                                                    else if (x=='a-Synuclein') {'AS'}
                                                    else x) 
    mylevels <- colnames(df__all)

    # rearrange matrix and aggregate mean values
    mc_all_sub <- mc_all[(mc_all$pp=='pre') & (grepl(sel_sample, mc_all$sample)), 'mc'] # specify mc that corresponds to the df you loaded
    df_clus_all <- cbind.data.frame(mc_all_sub, df__all)
    aa <- aggregate(df_clus_all[, -1], by=list(df_clus_all$mc_all_sub), FUN=mean)
    aaa <- aa[, sort(colnames(aa[, -1]))]
    rownames(aaa) <- aa[, 1]
    original <- aaa

    ########
    # get the predicted sample expressions obtained from clustering using leave one out method #
    ########
    filname <- paste0('presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_no_', sel_sample, '_sess_1_for_', sel_sample, '.csv')
    filname <- paste(strsplit(filname, '_')[[1]][5:15], collapse='_')
    # load clusters
    mc <- as.data.frame(fread(paste0('R_py_exchange/', filname)))
    mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')

    df_region <- list()
    for (region in c('BA9', 'DLCau', 'Hipp')) {
        # loading the presynaptic data with subtraction from the postsynaptic mean
        df_pre <- list()
        # load presynaptic data
        fcs_path <- '../raw_data/max_events/fcs/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        file_list <- file_list[grepl(sel_sample, file_list)]
        df_pre[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group) & grepl(sel_sample, mc$sample), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
            functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                        'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
            surfacePro <- sort(colnames(frame)[!colnames(frame) %in% functionalPro], decreasing=TRUE)
            df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
        }
        # attaching cluster_id and sample_id
        for (group in names(df_pre)) {
            df_pre[[group]] <- cbind.data.frame(mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group) & grepl(sel_sample, mc$sample), c('mc', 'sample')], df_pre[[group]])
        }
        df_region[[region]] <- df_pre
        # remove low quality marker
        df_pre[[group]]['NET'] <- NULL
    }
    for (region in c('BA9', 'DLCau', 'Hipp')){
        df_LowNo <- df_region[[region]][['LowNo']]
        df_region[[region]][['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    }
    df_ <- df_region
    df__ <- data.frame()
    for (region in names(df_)) {
        df__ <- rbind.data.frame(df__, df_[[region]][[group]][, surfacePro])
    }
    colnames(df__) <- sapply(colnames(df__), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                    else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                    else if(x=='PrP_CD230') {'PrP'}
                                                    else if(x=='GATM_1') {'GATM'}
                                                    else if(x=='GAMT_2') {'GAMT'}
                                                    else if(x=='PARKIN') {'Parkin'}
                                                    else if (x=='a-Synuclein') {'AS'}
                                                    else x) 
    mylevels <- colnames(df__)
    mc_sub <- mc[,'mc']
    df_clus <- cbind.data.frame(mc_sub, df__)
    aa <- aggregate(df_clus[, -1], by=list(df_clus$mc_sub), FUN=mean)
    aaa <- aa[, sort(colnames(aa[, -1]))]
    rownames(aaa) <- aa[, 1]
    new <- aaa


    # calculate assignment and differences ---------------------------------------------------------
    a <- (scale((new)))
    b <- (scale((original)))

    # calculate assignment
    diff <- sapply(seq(1, nrow(a)), function(ii) colSums(abs(a[ii, ] - t(b))))
    # match the new cluster to the closest original's cluster
    assignment <- apply(diff, 1, which.min)

    # handle cases with assignment duplicates or if new one have >15 clusters
    if (any(duplicated(assignment)) | nrow(new)>15) {
        leftover <- seq(1:15)[!seq(1:15) %in% assignment]
        cc <- assignment[duplicated(assignment)]
        all_dup <- assignment[assignment %in% cc]
        assigned <- c()
        for (dup in cc) {
            lowest <- which.min(colSums(abs(a[dup, ])*(abs(a[dup, ] - t(b)))))
            assignment[lowest] <- dup
            assigned <- c(assigned, names(assignment[lowest]))
        }
        used_leftover <- c()
        for (l in names(all_dup)[!names(all_dup) %in% assigned]) {
            lowest <- leftover[which.min(diff[l, leftover])]
            assignment[l] <- lowest
            used_leftover <- c(used_leftover, lowest)
        }
        if (length(c(leftover[!leftover %in% used_leftover])) > 0) {
            for (l in c(leftover[!leftover %in% used_leftover])) {
                diff_in <- colSums(sapply( seq(1, nrow(a)), function(ii) (abs(a[ii, ] - a[l, ])) ))
                closest_cl <- which.min(diff_in[diff_in != 0])
                new[closest_cl, ] <- (new[closest_cl, ] + new[l, ])*0.5
                new <- new[!rownames(new) %in% l, ]
            }
        }
    }
    a <- (scale((new)))
    new <- new[as.character(assignment), ]

    # calculate R by cluster
    if (sel_sample==all_control_samples[1]) {
        clnames_res <- names(assignment)
        mnames_res <- colnames(new)
    } else {
        if (any(names(assignment) != clnames_res) | any(colnames(new) != mnames_res) | any(colnames(original) != mnames_res)) {
            print('shitttttt')
        }
    }
    r_cl[sel_sample, ] <- sapply(names(assignment), function(cl) cor(as.numeric(new[as.character(assignment[cl]), ]), as.numeric(original[cl, ])))

    # calculate R by marker
    r_m[sel_sample, ] <- sapply(colnames(new), function(m) cor(as.numeric(new[, m]), as.numeric(original[, m])))
    new_list[[sel_sample]] <- new
}

# name the result by cluster/marker names
colnames(r_cl) <- names(assignment)
colnames(r_m) <- colnames(new)

# get the mean of the marker results
mean_mat <- Reduce("+", new_list)/length(new_list)
rownames(mean_mat) <- names(assignment)
mean_mat <- mean_mat[c(rownames(mean_mat)[1:5], rownames(mean_mat)[8:15], rownames(mean_mat)[6:7]), ]

# plot out average heat map results
redgreen <- c("black", "#20114bff", "#1c1043ff", '#56147dff', '#d6446dff', "#fea16eff", '#fcf7b9ff')
pal <- colorRampPalette(redgreen)(100)
pdf(paste0('figures/heatmaps/meanAll.pdf'))
a <- heatmap.2(as.matrix(t(mean_mat)), Rowv=FALSE, Colv=FALSE, col=pal,
          scale='row', trace='none', srtCol=90, adjCol = c(0.75, 0.5),
          lhei=c(2, 12), lwid=c(2, 5), key.par=list(cex=0.622), density.info="none", 
          key.title='', keysize=0.1, margins=c(3,7))
dev.off()

# evaluate results
# NMI
sapply(colnames(new), function(m) NMI(as.numeric(new[, m]), as.numeric(original[, m])))
# effect size
eff_size <- sapply(names(assignment), function(cl) ((abs(new[as.character(assignment[cl]), ] - original[cl, ]))/t(abs(new[as.character(assignment[cl]), ] + original[cl, ]))*2))
eff_size <- matrix(unlist(eff_size), ncol=ncol(eff_size), byrow=FALSE, dimnames=list(rownames(eff_size), colnames(eff_size)))
avg_eff_size <- colSums(as.matrix(eff_size))/nrow(eff_size)*100