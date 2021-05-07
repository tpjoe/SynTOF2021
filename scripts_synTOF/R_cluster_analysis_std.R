library(flowCore)
library(FlowSOM)
library(ggpubr)
library(clue)
library(aricode)
library(data.table)

# metaclustering----------------------------------------------------------------------------
file_list <- c(
#                'presynTOF_predLowNo_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv',
#                'presynTOF_predLBD_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv',
#                'presynTOF_predPHAD_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv',
#                'postsynTOF_predLowNo_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv',
#                'postsynTOF_predLBD_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv',
                'presynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv',
                'presynTOF_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv',
                'presynTOF_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv',
                'postsynTOF_AdamMegaAEpredLowNo_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv',
                'postsynTOF_AdamMegaAEpredLBD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv',
                'postsynTOF_AdamMegaAEpredPHAD_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv'
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
write.csv(mc, paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv'), row.names=FALSE)
# write.csv(mc, paste0('R_py_exchange/test', 'DLCau', '.csv'), row.names=FALSE)


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



# plot hidden representaion with clusters assignment ---------------------------------------------------------
library(uwot)
library(ggpubr)
library(Rtsne)
library(data.table)


# load individual clusters---
mc <- as.data.frame(fread('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv'))
hidden <- as.data.frame(fread('R_py_exchange/hidden_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv', nrows=100000)[-1])
# hidden_scaled <- cbind.data.frame(scale(hidden[, -ncol(hidden)]), hidden[, ncol(hidden)])
# quantile <- apply(hidden[, 1:(ncol(hidden)-1)], 2, function(x) quantile(x, 0.9999))
# hidden_quantile <- t(t(hidden[, 1:(ncol(hidden)-1)])/quantile)

# sample for visual
sampled <- sort(sample(1:nrow(hidden), 30000, replace=FALSE))
hidden_ <- hidden[sampled, -ncol(hidden)]
mc_ <- factor(mc[sampled, 1])
# mc_ <- factor(cl_mat[sampled, 1])
umap_xy <- uwot::umap(hidden_, n_neighbors=100, metric='cosine', min_dist=0.2, scale=TRUE)#, spread=0.5)#, 
plot_df <- cbind.data.frame(umap_xy, mc_)
colnames(plot_df) <- c('x', 'y', 'mc')

p_mc <- ggscatter(plot_df, x="x", y="y", color = "mc", size=0.1, xlab='', ylab='') + theme_bw(base_size=15) + 
        theme(legend.position="top", axis.text.x = element_blank(), axis.text.y = element_blank(), 
        axis.ticks = element_blank()) + #, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        guides(colour=guide_legend(override.aes = list(size=8), title="Predicted cluster"))
set_palette(p_mc, "springfield_simpsons")
ggsave('hidden.pdf', p_mc)




# load corrected pre-synaptic data ----------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)

# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')


# for adjustiing column names
new_label <- c('GAD65', 'Synapto', 'Mixed-Low', 'High', 'ApoE', 
                'SLC6A8', 'GAMT', 'Low', 'Parkin', 'GBA', 
                'LRRK2', 'DAT', 'Mixed-High', 'Calret.', 'SERT')
mc$mc <- sapply(mc$mc, function(x) new_label[x])


df_region <- list()
adjust_markers <- c('CD47', 'DAT', 'a-Synuclein', 'VGLUT', 'GAD65', 'VMAT2', 'Synaptobrevin2')

for (region in c('BA9', 'DLCau', 'Hipp')) {
    # load data to get post-synaptic reference
    df_post <- list()
    for (group in c('LowNo', 'LBD', 'PHAD')) {
        # load postsynaptic data
        fcs_path <- '../raw_data/max_events/fcs_post_synap/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        df_post[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- read.FCS(paste0(fcs_path, file_list[i]))
            df_post[[group]] <- rbind.data.frame(df_post[[group]], exprs(frame))
            file_name_splt <- strsplit(file_list[i], '_')[[1]]
        }
        df_post[[group]]['NET'] <- NULL
    }

    functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
                    #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1') #possible
    surfacePro <- sort(colnames(df_post[[group]])[!colnames(df_post[[group]]) %in% functionalPro], decreasing=TRUE)

    # attaching cluster_id and sample_id
    for (group in names(df_post)) {
        df_post[[group]] <- cbind.data.frame(mc[(mc$pp=='post') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_post[[group]])
    }

    df_post_ref <- list()
    # get mean of the post
    for (group in names(df_post)) {
        df_post_ref[[group]] <- df_post[[group]] %>% group_by(sample, mc) %>% summarise_at(vars(c(functionalPro, surfacePro)), list(mean=mean))
        df_post_ref[[group]] <- df_post_ref[[group]] %>% pivot_wider(names_from=mc, 
                                values_from=one_of(paste0(c(functionalPro, surfacePro), '_mean')))
        df_post_ref[[group]] <- cbind.data.frame(group=rep(group, nrow(df_post_ref[[group]])), df_post_ref[[group]])
    }

    # loading the presynaptic data with subtraction from the postsynaptic mean
    df_pre <- list()
    for (group in c('LowNo', 'LBD', 'PHAD')) {
        # load presynaptic data
        fcs_path <- '../raw_data/max_events/fcs/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        df_pre[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
            # subtraction begins here
            sample <- strsplit(file_list[i], '_')[[1]][4]
            mc_sample <- mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group) & (mc$sample==sample), ]
            post_strsplit <- strsplit(colnames(df_post_ref[[group]][df_post_ref[[group]][, 'sample']==sample, ]), '_')
            post_cl <- sapply(post_strsplit, function(x) x[length(x)])
            post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
            for (marker in adjust_markers) {
                uni_cl <- sort(unique(post_cl))
                # for (cl in uni_cl[!uni_cl %in% c('group', 'sample')]) {
                for (cl in new_label[c(7, 8)]) {
                    # ref <- df_post_ref[[group]][df_post_ref[[group]][, 'sample']==sample, (post_marker==marker) & (post_cl==cl)]
                    post_strsplit <- strsplit(colnames(df_post_ref[['LowNo']][df_post_ref[['LowNo']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_LowNo <- tryCatch(df_post_ref[['LowNo']][, (post_marker==marker) & (post_cl==cl)],
                                          error=function(x) {0})
                    if (length(ref_LowNo)==0) {ref_LowNo<-0}
                    post_strsplit <- strsplit(colnames(df_post_ref[['LBD']][df_post_ref[['LBD']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_LBD <- tryCatch(df_post_ref[['LBD']][, (post_marker==marker) & (post_cl==cl)],
                                        error=function(x) {0})
                    if (length(ref_LBD)==0) {ref_LBD<-0}
                    post_strsplit <- strsplit(colnames(df_post_ref[['PHAD']][df_post_ref[['PHAD']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_PHAD <- tryCatch(df_post_ref[['PHAD']][, (post_marker==marker) & (post_cl==cl)],
                                          error=function(x) {0})
                    if (length(ref_PHAD)==0) {ref_PHAD<-0}
                    mean_ref <- c()
                    for (r in c(list(ref_LowNo), list(ref_LBD), list(ref_PHAD))) {
                        # if (sum(!is.na(r))/length(r) > 0.5) {
                            mean_ref <- c(mean_ref, mean(r, na.rm=TRUE))
                        # }
                    }
                    # if (length(mean_ref) > 0){
                        ref <- mean(mean_ref, na.rm=TRUE)
                        if (!is.na(ref)) {
                            frame[mc_sample$mc==cl, marker] <- frame[mc_sample$mc==cl, marker] - ref
                            # frame[mc_sample$mc==cl, marker][frame[mc_sample$mc==cl, marker]<0] <- rep(0, length(
                                # frame[mc_sample$mc==cl, marker][frame[mc_sample$mc==cl, marker]<0])) # flooring
                        }
                    # }
                }
            }
            df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
        }
        df_pre[[group]]['NET'] <- NULL
    }
    # attaching cluster_id and sample_id
    for (group in names(df_pre)) {
        df_pre[[group]] <- cbind.data.frame(mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_pre[[group]])
    }
    df_region[[region]] <- df_pre
}


for (region in c('BA9', 'DLCau', 'Hipp')){
    print(dim(df_region[[region]][['LowNo']]))
    df_LowNo <- df_region[[region]][['LowNo']]
    df_region[[region]][['ODC']] <- df_LowNo[df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    df_region[[region]][['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    print(dim(df_region[[region]][['LowNo']]))
    print(dim(df_region[[region]][['ODC']]))
}

# saveRDS(df_region, 'df_pre_save/exp2mc_5_13_std.rds')



# Create freq. and exprs matrix ------------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)


for (exp_no in c(2)){ #c(0, seq(6))) {
    # grouping data by IDs
    # exp_no <- 2
    df_region <- readRDS(paste0('df_pre_save/exp', exp_no, 'mc_5_13_std.rds'))
    for (region in names(df_region)) {
        df_mean <- list()
        df_freq <- list()
        for (group in names(df_region[[region]])) {
        # reorder the column
        for_level <- sort(unique(df_region[[region]][[group]][, 'mc']))
        for_level <- c(for_level[1:5], for_level[8:15], for_level[6:7])
        df_region[[region]][[group]][, 'mc'] <- factor(df_region[[region]][[group]][, 'mc'], levels=for_level, ordered=TRUE)
        ncol_ <- ncol(df_region[[region]][[group]])
        if (group=='LowNo') {
            sd_lowno <- apply(df_region[[region]][['LowNo']][, c(3:ncol_)], 2, sd)
            df_region[[region]][[group]][, c(3:ncol_)] <- scale(df_region[[region]][[group]][, c(3:ncol_)], center=FALSE)
        } else {
            df_region[[region]][[group]][, c(3:ncol_)] <- scale(df_region[[region]][[group]][, c(3:ncol_)], center=FALSE, scale=sd_lowno)
        }
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
        write.csv(df_freq_region, paste0('R_py_exchange/df_freq_', region, '_Std.csv'))
        write.csv(df_mean_region, paste0('R_py_exchange/df_meanAllMarkers_', region, '_Std_exp', exp_no, 'mc_5_13.csv'))
    }
}



library('ggpubr')
library('mixdist')

region <- 'Hipp'
sel_marker <- 'GAD65'
sel_clus <- 'A2'
group <- 'LowNo'
aa <- df_region[[region]][[group]]
colnames(aa) <- gsub('-', '.', colnames(aa))


plt <- gghistogram(aa[(aa$mc==sel_clus), c('mc', 'sample', sel_marker), ], x=sel_marker, facet.by='sample',
        add = "mean", rug = FALSE, add_density = FALSE) + ggtitle(paste(group, 'cluster: ',sel_clus))


ggsave('test2.pdf', plt)


# calculate mid local minima
values <- c()
for (samp in unique(aa$sample)) {
    bb <- aa[aa$sample==samp, ]
    
    # find local minima point
    his <- hist(bb[, sel_marker], breaks=101)  
    df <- data.frame(mid=his$mids, cou=his$counts)  
    fitpro <- mix(as.mixdata(df), mixparam(mu=c(0.8, 2.1), sigma=c(1, 1)), dist='norm')
    density <- density(bb[, sel_marker])
    local_min <- min(density$y[(density$x>fitpro$parameter[1, 'mu']) & (density$x<fitpro$parameter[2, 'mu'])])
    local_min_pos <- density$x[which.min(abs(density$y-local_min))]

    result <- mean(bb[bb[, sel_marker] > local_min_pos, sel_marker])
    values <- c(values, result)
}
mean(values)



# load only post-synaptic data ----------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)


# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')

# for adjustiing column names
new_label <- c('GAD65', 'Synapto', 'Mixed-Low', 'High', 'ApoE', 
                'SLC6A8', 'GAMT', 'Low', 'Parkin', 'GBA', 
                'LRRK2', 'DAT', 'Mixed-High', 'Calret.', 'SERT')
mc$mc <- sapply(mc$mc, function(x) new_label[x])

# for reordering colum names later
desired_order <- sort(new_label)[c(4, 13, 8, 3, 5, 14, 12, 1, 7, 10, 11, 9, 15, 6, 2)]

df_mean <- list()
df_freq <- list()
for (region in c('BA9', 'DLCau', 'Hipp')) {
    # load data to get post-synaptic reference
    df_post <- list()
    for (group in c('LowNo', 'LBD', 'PHAD')) {
        # load postsynaptic data
        fcs_path <- '../raw_data/max_events/fcs_post_synap/'
        file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
        df_post[[group]] <- data.frame()
        # check that sample order align with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- read.FCS(paste0(fcs_path, file_list[i]))
            df_post[[group]] <- rbind.data.frame(df_post[[group]], exprs(frame))
            file_name_splt <- strsplit(file_list[i], '_')[[1]]
        }
        df_post[[group]]['NET'] <- NULL
    }

    functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
                    #    'PARKIN', 'TMEM230_C20orf30', 'DJ-1_PARK7', 'GBA1') #possible
    surfacePro <- sort(colnames(df_post[[group]])[!colnames(df_post[[group]]) %in% functionalPro], decreasing=TRUE)

    # attaching cluster_id and sample_id
    for (group in names(df_post)) {
        df_post[[group]] <- cbind.data.frame(mc[(mc$pp=='post') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_post[[group]])
    }
    
    # removing fake LowNo
    df_LowNo <- df_post[['LowNo']]
    df_post[['ODC']] <- df_LowNo[df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    df_post[['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]

    df_post_ref <- list()
    df_freq <- list()
    # get mean of the post
    for (group in names(df_post)) {
        df_post_ref[[group]] <- df_post[[group]] %>% group_by(sample, mc) %>% summarise_at(vars(c(functionalPro, surfacePro)), list(mean=mean))
        df_post_ref[[group]] <- df_post_ref[[group]] %>% pivot_wider(names_from=mc, 
                                values_from=one_of(paste0(c(functionalPro, surfacePro), '_mean')))
        df_post_ref[[group]] <- cbind.data.frame(group=rep(group, nrow(df_post_ref[[group]])), df_post_ref[[group]])
        df_freq[[group]] <- df_post[[group]] %>% group_by(sample, mc) %>% tally() %>% pivot_wider(names_from=mc, values_from=n)
        nonNA <- which(sapply(1:ncol(df_freq[[group]]), function(x) any(is.na(df_freq[[group]][, x]))))
        # df_freq[[group]][, c(-1, -nonNA)] <- df_freq[[group]][, c(-1, -nonNA)]/rowSums(df_freq[[group]][, c(-1, -nonNA)])
        df_freq[[group]] <- cbind.data.frame(group=rep(group, nrow(df_freq[[group]])), df_freq[[group]])
    }
    df_freq <- do.call(bind_rows, df_freq)
    # df_freq <- df_freq[, c('group', 'sample', sort(as.numeric(colnames(df_freq)[3:1000L])))]
    df_freq[, 3:ncol(df_freq)] <- t(apply(df_freq[, 3:ncol(df_freq)], 1, function(x) x/sum(x, na.rm=TRUE)))
    mean_freq <- apply(df_freq[, 3:ncol(df_freq)], 2, function(x) mean(x, na.rm=TRUE)*100)
    df_mean_region <- do.call(bind_rows, df_post_ref)

     # reorder the column
    for_level <- colnames(df_freq)[!colnames(df_freq) %in% c('group', 'sample')]
    missing_columns <- new_label[!new_label %in% for_level]

     # filling missing columns (NA for every row) with new columns of 0's
    df_freq <- cbind.data.frame(df_freq, setNames(data.frame(matrix(0, nrow(df_freq), length(missing_columns))), missing_columns))
    df_freq <- df_freq[, c(colnames(df_freq[1:2]), desired_order)]

    # write.csv(df_freq, paste0('R_py_exchange/df_freqPost_', region, '_Std.csv'))
    # write.csv(df_mean_region, paste0('R_py_exchange/df_meanPostMarkers_', region, '_Std.csv'))
}



DAT_col <- grepl('DAT', colnames(df_mean_region))
df_mean_region[, DAT_col]




# Expression visualization and heatmaps -------------------------------------------------------------------------
# load data 
library(dplyr)
library(tidyr)
library(ggplot2)
library(flowCore)
library(data.table)
library(reshape2)
library(gplots)


mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_Std_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
mc[mc$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), 'group'] <- 'ODC'
df_region <- readRDS(paste0('df_pre_save/exp', 2, 'mc_5_13_std.rds'))

df_ <- df_region
group <- 'LowNo'
# region <- 'Hipp'

functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
surfacePro <- sort(colnames(df_[['BA9']][[group]])[!colnames(df_[['BA9']][[group]]) %in% 
                   c('mc', 'sample', functionalPro)], decreasing=TRUE)

df__ <- data.frame()
for (region in names(df_)) {
    df__ <- rbind.data.frame(df__, df_[[region]][[group]][, surfacePro])
}
# df__ <- df_[[region]][[group]][, surfacePro]
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


# mc_sub <- mc[(mc$pp=='pre') & (mc$group==group) & (mc$region==region), 'mc'] # specify mc that corresponds to the df you loaded
mc_sub <- mc[(mc$pp=='pre') & (mc$group==group), 'mc'] # specify mc that corresponds to the df you loaded


# for adjustiing column names
new_label <- c('GAD65', 'Synapto.', 'Mixed-Low', 'High', 'ApoE', 'SLC6A8', 'GAMT', 'Low', 
               'Parkin', 'GBA', 'LRRK2', 'DAT', 'Mixed-High', 'Calret.', 'SERT')
mc_sub <- sapply(mc_sub, function(x) new_label[x])
df_clus <- cbind.data.frame(mc_sub, df__)

aa <- aggregate(df_clus[, -1], by=list(df_clus$mc_sub), FUN=mean)
aaa <- aa[, sort(colnames(aa[, -1]))]
rownames(aaa) <- aa[, 1]

redgreen <- c("#13192b", "#1c2541", "#3a506b", '#5bc0be', '#e4fde1', 'white')
pal <- colorRampPalette(redgreen)(100)
aaa <- aaa[c("High", "Mixed-High", "Low", "Mixed-Low", "ApoE", "Calret.", "DAT", "GAD65", "GAMT", "GBA", "LRRK2", "Parkin", "SERT", "SLC6A8", "Synapto."), ]

# pdf(paste0('figures/heatmaps/all', '_', group, '_', region, '_meanNoStd.pdf'))
pdf(paste0('test.pdf'))
heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=as.dendrogram(hclust(dist((as.matrix(aaa))))), col=pal,
# heatmap.2(as.matrix(t(aaa)), Rowv=FALSE, Colv=FALSE, col=pal, main="Hippocampus",
          scale='row', trace='none', srtCol=90, adjCol = c(0.75, 0.5),
          lhei=c(2, 12), lwid=c(2, 5), key.par=list(cex=0.622), density.info="none",
          key.title='', keysize=0.1, margins=c(5,7)) #+ ggtitle('Caudate')
dev.off()






# quick description of each cluster -----------------------------------------
aa <- aggregate(df_clus[, -1], by=list(df_clus$mc_sub), FUN=median)
bb <- cbind.data.frame(Cluster=rep(unique(aa[, 1]), times=ncol(aa[, -1])), 
                    Marker=rep(colnames(aa[, -1]), each=nrow(aa[, -1])), value=c(as.matrix(aa[, -1])))
bb <- bb[order(bb$Marker, decreasing=FALSE), ]

labels_ <- bb[bb$value>0, ]
labels_ <- labels_ %>% group_by(Cluster) %>% mutate(y=paste(Marker, collapse = "+ "))
labels_ <- labels_[!(duplicated(labels_$Cluster)), ]
labels_[order(labels_$Cluster), -c(2, 3)]



# create venn diagram
library(VennDiagram)
# Writing to file
pdf("figures/others/Triple_Venn_diagram.pdf")
venn.plot <- draw.triple.venn(area1 = length(label_map[['BA9']]),
                              area2 = length(label_map[['Hipp']]),
                              area3 = length(label_map[['DLCau']]),
                              n12 = length(intersect(label_map[['BA9']], label_map[['Hipp']])),
                              n23 = length(intersect(label_map[['Hipp']], label_map[['DLCau']])),
                              n13 = length(intersect(label_map[['BA9']], label_map[['DLCau']])),
                              n123 = length(intersect(intersect(label_map[['BA9']], label_map[['Hipp']]), label_map[['DLCau']])),
                              category = c("BA9", "Hipp", "DLCau"),
                              fill = c("#485696", "#FC5130", "#E0A458"),
                              # lty = "blank",
                              cex = 3.5,
                              cat.cex = 2.5,
                              cat.col = c("#485696", "#FC5130", "#E0A458"))
grid.draw(venn.plot)
dev.off()






library(ggplot2)

cutoff <- read.csv('../raw_data/cutoff_gate.csv')
for (col in colnames(df)){
    min_value <- min(df[col])
    cutoff_ <- cutoff[cutoff[, 'SynTOF.Antibody.Panel']==col, 'arcsinh0.2.cutoff']
    df[df[col] < cutoff_, col] = min_value
}
df['NET'] <- NULL

df_scaled <- scale(df[!aa, ])

# Histogram of feature values
df_ <- do.call(rbind.data.frame, df_scaled)
sampled <- sample((1):1000000, 100000)
df_sub <- df_[sampled, ]
df_plot <- as.data.frame(df_sub) %>% gather()

pdf('figures/histograms/Hipp_cutoff_no0_LowNoSampled100k.pdf') 
ggplot(gather(df_plot), aes(value)) +  
geom_histogram(bins = 100) +  
facet_wrap(~key, scales = 'free_x') 
dev.off() 



