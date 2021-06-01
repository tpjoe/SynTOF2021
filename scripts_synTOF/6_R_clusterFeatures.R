"""
This script generate the final feature matrix by loading all previous data and reformat to a proper table format
"""

# Create freq. and exprs matrix ------------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)


for (exp_no in c(2)){ #c(0, seq(6))) {
    # grouping data by IDs
    # exp_no <- 0
    df_region <- readRDS(paste0('df_pre_save/exp', exp_no, 'mc_5_13.rds'))
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
        write.csv(df_freq_region, paste0('R_py_exchange/df_freq_', region, '_noStd.csv'))
        write.csv(df_mean_region, paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd_exp', exp_no, 'mc_5_13.csv'))
    }
}




# load only post-synaptic data ----------------------------------------------------------------------
# load data
library(dplyr)
library(tidyr)
library(flowCore)
library(data.table)


# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')

# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
mc$mc <- sapply(mc$mc, function(x) new_label[x])

# for reordering colum names later
desired_order <- sort(new_label)[c(1:5, 8:15, 6:7)]

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

    write.csv(df_freq, paste0('R_py_exchange/df_freqPost_', region, '_noStd.csv'))
    write.csv(df_mean_region, paste0('R_py_exchange/df_meanPostMarkers_', region, '_noStd.csv'))
}

