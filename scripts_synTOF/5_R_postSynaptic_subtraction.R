"""
This script performs post-processing of the clustered data by subtracting signals we know should not exist in the
post-synaptic out from both post- and per- synaptic expression values.
"""

# load corrected pre-synaptic data ----------------------------------------------------------------------
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

# define list of markers we know should not existing post-synaptic events
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
        # check that sample order aligns with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        # laod fcs files
        for (i in seq(length(file_list))) {
            frame <- read.FCS(paste0(fcs_path, file_list[i]))
            df_post[[group]] <- rbind.data.frame(df_post[[group]], exprs(frame))
            file_name_splt <- strsplit(file_list[i], '_')[[1]]
        }
        # remove low quality marker
        df_post[[group]]['NET'] <- NULL
    }

    # define type of markers
    functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
    surfacePro <- sort(colnames(df_post[[group]])[!colnames(df_post[[group]]) %in% functionalPro], decreasing=TRUE)

    # attaching cluster_id and sample_id
    for (group in names(df_post)) {
        df_post[[group]] <- cbind.data.frame(mc[(mc$pp=='post') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_post[[group]])
    }

    df_post_ref <- list()
    # get mean expression from the post synaptic events
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
        # check that sample order aligns with mc
        df_file_name <- as.data.frame(file_list) %>% separate(file_list, c("region", "group", "batch", "sample"), sep='_')
        if (any(unique(mc[(mc$region==region) & (mc$group==group), 'sample']) != df_file_name[, 'sample'])) {
            print('shitttttttttttttttttttttttttttttttt')
        }
        for (i in seq(length(file_list))) {
            frame <- exprs(read.FCS(paste0(fcs_path, file_list[i])))
            """
            ## subtraction begins here ##
            """
            sample <- strsplit(file_list[i], '_')[[1]][4]
            mc_sample <- mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group) & (mc$sample==sample), ]
            post_strsplit <- strsplit(colnames(df_post_ref[[group]][df_post_ref[[group]][, 'sample']==sample, ]), '_')
            post_cl <- sapply(post_strsplit, function(x) x[length(x)])
            post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
            for (marker in adjust_markers) {
                uni_cl <- sort(unique(post_cl))
                # for (cl in uni_cl[!uni_cl %in% c('group', 'sample')]) {
                for (cl in new_label[c(5, 13)]) { # cluster 5 and 13 are basically B1 and B2, and are essentially only ones postsynaptic events classified to
                    # get post-synaptic subtraction values (ref) of that marker for each cluster in LowNo samples
                    post_strsplit <- strsplit(colnames(df_post_ref[['LowNo']][df_post_ref[['LowNo']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_LowNo <- tryCatch(df_post_ref[['LowNo']][, (post_marker==marker) & (post_cl==cl)],
                                          error=function(x) {0})
                    if (length(ref_LowNo)==0) {ref_LowNo<-0}

                    # get post-synaptic subtraction values (ref) of that marker for each cluster in LBD samples
                    post_strsplit <- strsplit(colnames(df_post_ref[['LBD']][df_post_ref[['LBD']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_LBD <- tryCatch(df_post_ref[['LBD']][, (post_marker==marker) & (post_cl==cl)],
                                        error=function(x) {0})
                    if (length(ref_LBD)==0) {ref_LBD<-0}

                    # get post-synaptic subtraction values (ref) of that marker for each cluster in PHAD samples
                    post_strsplit <- strsplit(colnames(df_post_ref[['PHAD']][df_post_ref[['PHAD']][, 'sample']==sample, ]), '_')
                    post_cl <- sapply(post_strsplit, function(x) x[length(x)])
                    post_marker <- sapply(post_strsplit, function(x) if (length(x)==4) {paste(c(x[1], x[2]), collapse='_')} else {x[1]})
                    ref_PHAD <- tryCatch(df_post_ref[['PHAD']][, (post_marker==marker) & (post_cl==cl)],
                                          error=function(x) {0})
                    if (length(ref_PHAD)==0) {ref_PHAD<-0}

                    # get mean values in each group
                    mean_ref <- c()
                    for (r in c(list(ref_LowNo), list(ref_LBD), list(ref_PHAD))) {
                        mean_ref <- c(mean_ref, mean(r, na.rm=TRUE))
                    }
                    # get mean in all groups
                    ref <- mean(mean_ref, na.rm=TRUE)
                    # perform subtraction
                    if (!is.na(ref)) {
                        frame[mc_sample$mc==cl, marker] <- frame[mc_sample$mc==cl, marker] - ref
                    }
                }
            }
            df_pre[[group]] <- rbind.data.frame(df_pre[[group]], frame)
        }
        # remove low quality marker
        df_pre[[group]]['NET'] <- NULL
    }
    # attaching cluster_id and sample_id
    for (group in names(df_pre)) {
        df_pre[[group]] <- cbind.data.frame(mc[(mc$pp=='pre') & (mc$region==region) & (mc$group==group), c('mc', 'sample')], df_pre[[group]])
    }
    # done for 1 region
    df_region[[region]] <- df_pre
}

# remove samples that are found out later to not be LowNo
for (region in c('BA9', 'DLCau', 'Hipp')){
    print(dim(df_region[[region]][['LowNo']]))
    df_LowNo <- df_region[[region]][['LowNo']]
    df_region[[region]][['ODC']] <- df_LowNo[df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    df_region[[region]][['LowNo']] <- df_LowNo[!df_LowNo$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), ]
    print(dim(df_region[[region]][['LowNo']]))
    print(dim(df_region[[region]][['ODC']]))
}

# save final results
saveRDS(df_region, 'df_pre_save/exp2mc_5_13.rds')