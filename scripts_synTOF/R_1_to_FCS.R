# Libraries -----------------------------------------------------------------------
library(flowCore)
library(Biobase)


# extract groups, batch, and brain regions from row names -------------------------
ori_path <- '../raw_data/max_events/original_files_post_synap/'
file_list <- list.files(ori_path)

splt_str <- sapply(file_list, function(x) strsplit(x, '_|\\.'))
group <- unlist(sapply(splt_str, function(x) x[grepl('LowNo|LBD|PHAD', x)]))
batch <- factor(sapply(splt_str, function(x) substring(x[grepl('BC0', x)], 1, 100L)))
region <- unlist(sapply(splt_str, function(x) tail(x[grepl('BA9|DLCau|Hipp|VMCau', x)], n=1)))
techOrNot <- sapply(splt_str, function(x) any(grepl('TC1', x)))
sample_ID <- sapply(splt_str, function (x) if (any(grepl('TC1', x))) {x[grepl('TC1', x)]
} else {paste(x[which(grepl('LowNo|LBD|PHAD', x))-2], x[which(grepl('LowNo|LBD|PHAD', x))-1], sep='-')})


# Importing -----------------------------------------------------------------------
group_tech <- sapply(seq(length(techOrNot)), function(i) if (techOrNot[i]==TRUE) {'Tech'} else {group[i]})

# cell_counts <- c() #min cell counts excluding vmcaudate=51k
# Exporting to feather files for faster future readability and sample names changed
for (file_no in seq(length(file_list))) {
    print(paste('Working on file no.', file_no))
    if ((group_tech[file_no] %in% c('Tech')) | (region[file_no] %in% c('VMCau'))) {
        next
    }
    txt <- read.csv(paste0(ori_path, file_list[file_no]), sep='\t', 
                    header=FALSE, stringsAsFactors=FALSE)
    pro_names <- sapply(txt[2, ], function(x) strsplit(x, '_'))
    pro_names <- sapply(pro_names, function(x) if (length(x) == 1) {x[1]} 
                        else if (length(x) == 4) {paste(x[2], x[3], sep='_')}
                        else {x[2]})
    colnames(txt) <- pro_names
    txt <- data.matrix(txt[3:nrow(txt), 2:ncol(txt)])
    metadata <- data.frame(name=colnames(txt), desc=colnames(txt))
    metadata$minRange <- apply(txt, 2, min)
    metadata$maxRange <- apply(txt, 2, max)

    fcs <- new("flowFrame", exprs=txt, parameters=AnnotatedDataFrame(metadata))

    # export to feather
    sample_name <- paste(region[file_no], group_tech[file_no], batch[file_no],
                   sample_ID[file_no], sep='_')
    write.FCS(x=fcs, filename=paste0("../raw_data/max_events/fcs_post_synap/", sample_name, '.fcs'))
}





# For GFAP-EAAT1- -------------------------
ori_path <- '../raw_data/max_events/original_files_GFAPnegEAAT1neg/'
file_list <- list.files(ori_path)

splt_str <- sapply(file_list, function(x) strsplit(x, '_|\\.'))
group <- unlist(sapply(splt_str, function(x) x[grepl('LowNo|LBD|PHAD', x)]))
batch <- factor(sapply(splt_str, function(x) substring(x[grepl('BC0', x)], 1, 100L)))
region <- unlist(sapply(splt_str, function(x) tail(x[grepl('BA9|DLCau|Hipp|VMCau', x)], n=1)))
techOrNot <- sapply(splt_str, function(x) any(grepl('TC1', x)))
sample_ID <- sapply(splt_str, function (x) if (any(grepl('TC1', x))) {x[grepl('TC1', x)]
} else {paste(x[which(grepl('LowNo|LBD|PHAD', x))-2], x[which(grepl('LowNo|LBD|PHAD', x))-1], sep='-')})


# Importing -----------------------------------------------------------------------
group_tech <- sapply(seq(length(techOrNot)), function(i) if (techOrNot[i]==TRUE) {'Tech'} else {group[i]})

# cell_counts <- c() #min cell counts excluding vmcaudate=51k
# Exporting to feather files for faster future readability and sample names changed
for (file_no in seq(length(file_list))) {
    print(paste('Working on file no.', file_no))
    if ((group_tech[file_no] %in% c('Tech')) | (region[file_no] %in% c('VMCau'))) {
        next
    }
    txt <- read.csv(paste0(ori_path, file_list[file_no]), sep='\t', 
                    header=FALSE, stringsAsFactors=FALSE)
    pro_names <- sapply(txt[2, ], function(x) strsplit(x, '_'))
    pro_names <- sapply(pro_names, function(x) if (length(x) == 1) {x[1]} 
                        else if (length(x) == 4) {paste(x[2], x[3], sep='_')}
                        else {x[2]})
    colnames(txt) <- pro_names
    txt <- data.matrix(txt[3:nrow(txt), 2:ncol(txt)])
    metadata <- data.frame(name=colnames(txt), desc=colnames(txt))
    metadata$minRange <- apply(txt, 2, min)
    metadata$maxRange <- apply(txt, 2, max)

    fcs <- new("flowFrame", exprs=txt, parameters=AnnotatedDataFrame(metadata))

    # export to feather
    sample_name <- paste(region[file_no], group_tech[file_no], batch[file_no],
                   sample_ID[file_no], sep='_')
    write.FCS(x=fcs, filename=paste0("../raw_data/max_events/fcs_GFAPnegEAAT1neg/", sample_name, '.fcs'))
}


# # Generating fcs_info.csv and sample_info.csv for MetaCyto pacakge -----------------
# fcs_path <- '../raw_data/max_events/fcs/'
# file_list <- list.files(fcs_path)
# splt_str <- sapply(file_list, function(x) strsplit(x, '_|\\.'))
# region <- sapply(splt_str, function(x) x[1])
# group <- sapply(splt_str, function(x) x[2])
# batch <- sapply(splt_str, function(x) x[3])
# sample_ID <- sapply(splt_str, function(x) x[4])

# demo <- read.csv('../raw_data/demographics.csv')
# demo$sample_ID <- gsub(' ', '', demo$sample_ID)

# fcs_info_ <- cbind.data.frame(file_list, sample_ID, region, group, batch)
# fcs_samp_info <- merge(fcs_info_, demo, by='sample_ID')
# fcs_samp_info$file_list <- paste0('../raw_data/max_events/fcs/', fcs_samp_info$file_list)

# fcs_info <- fcs_samp_info[, c(2, 1)]
# colnames(fcs_info) <- c('fcs_files', 'study_id')
# sample_info <- fcs_samp_info[, -1]
# colnames(sample_info)[1] <- 'fcs_files'


# write.csv(fcs_info, paste0(fcs_path, 'fcs_info.csv'), row.names=FALSE)
# write.csv(sample_info, paste0(fcs_path, 'sample_info.csv'), row.names=FALSE)


