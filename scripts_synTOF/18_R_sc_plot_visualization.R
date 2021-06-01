"""
This script visualizes all scatter tsne plots in single cell levels (both colored by
cluster assignments and colored by expression values of each marker).
"""

#### libraries -----------------------------------------------------------------------
library(yarrr)
library(randomcoloR)
library(ggpubr)
library(data.table)
library(tidyr)
library(viridis)
library(dplyr)
library(parallel)

mcsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

get_fdrp <- function(r, n=24, log=TRUE) {
#     zY <- sqrt((n-3)/1.06)*atanh(r)
#   aaP <- p.adjust(2*pnorm(-abs(zY)), method="fdr")
    zY <- r*sqrt((n-2)/((r+1)*(1-r)))
    aaP <- 2*(1-pt(abs(zY), df=(n-2)))
#     aaP[aaP<0] <- 0
    aaP <- p.adjust(aaP, method="fdr")
    if (log==TRUE) {
        aaP <- -log10(aaP)
    }
    aaP
}

#########################################################################################
# plot single cells with color by cluster assignment #
#########################################################################################
#### Importing data ------------------------------------------------------------------
df_plot <- read.csv('R_py_exchange_afterCluster/umap_sc_randomLowNo_allRegions.csv')[, -1]

# for adjustiing column names
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
df_plot$c <- sapply(df_plot$c, function(x) new_label[x])
df_plot$c <- factor(df_plot$c, levels=c('A1', 'A2', 'B1', 'B2', paste0('C', 1:11)), ordered=TRUE)

p <- ggscatter(df_plot, x="x", y="y", color="c", size=0.05*40, palette=
    #  c('#128a6a', '#a909e8', '#09e87c', '#0948e8', 'black', '#e80966', '#124a8a', '#93e809', '#8a5012',
    #    '#e87c09', '#09c3e8', '#737072', '#4409e8', '#e8bf09', '#09e8bf')) + 
    c('#845EC2', '#D65DB1', '#FF6F91', '#FF9671', 
      'black', '#A8BB5C', '#5CA46E', '#20887A', '#1C6873',
      '#2F4858', '#716D6B', '#A5A09E', '#00AAAB', '#FFC75F', '#626E9C')) + 
     theme_classic() + ylab('') + xlab('') +
     theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=14*4), 
          axis.text=element_blank(), axis.title=element_text(size=14*4),
          axis.ticks=element_blank(), panel.border = element_rect(colour="black", fill=NA, size=2)) +
          guides(colour = guide_legend(override.aes = list(size=3*4)))

# ggsave(filename=paste0('figures/singlecell/cluster.pdf'), p)#, width=20, height=15)
ggsave(filename=paste0('figures/singlecell/cluster.png'), p, width=20.5, height=22)#, width=20, height=15)



#########################################################################################
# plot single cells with color by expression levels of each marker #
#########################################################################################

# import raw data and cluster
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08312020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
mc[mc$sample %in% c('HF14-017.fcs', 'HF14-083.fcs', 'HF14-025.fcs'), 'group'] <- 'ODC'
df_region <- readRDS(paste0('df_pre_save/exp', 2, 'mc_5_13.rds'))


for (region in names(df_region)) {
     for (group in names(df_region[[region]])) {
          print(all(mc[mc$pp=='pre' & mc$region==region & mc$group==group, 'sample'] == df_region[[region]][[group]][, 'sample']))
     }
}
regions <- c('BA9', 'DLCau', 'Hipp')
selected_groups <- c('LowNo')

# concatenate all groups and regions of median expression values
df_ <- list()
region_group_index <- matrix(NA, length(regions), length(selected_groups), dimnames=list(regions, selected_groups))
summ <- 0
for (region in regions) {
    df_[[region]] <- do.call(rbind.data.frame, df_region[[region]][which(names(df_region[[region]]) %in% selected_groups)])
    region_group_index[region, ] <- nrow(df_[[region]]) + summ    #####<<<<<< KNOW THAT THERE IS A BUG HERE IF REGION > BA9
    summ <- sum(region_group_index[region, ])
}
df_combined <- do.call(bind_rows, df_)


# load x and y coordinates for each cell (generated from python script)
df_xy <- read.csv('R_py_exchange_afterCluster/umap_sc_randomLowNo_allRegions.csv')[, -1]
df_xy <- df_xy[order(df_xy[, 4]), ]
sel_indices <- df_xy[, 4]

# define region for each cell and reorder
corresponding_region <- sapply(sel_indices, function(x) if (x<=region_group_index['BA9', 'LowNo']) {'BA9'} 
                                                        else if (x<=region_group_index['DLCau', 'LowNo']) {'DLCau'} else {'Hipp'})
sel_df <- df_combined[df_xy[, 4], ]

# define marker types
functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                    'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))#,
surfacePro <- sort(colnames(sel_df)[!colnames(sel_df) %in% c('mc', 'sample', functionalPro)], decreasing=TRUE)

# rescale all markers to a range of 0 to 1
df__ <- sel_df[, surfacePro]
df__ <- apply(df__, 2, function(x) range01(x))

# rename columns a bit
colnames(df__) <- sapply(colnames(df__), function(x) if(x=='DJ-1_PARK7') {'DJ1'} 
                                                   else if(x=='TMEM230_C20orf30') {'TMEM230'} 
                                                   else if(x=='PrP_CD230') {'PrP'}
                                                   else if(x=='GATM_1') {'GATM'}
                                                   else if(x=='GAMT_2') {'GAMT'}
                                                   else if(x=='PARKIN') {'Parkin'}
                                                   else if (x=='a-Synuclein') {'AS'}
                                                   else if (x=='p-Tau') {'PHF-tau'}
                                                   else if (x=='VGLUT') {'vGLUT'}
                                                   else x) 
df_clus <- cbind.data.frame(df_xy[, -4], df__)
new_label <- c('C1', 'C10', 'C3', 'C4', 'B1', 'C5', 'C11', 'A1', 
               'C2', 'C7', 'C9', 'C6', 'B2', 'A2', 'C8')
df_clus$c <- sapply(df_clus$c, function(x) new_label[x])

# pivot longer for plotting
df_plot <- df_clus %>% pivot_longer(!c('x', 'y', 'c'), names_to="marker", values_to="value")

# plot
png('figures/singlecell/expression.png', height=4000, width=4000)
# p <- ggplot(df_plot, aes(x=x, y=y, color=value)) + geom_point() + facet_grid(marker~x)
p <- ggscatter(data=df_plot, x='x', y='y', color='value', facet.by='marker', size=0.5, show.legend.text=FALSE,
               ggtheme=theme_bw()) +
     theme(legend.title=element_text(size=65), strip.text=element_text(size=65), legend.position=c(0.9, 0.1),
          axis.text.y=element_blank(), axis.text.x=element_blank(), axis.title=element_blank(), 
          legend.key.size=unit(1.2*5, "lines"), legend.text=element_text(size=65),
          axis.ticks.y=element_blank(), panel.border=element_rect(colour = "black", fill=NA, size=2)) +
     scale_color_viridis(option="magma", name = "Min-max normalized\nasinh-transformed value")
p
dev.off()