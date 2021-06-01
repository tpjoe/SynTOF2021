"""
This script plots out the radar plots used in Figure 1 of the manuscript.
"""


#import libraries
library(tidyverse)
library(ggplot2)
library(scales)


"""
This one is for positive %event plot
"""
# Tom's avg pre ----------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggradar)
suppressPackageStartupMessages(library(dplyr))
library(scales)

normalize <- function(x)
{
    # return((x- min(x)) /(max(x)-min(x)))
    return(x/max(x))
}

# create new coord : inherit coord_polar
coord_radar <- function(theta='x', start=0, direction=1){
    # input parameter sanity check
    match.arg(theta, c('x','y'))

    ggproto(
      NULL, CoordPolar, 
      theta=theta, r=ifelse(theta=='x','y','x'),
      start=start, direction=sign(direction),
      is_linear=function() TRUE)
}

# tom_df <- read.csv('../raw_data/TomsPlot.csv', check.names=FALSE)
tom_df <- read.csv('../raw_data/TomsPlot_new02072021.csv', check.names=FALSE)
marker_type <- tom_df[, 1]
tom_df <- tom_df[, -c(1, 6, 7)]
tom_df[is.na(tom_df)] <- 0
# tom_df[, 2:4] <- t(apply(tom_df[, 2:4], 1, normalize))
tom_df[is.na(tom_df)] <- 0

interested_region <- c('BA9', 'Hippo', 'Caudate')
plot_df <- list()


for (region in interested_region) {
    plot_df[[region]] <- as.data.frame(tom_df[, c(1, grep(region, colnames(tom_df)))])
    # plot_df[[region]][is.na(plot_df[[region]])] <- 0
    colnames(plot_df[[region]]) <- c('Marker', 'Value', 'SD')
    plot_df[[region]][, 'Marker'] <- factor(plot_df[[region]][, 'Marker'], levels=unique(plot_df[[region]][, 'Marker']), ordered=TRUE)
    if (region=='BA9') {
        plot_df[[region]][, 'Region'] <- rep('BA9', 34)
    } else if (region=='Caudate') {
        plot_df[[region]][, 'Region'] <- rep('DLCau', 34)
    } else {
        plot_df[[region]][, 'Region'] <- rep('Hipp', 34)
    }
}

plot_df[['50line']] <- plot_df[['Hippo']]
plot_df[['50line']][, 'Value'] <- 50

plot_df[['10line']] <- plot_df[['Hippo']]
plot_df[['10line']][, 'Value'] <- 10

plot_df[['25line']] <- plot_df[['Hippo']]
plot_df[['25line']][, 'Value'] <- 25


axis_color <- sapply(marker_type, function(x) if(x=='Neuron type') {'#727273'} else if (x=='AD related') {'#d95743'} else if 
                                                (x=='PD related') {'#3251a1'} else if (x=='Energy') {'#579660'} else {'#6f4a94'})

# start plotting
pdf('figures/radar/human_pre.pdf')
ggplot(aes(x=Marker, y=Value, color=Region, group=1), data=plot_df[['BA9']]) +
geom_polygon(data=plot_df[['BA9']], lwd=0.5, fill=NA) +
geom_polygon(data=plot_df[['Caudate']], lwd=0.5, fill=NA) +
geom_polygon(data=plot_df[['Hippo']], lwd=0.5, fill=NA) + 
geom_polygon(data=plot_df[['50line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 
geom_polygon(data=plot_df[['10line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 
geom_polygon(data=plot_df[['25line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 

geom_point(data=plot_df[['BA9']], size=1.5) +
geom_point(data=plot_df[['Caudate']], size=1.5) +
geom_point(data=plot_df[['Hippo']], size=1.5) +

geom_errorbar(aes(xend=Marker, ymin=Value-SD, ymax=Value+SD), data=plot_df[['BA9']], width=0.4, color=adjustcolor('#513B56', alpha.f=0.6)) +
geom_errorbar(aes(xend=Marker, ymin=Value-SD, ymax=Value+SD), data=plot_df[['Caudate']], width=0.4, color=adjustcolor('#348AA7', alpha.f=0.6)) +
geom_errorbar(aes(xend=Marker, ymin=Value-SD, ymax=Value+SD), data=plot_df[['Hippo']], width=0.4, color=adjustcolor('#eda413', alpha.f=0.6)) +

theme_classic() +
scale_color_manual(values=c("#513B56", "#348AA7", "#eda413")) + 
# scale_fill_manual(values=c(adjustcolor("#513B56", alpha.f=0.1), adjustcolor("#348AA7", alpha.f=0.1), adjustcolor("#eda413", alpha.f=0.1))) +
annotate(geom = "text", x=17.9, y =10+1, label = '10%', hjust = 0, vjust = 1, size = 4) +
annotate(geom = "text", x=17.8, y =25+1, label = '25%', hjust = 0, vjust = 1, size = 4) +
annotate(geom = "text", x=17.77, y =50+1, label = '50%', hjust = 0, vjust = 1, size = 4) +

labs(x = "", y = "Positive Events (%)", legend='Title') + scale_y_continuous(trans = log2_trans(),
    breaks = trans_breaks("log2", function(x) 2^x),
    labels = trans_format("log2", math_format(2^.x))) +
theme(legend.position='top', panel.grid.major.y=element_line(size=.1, color="gray"), panel.grid.major.x=element_line(size=.1, color="gray"), #panel.grid.major.x=element_blank(),
      axis.text.x=element_text(angle=0, colour=axis_color, size=9, face='bold'), legend.text=element_text(size=16), legend.title=element_blank(), 
      axis.text.y=element_text(size=16), axis.title.y=element_text(size=16), plot.margin = margin(0, 0, 0, 0, "cm")) +
coord_radar()
dev.off()





"""
This one is for positive pre/post ratios
"""
# Tom's avg pre/post new ----------------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(ggradar)
suppressPackageStartupMessages(library(dplyr))
library(scales)
library(grDevices)


# create new coord : inherit coord_polar
coord_radar <- function(theta='x', start=0, direction=1){
    # input parameter sanity check
    match.arg(theta, c('x','y'))

    ggproto(
      NULL, CoordPolar, 
      theta=theta, r=ifelse(theta=='x','y','x'),
      start=start, direction=sign(direction),
      is_linear=function() TRUE)
}


# import new data
tom_df_pre <- read.csv('../raw_data/TomsPlot_new02072021.csv', check.names=FALSE)
marker_type <- tom_df_pre[, 1]
tom_df_pre <- tom_df_pre[, -c(1, 6, 7)]

# import post data
tom_df <- read.csv('../raw_data/TomsPrePostControls.csv', check.names=FALSE)
tom_df <- tom_df[!is.na(tom_df[, 'Average PrePost']), ]
tom_df[, 'Marker'] <- rep(unique(tom_df[, 'Marker'])[unique(tom_df[, 'Marker'])!=""], each=3)
# reorder marker
tom_ori <- read.csv('../raw_data/TomsPlot.csv', check.names=FALSE)
marker_order <- factor(tom_ori[, 2], level=unique(tom_ori[, 2]), ordered=TRUE)
tom_df[, 'Marker'] <- factor(tom_df[, 'Marker'], levels=levels(marker_order), ordered=TRUE)
tom_df <- tom_df[order(tom_df[, 'Marker']), ]

# calculate new pre/post values
tom_df_pre <- tom_df_pre[tom_df_pre[[1]] %in% unique(tom_df$Marker), ]
as.character(tom_df_pre[[1]]) == as.character(tom_df[tom_df$Region==region, 'Marker']) # check the order matches
# start calculation
for (region in unique(tom_df$Region)) {
    if (region=='Caud') {region_<-'Caudate'} else {region_<-region}
    tom_df[tom_df$Region==region, 'Average PrePost'] <- tom_df_pre[region_]/tom_df[tom_df$Region==region, 'Average Post']
    tom_df[tom_df$Region==region, 'SEM PrePost'] <- 
                    tom_df_pre[paste0('sd_', region_)]/tom_df_pre[region_]*tom_df[tom_df$Region==region, 'Average PrePost']
}


# select relevant columns
marker_type <- tom_ori[, 1]
tom_df <- tom_df[, c('Marker', 'Region', 'Average PrePost', 'SEM PrePost')]
tom_df['Region'] <- gsub('Caud', 'DLCau', tom_df$Region)
tom_df['Region'] <- gsub('Hippo', 'Hipp', tom_df$Region)
# generate plot_df
plot_df <- list()
for (region in unique(tom_df[, 'Region'])) {
    plot_df[[region]] <- as.data.frame(tom_df[tom_df[, 'Region']==region, c('Marker', 'Region', 'Average PrePost', 'SEM PrePost')])
    colnames(plot_df[[region]]) <- c('Marker', 'Region', 'AveragePrePost', 'SEMPrePost')
}

# adjust marker type selection
all(as.character(unique(plot_df[['BA9']])[, 1]) == as.character(tom_ori[, 2][tom_ori[, 2] %in% unique(plot_df[['BA9']])[, 1]]))
marker_type <-marker_type[tom_ori[, 2] %in% unique(plot_df[['BA9']])[, 1]]


# contant lines
plot_df[['50line']] <- plot_df[['Hipp']]
plot_df[['50line']][, 'AveragePrePost'] <- 50

plot_df[['5line']] <- plot_df[['Hipp']]
plot_df[['5line']][, 'AveragePrePost'] <- 5

plot_df[['25line']] <- plot_df[['Hipp']]
plot_df[['25line']][, 'AveragePrePost'] <- 25

# define colors
axis_color <- sapply(marker_type, function(x) if(x=='Neuron type') {'#727273'} else if (x=='AD related') {'#d95743'} else if 
                                                (x=='PD related') {'#3251a1'} else if (x=='Energy') {'#579660'} else {'#6f4a94'})


# start plotting
# pdf('figures/radar/human_prepostRatio.pdf')

p <- ggplot(aes(x=Marker, y=AveragePrePost, color=Region, group=1), data=plot_df[['BA9']]) +
geom_polygon(data=plot_df[['BA9']], lwd=0.5, fill=NA) +
geom_polygon(data=plot_df[['DLCau']], lwd=0.5, fill=NA) +
geom_polygon(data=plot_df[['Hipp']], lwd=0.5, fill=NA) + 
geom_polygon(data=plot_df[['50line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 
geom_polygon(data=plot_df[['5line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 
geom_polygon(data=plot_df[['25line']], color='#747475', fill=NA, lwd=0.5, lty=2) + 

geom_point(data=plot_df[['BA9']], size=1.5) +
geom_point(data=plot_df[['DLCau']], size=1.5) +
geom_point(data=plot_df[['Hipp']], size=1.5) +

geom_errorbar(aes(xend=Marker, ymin=AveragePrePost-SEMPrePost, ymax=AveragePrePost+SEMPrePost), data=plot_df[['BA9']], width=0.4, color=adjustcolor('#513B56', alpha.f=0.6)) +
geom_errorbar(aes(xend=Marker, ymin=AveragePrePost-SEMPrePost, ymax=AveragePrePost+SEMPrePost), data=plot_df[['DLCau']], width=0.4, color=adjustcolor('#348AA7', alpha.f=0.6)) +
geom_errorbar(aes(xend=Marker, ymin=AveragePrePost-SEMPrePost, ymax=AveragePrePost+SEMPrePost), data=plot_df[['Hipp']], width=0.4, color=adjustcolor('#eda413', alpha.f=0.6)) +

theme_classic() +
scale_color_manual(values=c("#513B56", "#348AA7", "#eda413")) + 
# scale_fill_manual(values=c(adjustcolor("#513B56", alpha.f=0.1), adjustcolor("#348AA7", alpha.f=0.1), adjustcolor("#eda413", alpha.f=0.1))) +
annotate(geom="text", x=15.8, y =5+1, label='5x', hjust=0, vjust=1, size=4) +
annotate(geom="text", x=15.8, y =25+1, label='25x', hjust=0, vjust=1, size=4) +
annotate(geom="text", x=15.77, y =50+1, label='50x', hjust=0, vjust=1, size=4) +

labs(x = "", y = "Pre/Post Ratio", legend='Title') + scale_y_continuous(trans = log2_trans(),
    breaks = trans_breaks("log2", function(x) 2^x),
    labels = trans_format("log2", math_format(2^.x))) +
theme(legend.position='top', panel.grid.major.y=element_line(size=.1, color="gray"), panel.grid.major.x=element_line(size=.1, color="gray"), #panel.grid.major.x=element_blank(),
      axis.text.x=element_text(angle=0, colour=axis_color, size=9, face='bold'), legend.text=element_text(size=16), legend.title=element_blank(), 
      axis.text.y=element_text(size=16), axis.title.y=element_text(size=16), plot.margin = margin(0, 0, 0, 0, "cm")) +
coord_radar()

ggsave('figures/radar/human_prepostRatio.pdf', p)





