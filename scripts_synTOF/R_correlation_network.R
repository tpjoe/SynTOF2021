#### libraries -----------------------------------------------------------------------
source('R_utils_plots.R')
library(SDMTools) # for legend
library(Rtsne)
library(yarrr)
library(uwot)



#### Importing data ------------------------------------------------------------------
for (region in c('BA9', 'DLCau', 'Hipp')) {
    df_ <- read.csv(paste0('R_py_exchange/df_meanAllMarkers_', region, '_noStd.csv'))[, -1]
    if (region != 'BA9') {
        #check that samples are in the same order
        print(all(df[, 2] == df_[, 2]))
        colnames(df_) <- paste0(region, '_', colnames(df_))
        df <- cbind.data.frame(df, df_[, 3:ncol(df_)])
    } else {
        colnames(df_)[3:ncol(df_)] <- paste0(region, '_', colnames(df_[, 3:ncol(df_)]))
        df <- df_
    }
    print(dim(df_))
}
pair <- c('LowNo', 'PHAD')

# fill na and get rid of all-0 columns
for(i in 1:ncol(df)){
  df[is.na(df[,i]), i] <- mean(df[,i], na.rm = TRUE)
}

X <- df[df$group %in% pair, !colnames(df) %in% c('group', 'sample')]
X_all <- df[, !colnames(df) %in% c('group', 'sample')]
X <- X[, !(seq(ncol(X_all)) %in% which(colSums(X_all==0) == nrow(X_all)))]
X_all <- X_all[, !(seq(ncol(X_all)) %in% which(colSums(X_all==0) == nrow(X_all)))]
y <- df[df$group %in% pair, 'group']
y <- sapply(y, function(x) if (x==pair[1]) {0} else {1})

# calculate correlations of LowNo
corrX <- cor(X_all, method='spearman')
corrX[is.na(corrX)] <- 0

# corrY <- cor(X, y, method='spearman')
# corrY[is.na(corrY)] <- 0
# zY <- sqrt((nrow(df)-3)/1.06)*atanh(corrY)
# corrYP <- 2*pnorm(-abs(zY))
# nodeSize <- -log10(corrYP)

# layoutSample <- Rtsne(X_all, perplexity=5, check_duplicates=FALSE, pca=TRUE, theta=0)$Y
# layoutSample <- cbind.data.frame(layoutSample, df$group)
# names(layoutSample) <- c('x', 'y', 'group')
# ggplot(layoutSample, aes(x=x, y=y, color=group)) + geom_point(size=5)
# ggsave('test.pdf')

# Draw networks by r ------------------------------------------------------------------
set.seed(1)
layout <- Rtsne(asinh(corrX), perplexity=40, check_duplicates=FALSE, pca=FALSE, theta=0)$Y

network_graph <- plot_network(X=X, y=y, layout=layout, edgeby='p',  gradient='r')#,
                              # labels=TRUE, reduced=25)
maxlogP <- round((max(V(network_graph)$size)-1)/3, digits=2)
size <- c(round(maxlogP/4), round(maxlogP/3), round(maxlogP/2))

# Exporting graph
png(filename=paste0('figures/networks/allMeanAllMarkers_', pair[1], pair[2], '_noStd_mean.png'), width=4000, height=4000)
plot(network_graph, vertex.label=V(network_graph)$name)
###---- Legend and circles
pnts = cbind(x=c(0.80,0.85,0.80,0.85), y=c(-.97,-0.97,-0.70,-0.70))
legend.gradient(pnts, cols = colorRampPalette(c("blue", "white", "red"))(100),
                limits=pair, title="", cex=9)
legends <- legend(-1.05, 1.05, legend=size, cex=10, border='red', y.intersp=1.3,
            title=expression(paste('-log' ["10"], '(', bolditalic('P'), " value)")), bty='n')
xx <- (legends$text$x + legends$rect$left)/2
yy <- legends$text$y
text(xx[1] + 0.069, yy[3] + (yy[3]-yy[2]), paste0('max = ', round(maxlogP)), cex=9)
symbols(xx, yy, circles=size/60, inches=FALSE, add=TRUE, bg='light grey', fg='black')
title(main=paste0(pair[1], pair[2], '_mean'), cex.main=10)
dev.off()


### By regions -------------------------------------------------------------------
# Draw networks
set.seed(1)
layout <- Rtsne(asinh(corrX), perplexity=40, check_duplicates=FALSE, pca=FALSE, theta=0)$Y
# layout <- umap(corrX, min_dist=0.1)

# Getting stimulant name for each feature column
region_ <- table(sapply(colnames(X), function(x) strsplit(x, '_')[[1]][1]))
region_ <- c('BA9'=rep(1, region_[1]), 'DLCau'=rep(2, region_[2]), 'Hipp'=rep(3, region_[3]))
region_palette <- c('BA9'='#605856', 'DLCau'='#1C6E8C', 'Hipp'='#FF4A1C')

# Getting graph colored by stimulants
network_graph <- plot_network(X=X, y=y, layout=layout, edgeby='p', 
                             makeClus=TRUE, clusters=region_, cluster_palette=region_palette)

# Exporting graph
png(filename=paste0('figures/networks/allMeanAllMarkers_std_region.png'), width=4000, height=4000)
plot(network_graph, vertex.label=V(network_graph)$name, vertex.size=2.5)
###---- Legend
legends <- names(region_palette)
a <- legend(-1.1, 1.1, legend=legends,
            bty='n', cex=10, pt.cex=0.1) #2.5 for tsne 4 for umap
xx <- (a$text$x + a$rect$left)/2
yy <- a$text$y
symbols(xx, yy, circles=rep(4.5, length(region_palette))/200, inches=FALSE,
        add=TRUE, bg=as.character(region_palette),
        fg='black')
dev.off()


### By exprs types -------------------------------------------------------------------
# Draw networks
library(circlize)
set.seed(1)
layout <- Rtsne(asinh(corrX), perplexity=40, check_duplicates=FALSE, pca=FALSE, theta=0)$Y
# Getting stimulant name for each feature column
region_ <- table(sapply(colnames(X), function(x) strsplit(x, '_')[[1]][3]))
region_palette <- rand_color(35)
names(region_palette) <- seq(1:35)

# Getting graph colored by stimulants
network_graph <- plot_network(X=X, y=y, layout=layout, edgeby='p', label=FALSE,
                             makeClus=TRUE, clusters=region_, cluster_palette=region_palette)

# Exporting graph
png(filename=paste0('figures/networks/allMeanAllMarkers_exprs.png'), width=4000, height=4000)
plot(network_graph, vertex.label=V(network_graph)$name, vertex.size=2.5)
###---- Legend
legends <- names(region_palette)
a <- legend(-1.1, 1.1, legend=legends,
            bty='n', cex=10, pt.cex=0.1) #2.5 for tsne 4 for umap
xx <- (a$text$x + a$rect$left)/2
yy <- a$text$y
symbols(xx, yy, circles=rep(4.5, length(region_palette))/200, inches=FALSE,
        add=TRUE, bg=as.character(region_palette),
        fg='black')
dev.off()