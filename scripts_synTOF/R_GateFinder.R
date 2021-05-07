library(flowCore)
library(GateFinder)
library(ggpubr)
library(dplyr)
library(tidyr)
source('R_utils_preprocessing.R')

# load data -----------------------------------------------------------------------------------
region <- 'Hipp'

# load data
df <- list()
for (group in c('LowNo', 'LBD', 'PHAD')) {
    # get files
    fcs_path <- '../raw_data/max_events/fcs/'
    file_list <- list.files(fcs_path, pattern=paste0(region, '_', group, '*'))
    df[[group]] <- data.frame()
    for (i in seq(length(file_list))) {
        frame <- read.FCS(paste0(fcs_path, file_list[i]))
        df[[group]] <- rbind.data.frame(df[[group]], exprs(frame))
        file_name_splt <- strsplit(file_list[i], '_')[[1]]
    }
    cutoff <- read.csv('../raw_data/cutoff_gate.csv')
    for (col in colnames(df[[group]])){
        min_value <- min(df[[group]][col])
        cutoff_ <- cutoff[cutoff[, 'SynTOF.Antibody.Panel']==col, 'arcsinh0.2.cutoff']
        df[[group]][df[[group]][col] < cutoff_, col] <- min_value
    }
    df[[group]]['NET'] <- NULL
}

# define surface and functional markers
functionalPro <- sort(c('b-Amyloid_X40', 'b-Amyloid_X42', 'p-Tau', 'a-Synuclein_pS129',
                'EAAT1', 'GFAP', 'Casp3_Acti', '3NT', 'LC3B', 'K48-Ubiquitin'))
surfacePro <- sort(colnames(df[[1]])[!colnames(df[[1]]) %in% functionalPro])

df_scaled <- df



pro <- 'ApoE'
n_events <- 100000
dict_no <- 1
thres <- 0.86

targetpop <- as.matrix(df_scaled[[dict_no]][, pro] > thres)
##Subset the markers that should be considered for gating.
prop_markers <- c('CD56', 'CD47', 'SNAP25', 'GATM_1', 'a-Synuclein', 'VGLUT', 'BIN1', 'TMEM230_C20orf30')
#c('CD56', 'SNAP25', 'GATM_1', 'VGLUT', 'BIN1', 'a-Synuclein', 'CD47', 'TMEM230_C20orf30', 'GAD65', 'b-Amyloid_X42')
x <- df_scaled[[dict_no]][, prop_markers]
colnames(x) <- prop_markers
##Run GateFinder.
ans <- GateFinder(x, targetpop, max.iter=2, beta=1, nstart=10, subsample=n_events)

##Make the plots.
# png(paste0('figures/gate/', region, '_', pro, '_group', dict_no, '_nsub', n_events, '_thres', thres ,'.png'), 
          #  width=1000, height=1000, pointsize=24)
png('figures/gate/test.png', width=1000, height=1000, pointsize=30)
plot.GateFinder(x=x, y=ans, c(2, 2), targetpop=targetpop)#, xlim=c(-3, 3), ylim=c(-3, 3))
dev.off()





sum(df_scaled[[dict_no]][1:n_events, pro] > thres)
summary(df_scaled[[dict_no]][, pro])

quantile(df_scaled[[1]][, pro], c(0.999))
quantile(df_scaled[[3]][, pro], c(0.995))







library(flowCore)
data(LPSData)
##Select the target population. In this case cells with those with a pP38 expression (dimension 34) of higher than 3.5.
targetpop <- (exprs(rawdata)[,34]>3.5)
##Subset the markers that should be considered for gating.
x=exprs(rawdata)[,prop.markers]
colnames(x)=marker.names[prop.markers]
##Run GateFinder.
ans=GateFinder(x, targetpop)


quantile(exprs(rawdata)[,34], c(0.995))








# try gate finder ----------------------------------------------------
library(GateFinder)
library(gatepoints)
library(data.table)
library(GateFinder)
library(tidyr)


# load clusters
mc <- as.data.frame(fread(paste0('R_py_exchange/mcResultsDWH_allGroups_maxK40_allLowNoPresynaptic_LowNo_08062020Batch210_105_Adagradlr01_noStd_sess_1.csv')))
mc <- mc %>% separate(sample, c("pp", "region", "group", "batch", "sample"), sep='_')
df_region <- readRDS(paste0('df_pre_save/exp', 0, '.rds'))

# X=readRDS('~/CellMat') #matrix of cells by phenotypic markers (samples cells across files)
prop_markers <- c('CD56', 'CD47', 'GAD65', 'SNAP25', 'GATM_1', 'a-Synuclein', 'VGLUT', 'BIN1', 'TMEM230_C20orf30', 'GAMT_2', 'APP')
X <- df_region[['Hipp']][['LowNo']][5000:10100, prop_markers]
mc_ <- mc[5000:10100, ]

#plot this so you can select points
# temp2=cbind(umap_xy[,1],umap_xy[,2])
# X11()
# plot(umap_xy[,1],umap_xy[,2],cex=.3,pch=16)
# pop=fhs(temp2) #allows you to select points with cursor
# #which cells were in circled population
# goodInds=as.numeric(pop)
# B=rep(0,nrow(X))
# B[goodInds]=1
# B=as.logical(B)

B <- as.matrix(mc_$mc==2)
gf <- GateFinder(X, B, max.iter=2, nstart=10, subsample=1000) #apply GF

png('figures/gate/test2.png', width=1000, height=1000, pointsize=30)
plot.GateFinder(x=X, y=gf, c(2,3), targetpop=B) # plot gatefinder result
dev.off()







colnames(df_region[['Hipp']][['LowNo']][5000:10100, ])

























































library(Morpho)
library(Rtsne)
library(Seurat)

mc <- read.csv(paste0('R_py_exchange/mcResultsDWH_allGroups_05112020Batch215maxK40_', region, '.csv'))
hidden <- read.csv(paste0('R_py_exchange/hidden_', region, 'LowNo_05112020Batch215_105only_sess_1.csv'), nrows=3000)[-1]
hidden <- prcompfast(df[[1]][1:3000, ])$x[, 1:2]

hidden <- uwot::umap(df_scaled[[1]][1:3000, surfacePro])
aa <- slingshot(hidden, as.factor(mc[1:3000, 1]))


png('test.png')
plot(hidden, col=as.factor(mc[1:3000, 1]), pch=16, cex=0.5)
lines(aa, lwd=3, type = 'lineages')
dev.off()




# For reproducibility
RNGversion("3.5.0")
palette(brewer.pal(8, "Dark2"))
library(slingshot)
library(tradeSeq)

data <- df_scaled[[1]][1:3000, surfacePro]
colMin <- apply(data, 2, min)
data <- t(t(data) + abs(colMin))

mc <- read.csv(paste0('R_py_exchange/mcResultsDWH_allGroups_05112020Batch215maxK40_', region, '.csv'))

hidden <- uwot::umap(data)
crv <- slingshot(hidden, as.factor(mc[1:3000, 1]))
counts <- t(data)

pseudotime <- slingPseudotime(crv, na=FALSE)
cellWeights <- slingCurveWeights(crv)
sce <- fitGAM(counts=counts, pseudotime=pseudotime, cellWeights=cellWeights, nknots=600)

assoRes <- associationTest(sce)
head(assoRes)

startRes <- startVsEndTest(sce)
oStart <- order(startRes$waldStat, decreasing = TRUE)
sigGeneStart <- names(sce)[oStart[25]]
plotSmoothers(sce, counts, gene=sigGeneStart)
dev.off()
plotGeneCount(crv, counts, gene = sigGeneStart)
dev.off()

endRes <- diffEndTest(sce)
o <- order(endRes$waldStat, decreasing = TRUE)
sigGene <- names(sce)[o[1]]
plotSmoothers(sce, counts, sigGene)
dev.off()
plotGeneCount(crv, counts, gene = sigGene)
dev.off()

# Discovering genes with different expression patterns
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])
plotSmoothers(sce, counts, gene = rownames(patternRes)[oPat][1])
dev.off()
plotGeneCount(crv, counts, gene = rownames(patternRes)[oPat][1])
dev.off()




yhat <- predictCells(models=sce, gene="TMEM230_C20orf30")
ysmooth <- predictSmooth(models=sce, gene="TMEM230_C20orf30", nPoints=40)
library(clusterExperiment)
nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints=nPointsClus, genes=rownames(counts))#, clusterFunction='clara')

clusterLabels <- primaryCluster(clusPat$rsec)
library(cowplot)

cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
plots <- list()
for (xx in cUniq[1:4]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic()
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c("orange", "darkseagreen3"),
                       breaks = c("0", "1"))
  plots[[as.character(xx)]] <- p
}
plots$ncol <- 2
do.call(plot_grid, plots)






































