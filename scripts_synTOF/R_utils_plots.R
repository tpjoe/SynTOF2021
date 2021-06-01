# import library
library(Rtsne)
library(igraph)
library(reshape)
library(RColorBrewer)
library(ROCit) 



plot_auc <- function (y_preds, ys, xAUCat=c(0.3), yAUCat=c(-24.5)) {
    ROC_list <- list()
    for (i in 1:ncol(y_preds)) {
        y_pred <- c(y_preds[, i][!is.na(y_preds[, i])])
        y <- ys[, i][!is.na(ys[, i])]
        ROC_list[[i]] <- rocit(score=y_pred, class=y) 
    }
    svg('../dumpster/dumppic.svg')
    dev.control(displaylist="enable")
    par(mar=c(5.1, 5.8, 4.1, 2.1))
    plot(ROC_list[[1]]$TPR ~ ROC_list[[1]]$FPR, col=1, lwd=2, type="l", cex.axis=1.5, cex.lab=1.75, las=1,
    xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)", bty="n")
    abline(a=0, b=1, col='grey', lwd=2, lty=2) 
    mtext(bquote(bold(AUC) == .(round(ROC_list[[1]]$AUC, 2))), side=1, line=-14.5, cex=1.5)
    # Drawing additional ROC curve and typing out AUC values at specified x
    for (i in 2:length(ROC_list)) {
        lines(ROC_list[[i]]$TPR ~ ROC_list[[i]]$FPR, col=2, lwd=2, type="l", 
        xlab="1-Specificity (FPR)", ylab="Sensitivity (TPR)")
        mtext(at=xAUCat[i-1], bquote(bold(AUC) == .(round(ROC_list[[i]]$AUC, 2))), side=1, line=yAUCat[i-1], cex=1.5, col=2)
    }
    p <- recordPlot()
    dev.off()
    return(p)
}
 


plot_bar_lines <- function (barValue, xAxisLabels, lineValue=NA, constantLine=NA, yLabelBar, yLabelLine,
                            lineLim=c(0, 1)) {
    # This function plot the bar graphs (with line as second axis if needed).
    #
    # barValue  : The vector value for the bar graph
    # xAxisLabels   : x labels
    # lineValue : The vector value for the line graph (secondary axis)
    # constantLine: Add constant horizontal line to the first axis (bar graph) with
    #               the list containing name h = the vector of constant line at height h
    #               , color = the color of the corresponding line, and lineLabel = labels
    # yLabelBar: The y axis label
    # yLabelLine: The secondary y axis label
    # lineLim   : The secondary y axis limit
    # The function return the figure object

    svg('../dumpster/dumppic.svg')
    dev.control(displaylist="enable")
    par(mgp = c(0, 1, 0), mar=c(4.1, 4.1, 4.1, 4.1), xpd=FALSE) #from c(5.1, 4.1, 4.1, 2.1)
    barplot(barValue, las=2, names.arg=xAxisLabels)
    mp <- barplot(barValue, las=2, cex.axis=1.35, col='light grey')
    if (any(!is.na(constantLine))){
        for (i in 1:length(constantLine)) {
            abline(h=constantLine$h[i], col=constantLine$color[i], lwd=3, lty=2)
        }
        legend('topright', bty="n", legend=constantLine$lineLabel, col=constantLine$color, 
            lty=c(2,2), lwd=3, cex=1.25) #xpd=TRUE
    }
    mtext(yLabelBar, side=2, line=2.4, cex=1.5)
    axis(1, mp, labels=xAxisLabels, mgp=c(0, -1, 0), las=2, lty=0, lwd=0, hadj=0, cex.lab=1.25, cex.axis=1.5)
    if (any(!is.na(lineValue))) {
        par(new=TRUE, mgp = c(0, 1, 0))
        plot(mp, lineValue, axes=FALSE, xlim=c(min(mp)-0.5, max(mp)+0.5), ylim=lineLim,
        bty="n", xlab = "", ylab = "", type='b', col='red', lwd=3)
        axis(4, at= seq(0, 1 ,0.2), col='red', las=1, cex.axis=1.35)
        mtext(yLabelLine, side=4, line=3, cex=1.5)
    }
    p <- recordPlot()
    dev.off()
    return(p)
}



plot_box <- function(values, xAxisLabels=NA, yLabel, spaceBetween=NA, ylim=c(0, 1), jitter=TRUE) {
    svg('../dumpster/dumppic.svg')
    dev.control(displaylist="enable")
    par(mar=c(4.1, 4.3, 7.1, 1.1))
    ats <- sapply(1:ncol(values), function(x) if (all(is.na(values[, x]))) {NA} else {x})
    ats <- ats[!is.na(ats)]
    if (any(is.na(xAxisLabels))) {xAxisLabels <- 1:ncol(values)}
    boxplot(values, col="white", frame=F, las=2, ylab=yLabel, xaxt='n', ylim=ylim, outline=FALSE, cex.lab=1.5, cex.axis=1.5, boxwex=0.5)
    axis(1, at=ats, labels=xAxisLabels, cex.lab=1.5, cex.axis=1.5)
    if (jitter == TRUE){
        set.seed(420)
        stripchart(values, vertical=TRUE, add=TRUE, method="jitter", col='#404040', bg=rgb(red=0.25, green=0.25, blue=0.25, alpha=0.25), pch=21)
    }
    p <- recordPlot()
    dev.off()
    return(p)
}



get_corrMat <- function (X, y=NA, selected_col, type=c('p', 'r', 'none')) {
    # Separate X to that of class 0 and that of class 1
    X <- X[, selected_col]
    if (type=='p') {
        # By wilcox
        X0 <- X[y==0, ]
        X1 <- X[y==1, ]
        value <- sapply(seq(ncol(X0)), function(x) {wilcox.test(X0[, x], X1[, x])$p.value})
        value <- -log10(value)
    }
    if (type=='r') {
        # By spearman r
        corrY <- cor(scale(X), y, method='spearman')
        value <- corrY
        }
    if (type=='none') {
        # By mean values
        value <- apply(X, 2, mean)
        }
    # Get the stimulanting coniditions and signal repsonses responses
    stiVec <- sapply(colnames(X), function(x) {strsplit(x, '_')[[1]][1]})
    cellTypeVec <- sapply(colnames(X), function(x) {strsplit(x, '_|\\.\\.\\.')[[1]]})
    cellTypeVec <-  sapply(cellTypeVec, function(x) {x[(length(x)-1)]})

    # Binding them together
    plotDF <- data.frame(stiVec, cellTypeVec, value)

    # Create blank correlation matrix
    plotColNames <- as.character(unique(plotDF$stiVec))
    plotRowNames <- as.character(unique(plotDF$cellTypeVec))
    plotMat <- matrix(NA, length(plotRowNames), length(plotColNames))
    colnames(plotMat) <- plotColNames
    rownames(plotMat) <- plotRowNames

    # Assigning P value to the corresponding element of the matrix (by Sti and signal)
    for (i in plotRowNames) {
        for (j in plotColNames) {
        plotMat[i, j] <- plotDF[(plotDF$stiVec==j & plotDF$cellTypeVec==i), 'value']
        }
    }

    # The order that we want
    cellOrder <- c('B',
                    # NB Lymph
                    'CD4pT', 'Activated.CD4', 'Central.Memory.CD4', 'Effector.CD4', 'Effector.Memory.CD4', 'Naive.CD4', 'Tregs', #CD4
                    'CD8pT', 'Activated.CD8.cells', 'Central.Memory.CD8', 'Effector.CD8', 'Effector.Memory.CD8', 'Naive.CD8', #CD8
                    'CD4nCD8nT', 'CD4pCD8pT', #Other T Cells
                    'NKT', # others
                    'B.cells', 'IgApB.cells', 'IgDpMemory.B', 'IgDnCD27nB.cells', 'Naive.B.cells', 'Plasmablast', 'Switched.Memory.B.cells', 'Transitional.B.cells', #BCells
                    'NK.cells', 'CD16hi.NK', 'CD56.Bright.NK', 'CD56.dimnCD16.dim.NK', #NK
                    'DCs', 'mDC', 'pDC', #DC
                    # NB Mono
                    'Mono', 'CD16.hi.Mono', 'CD16.lo.Mono'
                    )
    # Rearrange the name to the desired one
    plotMat <- plotMat[match(cellOrder, rownames(plotMat)), ]
return(plotMat)
}



plot_network <- function(X, y, Xsize=NA, ysize=NA, gradient=c(NA, 'r'), reducedTo=NA, layout=NA, 
                        labels=FALSE, edgeby=c('r', 'p'), edgeCutoff=if(edgeby == 'r') {0.8} else {0.05},  
                        makeClus=FALSE, clusters=NA, cluster_palette=NA, individual_clusters=FALSE,
                        individual_clusters_output=NA, constantNodeSize=NA, nodeshape=NA) {
    # This function use, primarily, X and y arguments to contruct correlation network. Other arguments are optional
    # to tailor the graph as necessary. 
    #
    # Xsize, ysize  :Specifying the nodesize of the network (by spearman correlation). This flag is overriden if 
    #                constantNodeSize is specified. If X and y are either only have two classes or continuous, 
    #                X and y are taken as Xsize and ysize.
    # gradient      :Specifying the gradient color of the nodes using the continuous vector. The color scheme is 
    #                blue-white-red, this has to be fixed manually in the function if you wish to change the color.
    #                If the gradient is not specified, all nodes are colored black automatically.
    # reducedTo     :If specified, only the first max(abs(gradient))[1:reducedTo] will be colored, the rest is white.
    # layout        :Position of the nodes. The matrix/dataframe provided must be two columns of X and Y position.
    #                The order of the row in layout must correspond to the order of the column feature of X. If the
    #                layout is not provided, it is calculated using t-SNE on the Spearman's correlation matrix.
    # labels        :If TRUE, the feature names will be displayed on the node.
    # edgeby        :If the edge will be presented based on bonferroni correct P or Spearman r between features of X.
    # edgeCutoff    :The cutoff if the edge would present (lower than for corrected P, and higher than for r).
    # makeClus      :If specified the network will be colored by the provided list of cluster,
    #                this flag overrides the colors by gradient.
    # clusters      :A vector/list of cluster number each column of X belongs to.
    # cluster_palette: The color for the specified clusters, but be the same number as the provided unique clusters.
    #                  If not specified, it is generated automatically.
    # individual_clusters: If TRUE, images of individual cluster (equal to the number of unique clusters provided)
    #                      are generated and exported right away without being returned. If gradient and reducedTo 
    #                      is provided, the cluster color is overlayed by the reduced gradient colors.
    # individual_clusters_output: the output directory for the inidividual clusters.
    # constantNodeSize: Needs to specified if Xsize and ysize are not provided and if X and y are corresponded to 
    #                   multiple class factors. This flag set the size of all nodes, can be a vector of the same length
    #                   as the number of X columns, or a single value.
    #
    # Return an igraph graph object (with empty names of labels=FALSE)

    # Generate correlation matrix and its P values ---------------------------------------------------------------------
    corr <- cor(X, method="spearman")
    corr[is.na(corr)] <- 0
    z <- sqrt((nrow(df)-3)/1.06)*atanh(corr)
    corrP <- 2*pnorm(-abs(z))
    diag(corrP) <- 1
    testPair <- choose(ncol(corrP), 2)

    # Checking the given flags and possible conflicts ------------------------------------------------------------------
    if (any(is.na(layout))) {
        set.seed(9999)
        corr_layout <- Rtsne(corr, dims=2, perplexity=50, verbose=TRUE,
                            max_iter=1000, initial_dims=50, check_duplicates=FALSE, theta=0.1,
                            num_threads=60, partial_pca=TRUE)
        layout <- corr_layout$Y
    }
    if (any(!is.na(Xsize))) {
        if (any(colnames(X) != colnames(Xsize))) {
            stop('column names of Xsize must be the same as X.')
        }
    }
    if ((length(unique(y)) == 2)) {
        Xsize <- X
        ysize <- y
    }

    if (!is.na(reducedTo) & any(is.na(gradient))) {
        stop('A gradient must be "r" or provided (without NA) if reducedTo is specified.')
    }
    # Generating edges -------------------------------------------------------------------------------------------------
    # Edge by bonferroni P
    if (edgeby == 'p') {
        alpha <- edgeCutoff
        adjustedAlpha <- alpha/testPair
        newCorr <- corrP
        newCorr[newCorr>adjustedAlpha] <- 0
        newCorr[newCorr>0] <- 1
    }
    # Edge by r2
    if (edgeby == 'r') {
        r2cutoff <- edgeCutoff
        newCorr <- abs(corr)
        diag(newCorr) <- 0
        newCorr[newCorr<r2cutoff] <- 0
        newCorr[newCorr>r2cutoff] <- 1
    }

    # Melting the correlation networks to pairs of nodes ---------------------------------------------------------------
    meltedNewCorr <- melt(newCorr)
    edgeList <- meltedNewCorr[meltedNewCorr$value==1, c(1,2)]
    edgeList <- data.frame(edgeList)
    colnames(edgeList) <- c('from', 'to')

    # Converting cluster to corresponding palette color ----------------------------------------------------------------
    if (makeClus == TRUE) {
        if ((any(is.na(clusters)) | (length(clusters) != ncol(X)))) {
            stop('Clusters argument is not clearly defined')
        }
        if (any(is.na(cluster_palette))) {
            set.seed(7)
            qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
            col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
            set.seed(7)
            cluster_palette <- sample(col_vector, length(unique(clusters)))
        }
        if (length(cluster_palette) > length(unique(clusters))) {
            cluster_palette <- cluster_palette[1:length(unique(clusters))]
        }
    clusColor <- as.character(factor(clusters, labels=cluster_palette))
    clusColor <- cbind(names(clusters), clusColor)
    }

    # Getting node sizes -----------------------------------------------------------------------------------------------
    if (!is.na(constantNodeSize)) {
        nodeSize <- rep(constantNodeSize, ncol(X))
    } else if (!any(is.na(Xsize) | is.na(ysize))) {
        corrY <- cor(Xsize, ysize, method='spearman')
        corrY[is.na(corrY)] <- 0
        zY <- sqrt((nrow(df)-3)/1.06)*atanh(corrY)
        corrYP <- 2*pnorm(-abs(zY))
        nodeSize <- -log10(corrYP)
    } else {
        stop('Y has more than two classes or there is NA.')
    }
    
    # Combining nodes and edges ----------------------------------------------------------------------------------------
    nodes <- colnames(X)
    nodeList <- data.frame(nodes, layout)
    colnames(nodeList) <- c('nodes', 'x', 'y')

    # Assigning sizes and colors to the nodes --------------------------------------------------------------------------
    if (any(is.na(gradient))) {
        nodeList2 <- nodeList
        nodeList2[, 'gColor'] <- rep('black', ncol(X))
        nodeSize2 <- nodeSize
    } else {
        if (any(gradient == 'r')) {gradient <- c(corrY)}
        if (any(is.na(nodeshape))) {nodeshape <- rep('circle', length(corrY))}
        g <- gradient
        g[gradient > mean(gradient) + 4*sd(gradient)] <- mean(gradient) + 4*sd(gradient)
        g[gradient < mean(gradient) - 4*sd(gradient)] <- mean(gradient) - 4*sd(gradient)
        nodeList2 <- cbind(nodeList, gradient, g)
        nodeList2 <- nodeList2[order(abs(nodeList2[, 'gradient'])), ]
        nodeSize2 <- cbind(nodeSize, gradient)
        nodeSize2 <- nodeSize2[order(abs(nodeSize2[, 'gradient'])), ]
        nodeshape2 <- cbind.data.frame(nodeshape, gradient)
        nodeshape2 <- nodeshape2[order(abs(nodeshape2[, 'gradient'])), ]
        nodeList2[, 'gColor'] <- (nodeList2[, 'g'] - min(nodeList2[, 'g']))/(max(nodeList2[, 'g']) - min(nodeList2[, 'g']))
        ColorFunc <- colorRamp(c("blue", "white", "red"))
        nodeList2[, 'gColor'] <- rgb(ColorFunc(nodeList2[, 'gColor']), maxColorValue=255)
    }
    if (!is.na(reducedTo)) {
        # Reduced model only -------------------------
        upto <- reducedTo
        last <- length(nodeList2[, 'gColor'])-upto
        nodeList2[, 'gColor'][1:last] <- 'white'
        label <- rep('', 4200)
        label[nodeList2[, 'gColor'] != 'white'] <- as.character(nodeList2[, 'nodes'][nodeList2[, 'gColor'] != 'white'])
        labelNum <- rep('', 4200)
        labelNum[nodeList2[, 'gColor'] != 'white']  <- upto:1
    }

    # Generating the graph objects --------------------------------------------------------------------------------------
    if (makeClus == FALSE) {
        # Plotting just red-blue wt ---
        graph <- graph_from_data_frame(vertices=nodeList2, d=edgeList, directed=FALSE)
        E(graph)$arrow.mode <- 0
        E(graph)$color <- "light grey" #"#E8E8E8"
        E(graph)$curved <- 0
        V(graph)$size <- 3*nodeSize2[, 1] + 1
        V(graph)$color <- nodeList2[, 'gColor']
        if (any(!is.na(nodeshape))) {
            V(graph)$shape <- as.character(nodeshape2[, 1])
        }
        nodeLabel <- rep('', nrow(nodeList2))
        nodeLabel[nodeList2[, 'gColor'] != 'white'] <- as.character(nodeList2[nodeList2[, 'gColor'] != 'white', 1])
        V(graph)$label.cex <- 4
        V(graph)$label.color <- 'black'
        V(graph)$label.dist <- runif(length(corrY), -0.5, 0.5)
        if (is.na(reducedTo) == FALSE) {V(graph)$name[nodeList2[, 'gColor'] == 'white'] <- ''}
    } else if (makeClus == TRUE) {
        # Plotting clusters ---
        n_clusters <- length(unique(clusters))
        graph <- graph_from_data_frame(vertices=nodeList, d=edgeList, directed=FALSE)
        E(graph)$arrow.mode <- 0
        E(graph)$color <- "light grey" #"#E8E8E8"
        E(graph)$curved <- 0
        V(graph)$size <- nodeSize
        V(graph)$frame.color <- 'black'
        V(graph)$color <- clusColor[, 'clusColor'] 
        clusLabel <- rep('', length(clusters))
        clusLabelIndex <- sapply(seq(n_clusters), function(x) {max(which(clusters==x))})
        clusLabel[clusLabelIndex] <- 1:n_clusters
        V(graph)$name <- clusLabel
        V(graph)$label.cex <- 18
    }

    # Generating individual clusters -----------------------------------------------------------------------------------
    if (individual_clusters == TRUE) {
        if (is.na(individual_clusters_output)){
            stop('Need to specify output directory for individual cluster plots in (individual_clusters_output).')
        }
        nodeList3 <- nodeList2
        clusColor <- clusColor[match(nodeList3[, 'nodes'], clusColor[, 1]), ]
        if (!is.na(reducedTo)) {
            clusColor[, 'clusColor'][last:nrow(clusColor)] <- nodeList2[, 'gColor'][last:nrow(clusColor)]
        }
        clusters <- clusters[match(nodeList3[, 'nodes'], names(clusters))]
        nodeList3 <- cbind(nodeList2, clusColor=clusColor[, 2], nodeSize=nodeSize2[, 1], clusters)
        dir <- individual_clusters_output
        # ##---- Plotting individual clusters
        desiredClusters <- 1:n_clusters #c(16, 22, 7, 18)
        for (desiredCluster in desiredClusters) {
            desiredNodeList <- nodeList3[nodeList3['clusters']==desiredCluster, ]
            desiredEdgeList <- edgeList[(edgeList[, 'from'] %in% desiredNodeList[, 'nodes']) & (edgeList[, 'to'] %in% desiredNodeList[, 'nodes']), ]
            desiredNodeList['clusColor'] <- apply(desiredNodeList['clusColor'], 2, as.character)

            graphClus <- graph_from_data_frame(vertices=desiredNodeList, d=desiredEdgeList, directed=FALSE)
            E(graphClus)$arrow.mode <- 0
            E(graphClus)$color <- "grey"
            E(graphClus)$curved <- 0
            V(graphClus)$size <- 7*as.numeric(desiredNodeList[, 'nodeSize'])
            V(graphClus)$frame.color <- 'black'
            desiredNodeList[, 'clusColor'] <- as.character(desiredNodeList[, 'clusColor'])
            V(graphClus)$color <- desiredNodeList[, 'clusColor']

            png(filename=paste0(dir, '/cluster', desiredCluster, '.png'), width=4000, height=4000)
            plot(graphClus, vertex.label=NA)
            dev.off()
            }
        }

    if (labels == FALSE) {V(graph)$name <- ''}
    return(graph)
}

getSize <- function (pair=c(c('HC', 'AD'), c('HC', 'PD'), c('HC-II', 'DIAN'))) {
    # This function get the node size for the correlation network given a pair of interest
    # The number is based on the P value distribution in each pair. The pair has to be
    if (pair[2] == 'AD') {
      size <- c(1, 2, 3)
    }
    if (pair[2] == 'DIAN') {
      size <- c(0.5, 1.5, 2.5)
    }
    if (pair[2] == 'PD') {
      size <- c(0.5, 1.5, 2.5)
    }
    if ((pair[2] == 'HC_DIAN') | (pair[2] == 'HC-II'))  {
      size <- c(1, 2, 3)
    }
    return (size)
  }