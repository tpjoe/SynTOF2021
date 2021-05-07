library(sva)
library(Biobase)
# library(pvca)
library(vsn)
library(lme4)
library(ggplot2)

pvca_plt <- function(X, features, title, threshold=0.7, getValue=FALSE){
    # Plot for pvca: showing the weighted average of variance from different
    # X: a p x n matrix where n columns are number of samples
    # features: factors of covariates, for example batch, sex, etc.
    # threshold: is how many % variance PCA has to explain.
    # getValue: if you want the returns to be plot values, otherwise return actual plots.
    
    pct_threshold <- threshold
    if (is.matrix(X) == FALSE){
        X <- as.matrix(X)
    }
    
    # converting inputs to an Expression Set object
    ExpressionSet()
    expressionSet <- ExpressionSet(assayData=X)
    for (i in seq_along(features)){
        pData(expressionSet)[i] <- as.data.frame(features[i])
    }
    batch.factors <- names(pData(expressionSet))

    # run PVCA (without progress report)
    invisible(capture.output(pvcaObj <- pvcaBatchAssess(expressionSet, batch.factors, pct_threshold)))

    # plot data
    values <- as.numeric(pvcaObj[[1]])
    label <- as.data.frame(pvcaObj[[2]])
    plot_data <- cbind(values, label)
    colnames(plot_data) <- c('Value', 'Effect')
    
    # plot
    if (getValue==FALSE) {
        p <- ggplot(data=plot_data, aes(x=Effect, y=Value)) +
        ylab('Weighted average proportion variance') +
        geom_bar(stat="identity") +
        geom_text(aes(label=sprintf("%0.2f", round(values, digits=2))), vjust=-0.3, size=5) +
        ggtitle(title) + theme(plot.title=element_text(hjust=0.5)) +
        theme_bw(base_size = 25) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
        
        return(p)
    }
    # values
    if (getValue == TRUE) {
        return(plot_data)
    }
}

read_excel_allsheets <- function(filename, tibble=FALSE) {
    # This funciton import excel files that have multiple sheets into a single list of dataframes
    # filename: path to the excel file
    # I prefer straight data.frame but if you like tidyverse tibbles (the default with read_excel)
    # then just pass tibble = TRUE
    sheets <- readxl::excel_sheets(filename)
    x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet=X, na=c('n/a', '')))
    if(!tibble) x <- lapply(x, as.data.frame)
    names(x) <- sheets
    return(x)
}

combat_correction <- function(df, correct_col='Batch', start_column=7, mod=NULL) {
    
    # This funciton perform batch correction
    # df: data frame of values to be corrected with a column of covariate
    # return combat corrected dataframe

    correct_col <- as.data.frame(df[[correct_col]])
    correct_col_numeric <- as.vector(sapply(correct_col, as.numeric)) 
    X <- t(scale(as.matrix(df[ , start_column:ncol(df)]))) # scale as suggested by ComBat author

    # using ComBat
    combat <- sva::ComBat(dat=X, batch=correct_col_numeric, par.prior=TRUE, mean.only=FALSE, mod=mod)
    combat <- t(combat) # transpose it back

    return(cbind(df[, 1:(start_column - 1)], combat))
}

checkPriorName <- function(nameVec, StiName, cellName, respName) {
  # This function makes sure that the names (nameVec) extracted from data's column names 
  # belong to one of the unique names (StiName, cellName, and respName) of the prior table.
  if (((nameVec[1] %in% StiName) & (nameVec[2] %in% cellName) & (nameVec[3] %in% respName)) == FALSE) {
    return('Bad')
  }
  else {
    return('Good')
  }
}

pvcaBatchAssess <- function (abatch, batch.factors, threshold)
{
    theDataMatrix <- exprs(vsn2(abatch, verbose=FALSE, minDataPointsPerStratum=10))
    dataRowN <- nrow(theDataMatrix)
    dataColN <- ncol(theDataMatrix)

    ########## Center the data (center rows) ##########
    theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, ncol = dataColN)
    theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
    theDataMatrixCentered = t(theDataMatrixCentered_transposed)

    ########## Compute correlation matrix &  Obtain eigenvalues ##########

    theDataCor <- cor(theDataMatrixCentered)
    eigenData <- eigen(theDataCor)
    eigenValues = eigenData$values
    ev_n <- length(eigenValues)
    eigenVectorsMatrix = eigenData$vectors
    eigenValuesSum = sum(eigenValues)
    percents_PCs = eigenValues /eigenValuesSum 

    ##===========================================
    ##	Getting the experimental information
    ##===========================================
    expInfo <- pData(abatch)[,batch.factors]
    exp_design <- as.data.frame(expInfo)
    expDesignRowN <- nrow(exp_design)
    expDesignColN <- ncol(exp_design)

    ########## Merge experimental file and eigenvectors for n components ##########

    my_counter_2 = 0
    my_sum_2 = 1
    for (i in ev_n:1){
        my_sum_2  = my_sum_2 - percents_PCs[i]
        if ((my_sum_2) <= threshold ){
            my_counter_2 = my_counter_2 + 1
        }

    }

    if (my_counter_2 < 3){
        pc_n  = 3
    }else {
        pc_n = my_counter_2 
    }

    ## pc_n is the number of principal components to model

    pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN*pc_n), ncol = 1)
    mycounter = 0
    for (i in 1:pc_n){
        for (j in 1:expDesignRowN){
            mycounter <- mycounter + 1
            pc_data_matrix[mycounter,1] = eigenVectorsMatrix[j,i]
        }
    }

    AAA <- exp_design[rep(1:expDesignRowN,pc_n),]
    Data <- cbind(AAA,pc_data_matrix)


    ####### Edit these variables according to your factors #######

    variables <-c (colnames(exp_design))
    for (i in 1:length(variables))
    {
        Data$variables[i] <- as.factor(Data$variables[i] )
    }


    ########## Mixed linear model ##########
    op <- options(warn = (-1)) 
    effects_n = expDesignColN  + choose(expDesignColN, 2) + 1
    randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)

    ##============================#
    ##	Get model functions
    ##============================#
    model.func <- c()
    index <- 1

    ##	level-1
    for (i in 1:length(variables))
    {	
        mod = paste("(1|", variables[i], ")",   sep="")
        model.func[index] = mod
        index = index + 1
    }

    ##	two-way interaction
    for (i in 1:(length(variables)-1))
    {	
        for (j in (i+1):length(variables))
        {
            mod = paste("(1|", variables[i], ":", variables[j], ")",   sep="")
            model.func[index] = mod
            index = index + 1
        }
    }

    function.mods <- paste (model.func , collapse = " + ")

    ##============================#
    ##	Get random effects	#
    ##============================#

    for (i in 1:pc_n){
        y = (((i-1)*expDesignRowN)+1)
        funct <- paste ("pc_data_matrix", function.mods, sep =" ~ ")
                Rm1ML <- lmer( funct ,
                              Data[y:(((i-1)*expDesignRowN)+expDesignRowN),],
                              REML = TRUE, verbose = FALSE, na.action = na.omit)
                randomEffects <- Rm1ML
                randomEffectsMatrix[i,] <- c(unlist(VarCorr(Rm1ML)),resid=sigma(Rm1ML)^2)
    }
    effectsNames <- c(names(getME(Rm1ML,"cnms")),"resid")

    ########## Standardize Variance ##########
    randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        mySum = sum(randomEffectsMatrix[i,])
        for (j in 1:effects_n){
            randomEffectsMatrixStdze[i,j] = randomEffectsMatrix[i,j]/mySum	
        }
    }

    ########## Compute Weighted Proportions ##########

    randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
    for (i in 1:pc_n){
        weight = eigenValues[i]/eigenValuesSum
        for (j in 1:effects_n){
            randomEffectsMatrixWtProp[i,j] = randomEffectsMatrixStdze[i,j]*weight
        }
    }

    ########## Compute Weighted Ave Proportions ##########

    randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
    randomEffectsSums <-colSums(randomEffectsMatrixWtProp)
    totalSum = sum(randomEffectsSums)
    randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, ncol = effects_n)

    for (j in 1:effects_n){
        randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum 	
    }
    return(list(dat=randomEffectsMatrixWtAveProp, label=effectsNames))
}

matrixToFCS <- function (matrix) {
    if (is.null(colnames(matrix))) {
        stop("Missing colnames in matrix")
    }
    metadata <- data.frame(name=colnames(matrix), desc=colnames(matrix))
    metadata$minRange <- apply(matrix, 2, min)
    metadata$maxRange <- apply(matrix, 2, max)
    matrix <- data.matrix(matrix)

    fcs <- new("flowFrame", exprs=matrix, parameters=AnnotatedDataFrame(metadata), description=list(markers=colnames(matrix)))
}






library(gtools);

##to improve the runtime, the for loop here and in the next function should be removed (bottleneck identified during code profiling). However, it's not trivial how because of the complicated structure of the last line. Using global variables (e.g., for "a") with sapply() makes things even slower.
cmpclust <- function(classes, clusters){
	C=max(classes);
	K=max(clusters);
	n=length(classes);
	a=matrix(0,C,K);
	for (i in 1:n){
		if (classes[i]==0 && clusters[i]==0){
			next;
		}
		if (classes[i]==0){
			next;
		}
		if (clusters[i]==0){
			tr=as.integer(runif(1,1,K+1));
			a[classes[i], tr]=a[classes[i], tr]+1;
			next;
		}
		a[classes[i], clusters[i]]=a[classes[i], clusters[i]]+1;
	}
	return(a);
}


softcmpclust <- function(classes, clusters){
	C=max(classes);
	K=length(clusters[1,]);
	n=length(classes);
	a=matrix(0,C,K);
	for (i in 1:n){
		if (classes[i]==0){
			next;
		}
		if (sum(clusters[i,])<0.1){
                  tr=as.integer(runif(1,1,K+1));
                        clusters[i,]=array(0,K);
                        clusters[tr]=1;
		}
		a[classes[i], ]=a[classes[i], ]+clusters[i,];
	}
	return(a);
}


#these functions are very light since a is usually smaller than 10X10
R <- function(a, i, j){
  sum <- 0;
  for (l in 1:length(a[i,])){
    sum <- sum+a[i,l];
  }
  if (sum==0)
    return (0);
  return(a[i,j]/sum);
}
P <- function(a, i, j){
  sum <- 0;
  for (l in 1:length(a[,j])){
    sum <- sum+a[l,j];
  }
  if (sum==0)
    return (0);
  return(a[i,j]/sum);
}
F <- function(a, i, j){
  if ((R(a,i,j)+P(a,i,j))==0){
    return (0);
  }
  if ((R(a,i,j)+P(a,i,j))==0)
      return (0);
  return(2*R(a,i,j)*P(a,i,j)/(R(a,i,j)+P(a,i,j)));
}

FMeasure <- function(classes, clusters){
  if (sum(clusters)==0)
    return(0);

  a<-cmpclust(classes, clusters);
  ans <- 0;
  for (i in 1:length(a[,1])){
    temp <- vector();
    for (j in 1:length(a[1,])){
      temp[j] <- F(a,i,j);
    }
    ans <- ans + sum(a[i,])/sum(a)*max(temp);
  }
  return(ans);
}

precisionRecall <- function(classes, clusters){
  if (sum(clusters)==0)
    return(0);

  a<-cmpclust(classes, clusters);
  Recall <- 0;
  Precision <- 0;
  for (i in 1:length(a[,1])){
    temp <- vector();
    for (j in 1:length(a[1,])){
      temp[j] <- F(a,i,j);
    }
    Recall <- Recall + R(a,i,which.max(temp))*sum(a[i,])/sum(a);
    Precision <- Precision + P(a,i,which.max(temp))*sum(a[i,])/sum(a);
  }
  return(list(Precision=Precision, Recall=Recall));
}

FMeasure <- function(classes, clusters){
  if (sum(clusters)==0)
    return(0);

  a<-cmpclust(classes, clusters);
  ans <- 0;
  for (i in 1:length(a[,1])){
    temp <- vector();
    for (j in 1:length(a[1,])){
      temp[j] <- F(a,i,j);
    }
    ans <- ans + sum(a[i,])/sum(a)*max(temp);
  }
  return(ans);
}

AllFMeasures <- function(classes, clusters){
  if (sum(clusters)==0)
    return(0);

  a<-cmpclust(classes, clusters);
  ans <- 0;
  individuals <- vector()
  indexes <- vector()
  for (i in 1:length(a[,1])){
    temp <- vector();
    for (j in 1:length(a[1,])){
      temp[j] <- F(a,i,j);
    }
    ans <- ans + sum(a[i,])/sum(a)*max(temp);
    individuals[i] <- max(temp)
    indexes[i] <- which.max(temp)
  }
  return(list(individuals=individuals,indexes=indexes));
}

SoftFMeasure <- function(classes, clusters){
  a<-softcmpclust(classes, clusters);
  ans <- 0;
  for (i in 1:length(a[,1])){
    temp <- vector();
    for (j in 1:length(a[1,])){
      temp[j] <- F(a,i,j);
    }
    ans <- ans + sum(a[i,])/sum(a)*max(temp);
  }
  return(ans);
}

SoftToHard <- function(clusters){
  ans <- vector();
  for (i in 1:length(clusters[,1])){
    ##if (sum(clusters[i,])<0.01){
    ##  ans[i]=as.integer(runif(1,1,length(clusters[1,])+1));      
    ##  next;
    ##}
    ans[i] <- which.max(clusters[i,]);
  }
  ans
}

homogenieity <- function(classes, clusters, a){
	C=max(classes);
	K=max(clusters);
	#n=length(classes);
        n=sum(a);
	sum=0;
	for (k in 1:K){
		for (c in 1:C){	
			if (a[c,k]==0){
				next
			}		
			sigmaC=0;
			for (i in 1:C){
				sigmaC=sigmaC+a[i,k];
			}
			sum=sum+(a[c,k]/n)*log(a[c,k]/sigmaC,2);
		}	
	}
	Hck=sum;

	sum=0;
	for (c in 1:C){
		sigmaK=0;
		for (i in 1:K){
			sigmaK=sigmaK+a[c,i];
		}
		if (sigmaK==0){
			next
		}	
		sum=sum+(sigmaK/n)*log(sigmaK/n,2);
	}	
	Hc=sum;

	Ans=0;
	Ans$Hck=Hck;
	Ans$Hc=Hc;
	if (Hc==0){
		Ans$Homogenieity=1;
	}else{
		Ans$Homogenieity=1-Hck/Hc;
	}
	return (Ans);
}

completeness <- function(classes, clusters, a){

        C=max(classes);
	K=max(clusters);
	#n=length(classes);
        n=sum(a);
	sum=0;
	for (c in 1:C){
		for (k in 1:K){			
			if (a[c,k]==0){
				next
			}		
			sigmaK=0;
			for (i in 1:K){
				sigmaK=sigmaK+a[c,i];
			}
			sum=sum+(a[c,k]/n)*log(a[c,k]/sigmaK,2);
		}	
	}
	Hkc=sum;
      
	sum=0;
	for (k in 1:K){
		sigmaC=0;
		for (i in 1:C){
			sigmaC=sigmaC+a[i,k];
		}
		if (sigmaC==0){
			next
		}	
		sum=sum+(sigmaC/n)*log(sigmaC/n,2);
	}	
	Hk=sum;

	Ans=1;
	Ans$Hkc=Hkc;
	Ans$Hk=Hk;

	if (Hk==0){
		Ans$Completeness=1;
	}else{
		Ans$Completeness=1-Hkc/Hk;
	}
	return (Ans);
}

VMeasure <- function(classes, clusters){
	a<-cmpclust(classes, clusters)
	h <- homogenieity(classes, clusters, a);
	H <- h$Homogenieity;
	c <- completeness(classes, clusters, a);
	C <- c$Completeness;
	VM=2*H*C/(H+C);
	Ans=1;
	Ans$H=H;
	Ans$C=C;
	Ans$VM=VM;
	return (Ans);
}
