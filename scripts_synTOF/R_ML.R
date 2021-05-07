#### libraries -----------------------------------------------------------------------
library(randomForest)
library(glmnet)
library(MLmetrics)
library(reticulate)
library(Rdimtools)


#### Importing data ------------------------------------------------------------------
regions <- c('BA9', 'DLCau', 'Hipp')
for (region in regions) {
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
}


pair <- c('LowNo', 'PHAD')
X <- df[df$group %in% pair, !colnames(df) %in% c('group', 'sample')]
X <- X[, !(seq(ncol(X)) %in% which(colSums(X==0) == nrow(X)))]
y <- df[df$group %in% pair, 'group']
y <- factor(sapply(y, function(x) if (x==pair[1]) {0} else {1}))

if (pair[2] == 'LBD'){
    for(i in 1:ncol(X)){
        X[is.na(X[, i]), i] <- mean(X[,i], na.rm = TRUE)
    }
    X <- X[, which(!(apply(X, 2, sd)==0))]
} else {
    X <- X[, !((grepl('Hipp', colnames(X))) & (grepl('p\\.Tau', colnames(X))))]
}

# X <- X[, grepl('Hipp', colnames(X))]


# leave one out -----------------------------------------------------------------------------------------
library(doParallel)
library(zetadiv)
cl <- makeCluster(20)
registerDoParallel(cl)

X <- scale(X)

y_pred <- foreach (i=seq(1:nrow(X)), .combine=rbind, .packages=c('glmnet', 'zetadiv', 'reticulate')) %dopar% {
# for (i in seq(nrow(X))) {
    sklearn <- import('sklearn')
    test_index <- i
    X_train <- as.matrix(X[-test_index, ])
    y_train <- y[-test_index]
    X_test <- t(as.matrix(X[test_index, ]))
    y_test <- y[test_index]
    # scaler <- sklearn$preprocessing$StandardScaler(with_mean=FALSE)
    # X_train <- scaler$fit_transform(X_train)
    # X_test <- scaler$transform(t(X_test))
    

    # lsvc <- sklearn$linear_model$LogisticRegression()$fit(X_train, y_train)
    # model <- sklearn$feature_selection$SelectFromModel(lsvc, prefit=TRUE)
    # X_train <- model$transform(X_train)
    # X_test <- model$transform(X_test)

    #feature selection ---
    # corrY <- cor(X_train, as.numeric(y_train), method='spearman')
    # corrY[is.na(corrY)] <- 0
    # zY <- sqrt((nrow(df)-3)/1.06)*atanh(corrY)
    # corrYP <- 2*pnorm(-abs(zY))
    # nodeSize <- -log10(corrYP)
    # X_train <- X_train[, nodeSize>1.30103]
    # X_test <- as.matrix(X_test[, nodeSize>1.30103])

    # clf <- cv.glmnet(X_train, y_train, family='binomial', nfolds=3, type.measure='deviance', standardize=FALSE, intercept=FALSE)
    # y_pred <- predict(clf, (X_test), type="response", s='lambda.min')
    # y_pred
    EN <- sklearn$linear_model$LogisticRegressionCV(fit_intercept=FALSE, solver='saga', max_iter=10000, penalty='elasticnet', l1_ratios=c(0, 0.25, 0.5, 0.75, 1))
    EN$fit(X_train, y_train)
    y_pred_ <- EN$predict_proba((X_test))[, 2]
    # clf <- randomForest(X_train, y_train, ntree=100)
    # y_pred <- c(y_pred, predict(clf, t(X_test), type="prob")[, 2])
}
AUC(y_pred, y)
wilcox.test(y_pred[y==0], y_pred[y==1])$p.value
# importance(clf)


colnames(X)[order(abs(EN$coef_), decreasing=FALSE)]
sort(abs(EN$coef_))




# leave one out with stack gen --------------------------------------------------------------------------
library(doParallel)
library(zetadiv)
cl <- makeCluster(20)
registerDoParallel(cl)


y_pred <- foreach(test_index=seq(1:nrow(X)), .combine=rbind, .packages=c('glmnet', 'zetadiv', 'reticulate')) %dopar% {
# for (test_index in c(10)) {
    sklearn <- import('sklearn')
    X_train_val <- as.matrix(X[-test_index, ])
    y_train_val <- y[-test_index]
    X_test <- as.matrix(X[test_index, ])
    y_test <- y[test_index]

    y_pred_val <- matrix(NA, nrow(X_train_val), 3)
    colnames(y_pred_val) <- c('BA9', 'DLCau', 'Hipp')
    coefs <- list()
    intercepts <- list()
    for (region in c('BA9', 'DLCau', 'Hipp')) {
        X_region <- X_train_val[, grepl(region, colnames(X_train_val))]
        coefs[[region]] <- matrix(NA, ncol(X_region), nrow(X_region))
        intercepts[[region]] <- matrix(NA, 1, nrow(X_region))
        for (val_index in seq(nrow(X_region))) {
            X_train <- as.matrix(X_region[-val_index, ])
            y_train <- y_train_val[-val_index]
            X_val <- t(as.matrix(X_region[val_index, ]))
            y_val <- y_train_val[val_index]
            # clf <- cv.glmnet(X_train, y_train, family='binomial', nfolds=nrow(X_train), type.measure='deviance', standardize=FALSE)
            # y_pred_val[val_index, region] <- predict(clf, (X_val), type="response", s='lambda.min')
            clf <- sklearn$linear_model$LogisticRegressionCV(fit_intercept=FALSE, solver='saga', max_iter=10000, penalty='elasticnet', l1_ratios=c(0, 0.25, 0.5, 0.75, 1))
            clf$fit(X_train, y_train)
            y_pred_val[val_index, region] <- clf$predict_proba(X_val)[, 2]
            coefs[[region]][, val_index] <- t(clf$coef_) #coef(clf)[-1]
            intercepts[[region]][, val_index] <- 0 #coef(clf)[1]
        }
        coefs[[region]] <- apply(coefs[[region]], 1, mean)
        intercepts[[region]] <- apply(intercepts[[region]], 1, mean)
    }
    # clf_meta <- cv.glmnet(y_pred_val, y_train_val, intercept=FALSE, lower.limits=0, family='binomial', 
    #                       nfolds=3, type.measure='deviance', alpha=0)
    clf_meta <- glm(y_train_val ~ y_pred_val + 0, family='binomial')
    # clf_meta <- glm.cons(y_train_val ~ y_pred_val + 0, family='binomial', cons=1)
    predict_prob <- c()
    for (region in c('BA9', 'DLCau', 'Hipp')) {
        X_test_region <- X_test[grepl(region, rownames(X_test)), ]
        predict_prob <- c(predict_prob, plogis(X_test_region %*% coefs[[region]] + intercepts[[region]]))
    }
    # test <- 
    plogis(t(predict_prob) %*% (coef(clf_meta)))
    # predict(clf_meta, t(predict_prob), type="response", s='lambda.min')
    # coef(clf_meta)
}

AUC(EN$predict_proba(y_pred_val)[, 1], y_train_val)


aa = predict(clf_meta, as.data.frame(y_pred_val), type="response")#, s='lambda.min')
AUC(aa, y_train_val)
AUC(y_pred_val[, 3], y_train_val)

AUC(y_pred, y)



    # clf <- cv.glmnet(X_train, y_train, family='binomial', nfolds=18, type.measure='deviance')
    # y_pred <- c(y_pred, predict(clf, (X_test), type="response", s='lambda.min'))
    # EN <- sklearn$linear_model$LogisticRegressionCV( fit_intercept=FALSE, solver='saga', max_iter=10000)
    # EN$fit(X_train, y_train)
    # y_pred <- c(y_pred, EN$predict_proba(t(X_test))[, 1])
    # set.seed(420)
    # clf <- randomForest(X_train, y_train, ntree=500)
    # y_pred <- c(y_pred, predict(clf, t(X_test), type="prob")[, 2])