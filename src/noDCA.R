library(xgboost)
library(ROCR)
#library(e1071)


# Input variables
smoothedNetworkFile = "/workspace/data/msigdb/smoothedNetworks/Human.GO.TM.H.Loose.Smoothed.KeptTop20.rds"
possetFile = "/workspace/data/msigdb/setlists/combo.txt"
nfolds = 5
nXGBthreads = 8

smoothedNetwork = readRDS(smoothedNetworkFile)
geneNames = rownames(smoothedNetwork)
featureNames = colnames(smoothedNetwork)
numGenes = length(geneNames)

dFullTrain = xgb.DMatrix(smoothedNetwork)
possets = as.character(read.table(possetFile)[,1])
blankvec = structure(rep(0,length(geneNames)), names = geneNames)

for (queryFile in possets){
    posset = tail(unlist(strsplit(queryFile, "/")),1)
    posset = gsub(".txt","",posset)

    if(file.info(queryFile)$size == 0){show(c("Empty file ", queryFile)); next}

    query_gs = read.table(queryFile)
    if(dim(query_gs)[2]<2){
        query_gs = cbind(as.character(query_gs[,1]), rep(1, length(query_gs[,1])) )
    }

    rownames(query_gs) = as.character(query_gs[,1])
    queryIDs = sort(intersect(geneNames, as.character(unique(query_gs[,1]))))
    nquery = length(queryIDs)

    set.seed(041185)
    folds = sample(cut(seq(1,nquery),breaks=nfolds,labels=FALSE))

    for(iter in 1:nfolds) {
        numFeatures = length(featureNames)
        train_idxs = which(folds!=iter)
        test_idxs = which(folds==iter)
        train_nidxs = queryIDs[train_idxs]
        ntrain = length(train_nidxs)
        test_nidxs = queryIDs[test_idxs]
        ntest = length(test_nidxs)
        testuni = as.character(setdiff(geneNames, train_nidxs))
        trainuni = as.character(setdiff(geneNames, test_nidxs))
        show(posset)

        train_labels = blankvec
        midxs = match(train_nidxs, geneNames)
        train_labels[midxs] = as.numeric(query_gs[train_nidxs, 2])

        test_labels = blankvec
        midxs = match(test_nidxs, geneNames)
        test_labels[midxs] = as.numeric(query_gs[test_nidxs,2])

        setinfo(dFullTrain, 'label', train_labels)

        #dtest <- xgb.DMatrix(smoothedNetwork[testuni, ], label=test_labels[testuni])
        param <- list(eta=.05, gamma=.15, max_depth=15, subsample=.7, objective="binary:logistic", silent=1, scale_pos_weight=(numGenes - ntrain) / ntrain , eval_metric="auc")
        nround = 80
        base_score = .5
        # negWeight = .02
        # weights = blankvec
        # weights[testuni] = negWeight
        # weights[train_nidxs] = 1
        # setinfo(dFullTrain, "weight", weights)
        #watchlist <- list(eval = dtest, train = dFullTrain)

        #show("Getting top 100 features with XGB")
        ptm <- proc.time()[3]
        bst = xgb.train(param, dFullTrain, nthread=nXGBthreads, nround=nround, verbose=0,  base_score=base_score)

        pred <- predict(bst, smoothedNetwork)

        test_uni_idxs = match(testuni, geneNames)
        model = prediction(pred[test_uni_idxs], as.factor(test_labels)[test_uni_idxs])
        auc  = performance(model, "auc")
        test_aucval = round(as.numeric(slot(auc, "y.values")), 5)

        train_uni_idxs = match(trainuni, geneNames)
        model = prediction(pred[train_uni_idxs], as.factor(train_labels)[train_uni_idxs])
        auc  = performance(model, "auc")
        train_aucval = round(as.numeric(slot(auc, "y.values")), 5)
        show("XGBoost Test AUC, Train AUC")
        show(c(test_aucval, train_aucval))

        # featureImportance = xgb.importance(feature_names=featureNames, model=bst)
        # rm(bst)

        # topFeatures = featureImportance$Feature[1:min(100, length(featureImportance$Feature))]
        # show(proc.time()[3] - ptm)

        # show("Generating New Features")
        # ptm <- proc.time()[3]
        # featureCombs = combn(topFeatures, 2)
        # newFeatures = matrix(0, numGenes, dim(featureCombs)[2] * 4)

        # comb1Idxs = match(featureCombs[1,], featureNames)
        # comb2Idxs = match(featureCombs[2,], featureNames)

        # rm(featureCombs, featureImportance)

        # for(i in 1:numGenes) {
        #      vector1 = smoothedNetwork[i, comb1Idxs]
        #      vector2 = smoothedNetwork[i, comb2Idxs]
        #      newFeatures[i, ] = as.numeric(c(xor(vector1, vector2), vector1 & vector2, !vector1 & vector2, vector1 & !vector2))
        # }

        # rownames(newFeatures) = geneNames
        # comb1Names = featureNames[comb1Idxs]
        # comb2Names = featureNames[comb2Idxs]
        # newFeatureNames = c(mapply(paste, comb1Names, comb2Names, sep = " xor ", USE.NAMES = FALSE),
        #                     mapply(paste, comb1Names, comb2Names, sep = " and ", USE.NAMES = FALSE),
        #                     mapply(paste, comb1Names, comb2Names, sep = " not-and ", USE.NAMES = FALSE),
        #                     mapply(paste, comb1Names, comb2Names, sep = " and-not ", USE.NAMES = FALSE))

        # colnames(newFeatures) = newFeatureNames
        # rm(comb1Names, comb2Names, newFeatureNames, comb1Idxs, comb2Idxs)
        # show(proc.time()[3] - ptm)

        # show("Combining old and new features and creating DMatrix")
        # ptm <- proc.time()[3]
        # newSmoothedNetwork = cbind(smoothedNetwork, newFeatures)                # lots of memory
        # rm(newFeatures)
        # dCombtrain = xgb.DMatrix(newSmoothedNetwork, label=train_labels)        # lots of memory
        # numFeatures = ncol(newSmoothedNetwork)
        # show(proc.time()[3] - ptm)


        # param$max_depth = numFeatures
        # setinfo(dCombtrain, "weight", weights)
        # show("Training with Engineered Features")
        # ptm <- proc.time()[3]
        # bst = xgb.train(param, dCombtrain, nthread=nXGBthreads, nround=nround, verbose=0,  base_score=base_score)
        # show(proc.time()[3] - ptm)
        # rm(dCombtrain)

        # show("Running Prediction on Full Engineered Matrix")
        # ptm <- proc.time()[3]
        # pred <- predict(bst, newSmoothedNetwork)
        # show(proc.time()[3] - ptm)

        # rm(bst, newSmoothedNetwork)

        # model = prediction(pred[test_uni_idxs], as.factor(test_labels)[test_uni_idxs])
        # auc  = performance(model, "auc")
        # test_aucval = round(as.numeric(slot(auc, "y.values")), 5)
        # show("XGBoost Engineered Test AUC")
        # show(test_aucval)

        # model = prediction(pred[train_uni_idxs], as.factor(train_labels)[train_uni_idxs])
        # auc  = performance(model, "auc")
        # train_aucval = round(as.numeric(slot(auc, "y.values")), 5)
        # show("XGBoost Engineered Train AUC")
        # show(train_aucval)
        }
}
