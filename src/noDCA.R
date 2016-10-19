library(xgboost)
library(ROCR)
library(e1071)


# Input variables
smoothedNetworkFile = "/workspace/R-code/Human.GO.TM.H.Loose.Smoothed.Booleanized.rds"
possetFile = "/workspace/R-code/data/msigdb/setlists/combo.txt"
nfolds = 5
nXGBthreads = 8

smoothedNetwork = readRDS(smoothedNetworkFile)
geneNames = rownames(smoothedNetwork)
featureNames = colnames(smoothedNetwork)
numGenes = length(geneNames)
numFeatures = length(featureNames)

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
        param <- list(eta=.15, gamma=.15, max_depth=numFeatures, subsample=1, objective="binary:logistic", silent=1, eval_metric="auc")
        nround = 50
        base_score = .5
        negWeight = .02
        weights = blankvec
        weights[testuni] = negWeight
        weights[train_nidxs] = 1
        setinfo(dFullTrain, "weight", weights)
        #watchlist <- list(eval = dtest, train = dFullTrain)

        bst = xgb.train(param, dFullTrain, nthread=nXGBthreads, nround=nround, verbose=0,  base_score=base_score)

        pred <- predict(bst, smoothedNetwork)

        test_uni_idxs = match(testuni, geneNames)
        model = prediction(pred[test_uni_idxs], as.factor(test_labels)[test_uni_idxs])
        auc  = performance(model, "auc")
        test_aucval = round(as.numeric(slot(auc, "y.values")), 5)
        show("XGBoost Test AUC")
        show(test_aucval)

        train_uni_idxs = match(trainuni, geneNames)
        model = prediction(pred[train_uni_idxs], as.factor(train_labels)[train_uni_idxs])
        auc  = performance(model, "auc")
        train_aucval = round(as.numeric(slot(auc, "y.values")), 5)
        show("XGBoost Train AUC")
        show(train_aucval)

        # featureImportance = xgb.importance(feature_names=featureNames, model=bst)
        # topFeatures = featureImportance$Feature
        # featIdxs = match(topFeatures, featureNames)

        }
    }
}
