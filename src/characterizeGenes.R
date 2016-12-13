library(ROCR)


# Default input variables
# geneFeaturesFile = "C:/Users/Alex/workspace/data/test/testFeatures/testFeatures.rds"
# possetFile = "C:/Users/Alex/workspace/data/test/setlists/tests.txt"

geneFeaturesFile = "/workspace/data/smoothedNetworks/Human.GO.TM.H.Loose.Smoothed.KeptTop20.rds" # RDS file containing the input matrix, col = feature, row = sample
possetFile = "/workspace/data//setlists/msigdb_combo.txt"                                        # Text file with list of gene sets to test
resultsDir = "/workspace/results/"                                                               # Directory to put results into
n_folds = 5                                                                                      # Number of folds for k-fold cross validation
alg_params = list(algorithm = 1, eta = 0.05, gamma = 0.15, max_depth = 15,                       # Params for the alg that will affect the results. Alg ID followed by the params.
                  subsample = 0.7, objective = "binary:logistic", scale_pos_weight = TRUE, eval_metric="auc", nround = 80, base_score = 0.5)
extra_params = list(silent = 1, nthread = 8, verbose = 0)                                        # Extra params for the alg that doesn't affect results (printing, threads, etc...)


characterize <- function(featuresFile = geneFeaturesFile, possetsFile = possetFile, nfolds = n_folds, algParams = alg_params, extraParams = extra_params) {
    # Quick check if given a valid algorithm value = 1 or 2
    if (!algParams$algorithm %in% c(1, 2)) {
        show("Invalid Algorithm")
        return(1) 
    }

    geneFeatures = readRDS(featuresFile) # Read each gene's features from the input file
    geneNames = rownames(geneFeatures) # Get names of each gene from row names
    featureNames = colnames(geneFeatures) # Get names of each feature type from col names
    numGenes = length(geneNames)
    numFeatures = length(featureNames)
    
    featuresFileName = gsub(".rds", "", tail(unlist(strsplit(featuresFile, "/")), 1))
    possetFileName = gsub(".txt", "", tail(unlist(strsplit(possetsFile, "/")), 1))
    
    resultsTable = NULL
    resultsFile = paste(sep = "", resultsDir, featuresFileName, "_", possetFileName, "_", gsub(":", "", paste(algParams, collapse = "_")), "_.csv")

    if (file.exists(resultsFile)) {
        show(c("Results File Already Exists", resultsFile))
        return(1)
    } else { show(resultsFile) }
    
    possets = as.character(read.table(possetsFile)[,1]) # Read set lists from posset file
    blankvec = structure(rep(0, length(geneNames)), names = geneNames) # General use vector of all genes with all values set to 0
    
    # Perform preliminary algorithm specific tasks
    if (algParams$algorithm == 1) {
        library(xgboost) # Load the XGBoost library
        dFullTrain = xgb.DMatrix(geneFeatures) # Convert input matrix to DMatrix for XGBoost
    } else if (algParams$algorithm == 2) {
        library(e1071) # Load the library for SVM
    }
    
    # Perform characterization for each queryFile
    for (queryFile in possets) {
        # Get the name of the query
        posset = tail(unlist(strsplit(queryFile, "/")),1)
        posset = gsub(".txt","",posset)
        if(file.info(queryFile)$size == 0){show(c("Empty file ", queryFile)); next}
        show(posset)
        
        # Read query file for positive labeled genes
        query_gs = read.table(queryFile)
        if(dim(query_gs)[2]<2){
            query_gs = cbind(as.character(query_gs[,1]), rep(1, length(query_gs[,1])) )
        }
        rownames(query_gs) = as.character(query_gs[,1])
        queryIDs = sort(intersect(geneNames, as.character(unique(query_gs[,1])))) # Row index of each query gene in geneFeatures
        nquery = length(queryIDs)
        
        set.seed(041185) # Default seed for repeatability
        
        folds = sample(cut(seq(1, nquery), breaks = nfolds, labels = FALSE)) # Split the query into n folds
        test_scores = NULL
        train_scores = NULL
        
        # Perform characterization on each fold of the query
        for(iter in 1:nfolds) {
            train_idxs = which(folds != iter)  # Indices into queryIDs for the training genes
            test_idxs = which(folds == iter)   # Indices into queryIDs for the test genes
            train_nidxs = queryIDs[train_idxs] # Names of genes in train set
            test_nidxs = queryIDs[test_idxs]   # Names of genes in test set
            ntrain = length(train_nidxs)
            ntest = length(test_nidxs)
            
            testuni = as.character(setdiff(geneNames, train_nidxs)) # Names of all genes within the test universe (all genes minus the training genes from the query)
            trainuni = as.character(setdiff(geneNames, test_nidxs)) # Names of all genes within the train universe (all genes minus the testing genes from the query)
            
            # Create a gene vector with values set as the training labels
            train_labels = blankvec
            train_labels[match(train_nidxs, geneNames)] = as.numeric(query_gs[train_nidxs, 2])
            
            # Create a gene vector with values set as the testing labels
            test_labels = blankvec
            test_labels[match(test_nidxs, geneNames)] = as.numeric(query_gs[test_nidxs,2])
            
            # Perform XGBoost Characterization
            if (algParams$algorithm == 1) {
                setinfo(dFullTrain, 'label', train_labels) # Set the labels for training
                
                # Set XGBoost specific param list
                param <- list(eta = algParams$eta, gamma = algParams$gamma, max_depth = algParams$max_depth, subsample = algParams$subsample, 
                              objective = algParams$objective, eval_metric=algParams$eval_metric, silent = extraParams$silent)
                if (algParams$scale_pos_weight) {
                    param$scale_pos_weight = (numGenes - ntrain) / ntrain 
                }
                
                # Train XGBoost classifier
                bst = xgb.train(param, dFullTrain, nthread = extraParams$nthread, nround = algParams$nround,
                                verbose = extraParams$verbose,  base_score = algParams$base_score)
                
                # Use classifier to return probability of each gene being a part of the query
                pred <- predict(bst, geneFeatures)
            }
            
            # Calculate AUC score for the test universe
            test_uni_idxs = match(testuni, geneNames)
            model = prediction(pred[test_uni_idxs], as.factor(test_labels)[test_uni_idxs])
            auc  = performance(model, "auc")
            test_aucval = round(as.numeric(slot(auc, "y.values")), 5)
            test_scores = c(test_scores, test_aucval)
            
            pdf(paste(sep = "", resultsDir, "plots/", featuresFileName, "_", possetFileName, "_", gsub(":", "", paste(algParams, collapse = "_")),
                      "_", posset, "_", iter, ".pdf"))
            plot(performance(model, measure = "tpr", x.measure = "fpr"))
            dev.off()
            
            # Calculate AUC score for the train universe
            train_uni_idxs = match(trainuni, geneNames)
            model = prediction(pred[train_uni_idxs], as.factor(train_labels)[train_uni_idxs])
            auc  = performance(model, "auc")
            train_aucval = round(as.numeric(slot(auc, "y.values")), 5)
            train_scores = c(train_scores, train_aucval)
            
            show(c(test_aucval, train_aucval))
        } # end loop over folds
        row = c(posset, test_scores, mean(test_scores), train_scores, mean(train_scores))
        resultsTable = rbind(resultsTable, row)
    } # end loop over query files
    
    # TODO HARDCODED FOR 5 FOLDS!!
    colnames(resultsTable) = c("Query Name", "F1 Test", "F2 Test", "F3 Test", "F4 Test", "F5 Test", "Mean Test",
                               "F1 Train", "F2 Train", "F3 Train", "F4 Train", "F5 Train", "Mean Train")
    write.table(resultsTable, resultsFile, quote = FALSE, sep = ",", row.names = FALSE)
} # end function characterize

characterize()