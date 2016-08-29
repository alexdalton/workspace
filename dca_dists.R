vers = R.Version()$major
# .libPaths("/shared-mounts/sinhas/softwares/R-3.1.3/library/")
# if(vers==2){
# 	.libPaths("/shared-mounts/sinhas/lib/R.2.15/")
# }
vers = paste(sep=".","R",R.Version()$major,R.Version()$minor)
show(vers)
library(Matrix)
library(ROCR)
library(glmnet)
library(e1071)
library(proxy)

.ls.objects <- function (pos = 1, pattern, order.by, decreasing=FALSE, head=FALSE, n=5) { # order objects by size
    napply <- function(names, fn) sapply(names, function(x)
        fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.prettysize <- sapply(obj.size, function(r) prettyNum(r, big.mark = ",") )
    obj.dim <- t(napply(names, function(x)
        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size,obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    out <- out[c("Type", "PrettySize", "Rows", "Columns")]
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (head)
        out <- head(out, n)
    out
}

lsos <- function(..., n=20) { # get top objects by size
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}     

dca_dist_analysis<- function(possetfile = "C:/Users/Alex/workspace/data/msigdb/setlists/combo.txt", unifile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/hsap_universe.txt", DCANetFile = "C:/Users/Alex/workspace/data/GO_DCA_500.emb", weight = "weighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 1, nfeatures = 500, st2keep = 50, property_types = c("go_curated_evidence", "go_inferred_evidence", "pfam_domain"), writepreds = 0, outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_"){
    possetname = tail(unlist(strsplit(possetfile, "/")),1)
    possetname = gsub(".txt","",possetname)
    
    # read gene universe file
    universe = read.table(unifile)
    rownames(universe) = as.character(universe[,1])
    if(dim(universe)[2]<2){
        universe = cbind(as.character(universe[,1]), rep(1, length(universe[,1])) )
    }
    
    uniIDs = sort(as.character(unique(universe[,1])))
    dcaMatrix <- as.matrix(read.table(DCANetFile, row.names=1))
    nodeNames = rownames(dcaMatrix)
    nnodes = length(nodeNames)
    geneNodeNames = nodeNames[grep("ENSG.*", nodeNames)]
    propertyNodeNames = setdiff(nodeNames, geneNodeNames)

    # read positive set names
    possets = as.character(read.table(possetfile)[,1])

    rtable = NULL
    blankvec = structure(rep(0,length(geneNodeNames)), names=geneNodeNames)
    
    for (queryfile in possets){
        
        posset = tail(unlist(strsplit(queryfile, "/")),1)
        posset = gsub(".txt","",posset)
        
        if(file.info(queryfile)$size == 0){show(c("Empty file ", queryfile)); next}
        
        query_gs = read.table(queryfile)
        if(dim(query_gs)[2]<2){
            query_gs = cbind(as.character(query_gs[,1]), rep(1, length(query_gs[,1])) )
        }
        rownames(query_gs) = as.character(query_gs[,1])
        queryIDs = sort(intersect(geneNodeNames, as.character(unique(query_gs[,1]))))
        nquery = length(queryIDs)
        
        set.seed(041185)
        folds = sample(cut(seq(1,nquery),breaks=nfolds,labels=FALSE))
        
        for(iter in 1:nfolds){
            ## separate training and testing
            train_idxs = which(folds!=iter)
            test_idxs = which(folds==iter)
            train_nidxs = queryIDs[train_idxs]
            ntrain = length(train_nidxs)
            test_nidxs = queryIDs[test_idxs]
            ntest = length(test_nidxs)
            testuni = as.character(setdiff(geneNodeNames, train_nidxs))
            ntestuni = length(testuni)
            show(queryfile)
            
            similarityMatrix = simil(dcaMatrix, dcaMatrix[train_nidxs, ], method="cosine")
            averageSimilarities = sort(rowSums(similarityMatrix) / length(train_nidxs), decreasing=TRUE)
            selectedFeatureNames = names(averageSimilarities[1:nfeatures])
            
            # run genlasso on training set and
            y = blankvec[geneNodeNames]
            midxs = match(train_nidxs, geneNodeNames)
            y[midxs] = as.numeric(query_gs[train_nidxs,2])

            selectedSimilarityMatrix = simil(dcaMatrix[geneNodeNames, ], dcaMatrix[selectedFeatureNames, ], method="cosine")
            selectedSimilarityMatrix = selectedSimilarityMatrix[geneNodeNames, selectedFeatureNames]
            
            weights = c(nrow(selectedSimilarityMatrix) / length(train_nidxs), 1)
            names(weights) = c(1, 0)
            
            svc = svm(x = selectedSimilarityMatrix, as.factor(y), kernel='radial', gamma = 1 / length(midxs), probability = TRUE, class.weights = weights)
            pred = predict(svc,  selectedSimilarityMatrix, probability = TRUE)
            
            y = blankvec[geneNodeNames]
            midxs = match(test_nidxs, geneNodeNames)
            y[midxs] = as.numeric(query_gs[test_nidxs,2])
            model = prediction(attr(pred, "probabilities")[,2][testuni], as.factor(y)[testuni])
            auc  = performance(model, "auc")
            perf = performance(model,"tpr","fpr")
            test_aucval = round(as.numeric(slot(auc, "y.values")),3)
            show(test_aucval)
            
            y = blankvec[geneNodeNames]
            midxs = match(train_nidxs, geneNodeNames)
            y[midxs] = as.numeric(query_gs[train_nidxs,2])
            trainuni = as.character(setdiff(geneNodeNames,test_nidxs))
            model = prediction(attr(pred, "probabilities")[,2][trainuni], as.factor(y)[trainuni])
            auc  = performance(model, "auc")
            perf = performance(model,"tpr","fpr")
            train_aucval = round(as.numeric(slot(auc, "y.values")),3)
            show(train_aucval)
            
            rtable = rbind(rtable, c(queryfile, numFeatures, test_aucval, train_aucval))
        } #end iter
    } #end queryset
    colnames(rtable) = c("Query File", "Number of Features", "Test AUC", "Train AUC")
    write.table(rtable, "out2.txt", quote=F, row.names=F, sep=",")
} #end function


dca_dist_analysis(possetfile = "C:/Users/Alex/workspace/data/msigdb/setlists/combo.txt", unifile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/hsap_universe.txt", weight = "unweighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 3, st2keep = 50, nfeatures = 250, property_types = c("gene_ontology"), writepreds = 1, outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_")