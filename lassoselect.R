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

threeCol2MaxMat<- function(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4), sym = 0){
    if(length(a)==length(b) & length(a)==length(v)){
        a = as.character(a)
        b = as.character(b)
        v = as.numeric(v)
        avals = sort(unique(a))
        bvals = sort(unique(b))
        if(sym){
            avals = sort(unique(c(avals,bvals)))
            bvals = avals
        }
        nrows = length(avals)
        ncols = length(bvals)

        aidxs = match(a, avals)
        bidxs = match(b, bvals)

        vidxs = (bidxs-1)*nrows+aidxs

        retMat = Matrix(0, nrows, ncols, dimnames = list(avals,bvals), sparse=T)
        retMat[vidxs] = v
        if(sym){
            vidxs = (aidxs-1)*ncols+bidxs
            #	write over duplicates
            #			betteridxs = which(retMat[vidxs] < v)
            #			retMat[vidxs[betteridxs]] = v[betteridxs]
            retMat[vidxs] = v
        }
        return(retMat)
    }
    return(-1)
}

threeCol2listMat<- function(a = c("a","b","c","c"), b = c("a","b","b","b"), v = c(1,2,3,4), sym = 0){
    if(length(a)==length(b) & length(a)==length(v)){
        a1 = as.character(a)
        b1 = as.character(b)
        #		a1 = as.integer(a)
        #		b1 = as.integer(b)
        v1 = as.numeric(v)

        if(sym){
            show("appending...")
            a1 = c(b1,a1)
            b1 = c(as.character(a),b1)
            #			b1 = c(as.integer(a), b1)
            v1 = c(v1,v1)
        }

        show("sorting...")
        ss = sort(a1, index.return=T)
        a2 = a1[ss$ix]
        b2 = b1[ss$ix]
        v2 = v1[ss$ix]

        avals = unique(a2)
        bvals = sort(unique(b2))
        #		nvals = 10000  # for testing limit
        nvals = length(a2)
        nrows = length(avals)
        ncols = length(bvals)

        listm = list()
        colsums = 1:length(bvals)*0
        show("creating list...")
        curval = a2[1]
        idxstart = 1
        for (i in 2:nvals){
            if(i%%100000==1){show(i)}

            if(a2[i] != curval){
                # process values for row
                cols = b2[idxstart:(i-1)]
                vals = v2[idxstart:(i-1)]
                # arbitrarily select one value
                listidxs = match(unique(cols), cols)
                colidxs = match(cols[listidxs],bvals)
                listm[[as.character(curval)]] = list("colidxs" = colidxs,"vals" = vals[listidxs])
                colsums[colidxs] = colsums[colidxs] + vals[listidxs]
                # update pointers
                curval = a2[i]
                idxstart = i
            }
        }
        # process values for end
        if(a2[i]==a2[i-1]){
            cols = b2[idxstart:i]
            vals = v2[idxstart:i]
            # arbitrarily select one value
            listidxs = match(unique(cols), cols)
            colidxs = match(cols[listidxs],bvals)
            listm[[as.character(curval)]] = list("colidxs" = colidxs,"vals" = vals[listidxs])
            colsums[colidxs] = colsums[colidxs] + vals[listidxs]
        }
        return(list("listm" = listm, "avals" = avals, "colsums"=colsums))
    }
    return(-1)
}

SMOTE <- function(sampleSet, N, k) {
    T = nrow(sampleSet)
    N = N / 100

    syntheticSamples = NULL

    for (i in 1:T) {
        currentSample = sampleSet[i, ]
        distances = vector(mode = "numeric", length = T)
        for (j in 1:T) {
            otherSample = sampleSet[j, ]
            distances[j] = sqrt(sum((currentSample - otherSample)^2))
            if (j == i) {
                distances[j] = 10000000
            }
        }
        kNearest = order(distances)[1:k]
        for (randomNearest in sample(kNearest, N)) {
            neighborSample = sampleSet[randomNearest, ]
            syntheticSample = vector(mode = "numeric", length=length(neighborSample))
            for (featureIndex in 1:length(neighborSample)) {
                if (runif(1, 0, 1) > .50) {
                    syntheticSample[featureIndex] = neighborSample[featureIndex]
                }
                else {
                    syntheticSample[featureIndex] = currentSample[featureIndex]
                }
            }
            syntheticSamples = rbind(syntheticSamples, syntheticSample)
        }
    }
    rownames(syntheticSamples) = paste("syntheticSample", as.character(1:nrow(syntheticSamples)), sep="")
    return(syntheticSamples)
}

RWR<- function(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh){
    damping = 1-restart
    query = query / sum(query)
    vec = startvec
    mycolsum = NULL
    mylistm = NULL
    mismatch = 0

    if(!boolSparceMat) {
        mycolsums = transmat$colsums
        mylistm = transmat$listm
        mismatch = sum(names(query)!=names(mylistm))
    } else {
        mismatch = sum(names(query)!=colnames(transmat))
    }

    if(mismatch>0){
        show("data names do not match")
        return -1
    }

    for (i in 1:maxiters){
        newvec = NULL
        if(boolSparceMat) {
            newvec = damping*transmat%*%vec+restart*query
        } else {
            rr = sapply(mylistm,function(x){vec[x$colidxs]%*%(x$vals/mycolsums[x$colidxs])})
            newvec = damping*rr+restart*query
        }
        diff = sum(abs(newvec-vec))
        vec = newvec
        #		show(c(i,signif(diff,3)))
        if(diff < thresh) { break }
    }
    vec = as.vector(vec)
    names(vec) = names(query)
    return(list("iter" = i, "diff" = diff, "vec"=vec))
}

# show(lsos())

possetfile = "/workspace/R-code/data/msigdb/setlists/combo.txt"
unifile = "/workspace/R-code/data/msigdb/gene_sets/hsap_universe.txt"
networkfile = "/workspace/R-code/data/msigdb/networks/1sp_many/Human.GO.TM.Loose.txt"
outdir = "/workspace/R-code/data/msigdb/results/1st_"
weight = "weighted"
normalize = "type"
restarts = .7
maxiters = 50
thresh = 0.001

st2keep = 50
writepreds = 1

# queryfile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/msigdb_combo/ZHANG_BREAST_CANCER_PROGENITORS.CB.txt"
# restart = .7
# iter = 1
forcelarge = 0

property_types = c("go_curated_evidence", "go_inferred_evidence", "pfam_domain")

rwr_2stage<- function(possetfile = "/workspace/R-code/data/msigdb/setlists/combo.txt", unifile = "/workspace/R-code/data/msigdb/gene_sets/hsap_universe.txt", networkfile = "C:/Users/Alex/workspace/data/msigdb/networks/1sp_many/Human.GO.TM.Loose.txt", weight = "weighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 1, nfeatures = 500, st2keep = 50, property_types = c("go_curated_evidence", "go_inferred_evidence", "pfam_domain"), writepreds = 0, outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_"){

    uni = tail(unlist(strsplit(unifile, "/")),1)
    uni = gsub(".txt","",uni)
    network = tail(unlist(strsplit(networkfile, "/")),1)
    network = gsub(".edge","",network)
    network = gsub(".txt","",network)
    possetname = tail(unlist(strsplit(possetfile, "/")),1)
    possetname = gsub(".txt","",possetname)

    restable = NULL
    resfile = paste(sep="", outdir, uni, ".", network, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", st2keep, ".", paste(collapse="_", restarts), ".", possetname, ".stats")

    #if(file.exists(resfile)){show(c("Already exists ", resfile)); return(1)}
    write.table(restable, resfile, quote=F, sep="\t", row.names=F)
    show(resfile)

    # read in edge file
    edges = read.table(networkfile)
    colnames(edges) = c("source", "target", "weight","type")
    edges$weight <- as.numeric( as.character( edges$weight ) )

    # check for strange values, 0, negs, max
    wmin = min(edges$weight)
    wmins = edges[which(edges$weight==wmin),]
    wmins[1:min(5,length(wmins[,1])),]

    wmax = max(edges$weight)
    wmaxs = edges[which(edges$weight==wmax),]
    wmaxs[1:min(5,length(wmaxs[,1])),]
    rm(wmaxs)
    rm(wmins)

    #Stage 0

    # convert to unweighted
    if(weight == "unweighted"){
        show("Removing Edge Weights")
        edges[,"weight"] = 1
    }

    # normalize among types
    typetable = NULL
    if(normalize == "type"){
        show("Original Edge Values")
        show(edges[1:5,])
        show("Edge Type Summary")
        typetable = aggregate(weight~type, data = edges, FUN=sum)
        show(typetable)
        rownames(typetable) = as.character(typetable[,1])
        edges[,"weight"] = as.numeric(edges[,"weight"]) / as.numeric(typetable[as.character(edges[,"type"]),"weight"])
        show("Normalized Edge Values")
        show(edges[1:5,])
    }

    #output normalized matrix and stop
    #write.table(edges, paste(edgefile, ".norm.txt", sep=""), 	quote=F, sep="\t", col.names=F, row.names=F)
    #return(0)

    # collect info about feature nodes
    all_etypes = as.character(unique(edges$type))
    prop_etypes = intersect(all_etypes, property_types)
    ntypes = length(prop_etypes)
    features = unique(edges[which(edges$type %in% prop_etypes),c("source","type")])
    featnodes = as.character(features[,"source"])
    nfeats = length(featnodes)
    nkeep = min(st2keep*ntypes,nfeats)

    transmat = NULL
    nodenames = NULL
    colsum = NULL

    # covert edgelist to sparce matrix or personal list format
    # check is space matrix will work
    node_estimate = length(unique(c(as.character(edges[,1]),as.character(edges[,2]))))
    boolSparceMat = (node_estimate^2 < .Machine$integer.max) * (1-forcelarge)   # will fit in sparce mat

    if(boolSparceMat) {
        # convert edgelist to matrix
        n1Matrix = threeCol2MaxMat( as.character(edges[,"source"]),  as.character(edges[,"target"]), as.numeric(edges[,"weight"]), 1)
        nodenames = rownames(n1Matrix)

        # column normalize
        colsum = colSums(n1Matrix)
        transmat = t(t(n1Matrix)/colsum)
        rm(n1Matrix)
    } else {  # must use personal list format
        ll = threeCol2listMat( as.character(edges[,"source"]),  as.character(edges[,"target"]), as.numeric(edges[,"weight"]), 1)
        transmat = ll
        nodenames = ll$avals
        colsum = ll$colsums
        rm(ll)
    }
    nnodes = length(nodenames)

    # read gene universe file
    universe = read.table(unifile)
    rownames(universe) = as.character(universe[,1])
    if(dim(universe)[2]<2){
        universe = cbind(as.character(universe[,1]), rep(1, length(universe[,1])) )
    }
    uniIDs = sort(intersect(nodenames, as.character(unique(universe[,1]))))
    propIDs = setdiff(colnames(transmat), uniIDs)

    if(length(uniIDs)<1){ return(-1)}

    for (restart in restarts){
        blankvec = structure(rep(0,nnodes), names = nodenames)
        startvec = blankvec + 1 / nnodes

        smoothedNetFile = paste(sep="", outdir, uni, ".", network, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", st2keep, ".", paste(collapse="_", restarts), ".smoothed")
        if(!file.exists(smoothedNetFile)) {
            smoothedNetwork = matrix(0, nrow=length(uniIDs), ncol=length(propIDs))
            colnames(smoothedNetwork) = propIDs
            rownames(smoothedNetwork) = uniIDs
            i = 0
            for(uniID in uniIDs){
                query = blankvec
                query[uniID] = 1
                rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)
                smoothedNetwork[uniID, ] = rwr_res$vec[propIDs]
                i = i + 1
                show(i / length(uniIDs))
            }
            #write.table(smoothedNetwork, file=smoothedNetFile, row.names=TRUE, col.names=TRUE)
        } else {
            smoothedNetwork = as.matrix(read.table(smoothedNetFile, header=TRUE, row.names = 1))
        }

        # read positive set names
        possets = as.character(read.table(possetfile)[,1])
        #transmat = as.matrix(transmat[uniIDs, uniIDs])

        rtable = NULL

        for (queryfile in possets){

            posset = tail(unlist(strsplit(queryfile, "/")),1)
            posset = gsub(".txt","",posset)

            if(file.info(queryfile)$size == 0){show(c("Empty file ", queryfile)); next}

            # Stage 1
            # read query file
            query_gs = read.table(queryfile)
            if(dim(query_gs)[2]<2){
                query_gs = cbind(as.character(query_gs[,1]), rep(1, length(query_gs[,1])) )
            }
            rownames(query_gs) = as.character(query_gs[,1])
            queryIDs = sort(intersect(uniIDs, as.character(unique(query_gs[,1]))))
            nquery = length(queryIDs)

            set.seed(041185)
            folds = sample(cut(seq(1,nquery),breaks=nfolds,labels=FALSE))

            for(iter in 1:1){

                ## separate training and testing
                train_idxs = which(folds!=iter)
                test_idxs = which(folds==iter)
                train_nidxs = queryIDs[train_idxs]
                ntrain = length(train_nidxs)
                test_nidxs = queryIDs[test_idxs]
                ntest = length(test_nidxs)
                testuni = as.character(setdiff(uniIDs,train_nidxs))
                ntestuni = length(testuni)
                show(queryfile)

                # syntheticSamples = SMOTE(smoothedNetwork[train_nidxs, ], k=5, N=300)
                # smoothedNetwork = rbind(smoothedNetwork, syntheticSamples)
                # y_synthetic = structure(rep(1,nrow(syntheticSamples)), names = rownames(syntheticSamples))

                # run genlasso on training set and
                y = blankvec[uniIDs]
                midxs = match(train_nidxs, uniIDs)
                y[midxs] = as.numeric(query_gs[train_nidxs,2])
                #y = c(y, y_synthetic)
                glmmod = glmnet(x=smoothedNetwork, y=as.factor(y), alpha=1, family='binomial', nlambda=400)
                coefs = coef(glmmod)
                selectedFeatures = names(which(coefs[,max(which(abs(glmmod$df - nfeatures) == min(abs(glmmod$df[2:length(glmmod$df)] - nfeatures))))] != 0) )
                #selectedFeatures = names(which(coefs[,ncol(coefs)] != 0))
                selectedFeatures = selectedFeatures[2:length(selectedFeatures)]  # get rid of "intercept"
                numFeatures = length(selectedFeatures)
                midxs = match(selectedFeatures, propIDs)
                if (length(midxs) == 0 || is.na(midxs[1])) {
                    #smoothedNetwork = smoothedNetwork[uniIDs, uniIDs]
                    next
                }
                weights = c(nrow(smoothedNetwork) / length(train_nidxs), 1)
                names(weights) = c(1, 0)

                #svc = svm(x = transmat[, midxs], as.factor(y), kernel='radial', gamma = 1 / length(uniIDs), probability = TRUE, class.weights = weights)
                #pred = predict(svc,  transmat[, midxs], probability = TRUE)
                svc = svm(x = smoothedNetwork[, midxs], as.factor(y), kernel='radial', gamma = 1 / length(midxs), probability = TRUE, class.weights = weights)
                pred = predict(svc,  smoothedNetwork[, midxs], probability = TRUE)

                #smoothedNetwork = smoothedNetwork[uniIDs, uniIDs]

                y = blankvec[uniIDs]
                midxs = match(test_nidxs, uniIDs)
                y[midxs] = as.numeric(query_gs[test_nidxs,2])
                model = prediction(attr(pred, "probabilities")[,2][testuni], as.factor(y)[testuni])
                auc  = performance(model, "auc")
                perf = performance(model,"tpr","fpr")
                test_aucval = round(as.numeric(slot(auc, "y.values")),3)
                show(test_aucval)

                y = blankvec[uniIDs]
                midxs = match(train_nidxs, uniIDs)
                y[midxs] = as.numeric(query_gs[train_nidxs,2])
                trainuni = as.character(setdiff(uniIDs,test_nidxs))
                model = prediction(attr(pred, "probabilities")[,2][trainuni], as.factor(y)[trainuni])
                auc  = performance(model, "auc")
                perf = performance(model,"tpr","fpr")
                train_aucval = round(as.numeric(slot(auc, "y.values")),3)
                show(train_aucval)

                rtable = rbind(rtable, c(queryfile, numFeatures, test_aucval, train_aucval))
            } #end iter
        } #end queryset

        colnames(rtable) = c("Query File", "Number of Features", "Test AUC", "Train AUC")
        write.table(rtable, "out3.txt", quote=F, row.names=F, sep=",")
    } #end restart
} #end function


rwr_2stage(possetfile = "/workspace/R-code/data/msigdb/setlists/combo.txt", unifile = "/workspace/R-code/data/msigdb/gene_sets/hsap_universe.txt", networkfile = "/workspace/R-code/data/msigdb/networks/network-v3beta/gene_ontology_kegg_pathway.edge", weight = "unweighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 3, st2keep = 50, nfeatures = 400, property_types = c("gene_ontology", "kegg_pathway"), writepreds = 1, outdir = "/workspace/R-code/data/msigdb/results/1st_")