vers = R.Version()$major
# .libPaths("/shared-mounts/sinhas/softwares/R-3.1.3/library/")
# if(vers==2){
# 	.libPaths("/shared-mounts/sinhas/lib/R.2.15/")
# }
vers = paste(sep=".","R",R.Version()$major,R.Version()$minor)
show(vers)
library(Matrix)
library(ROCR)

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

possetfile = "C:/Users/Alex/workspace/data/msigdb/setlists/combo.txt"
unifile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/hsap_universe.txt"
networkfile = "C:/Users/Alex/workspace/data/msigdb/networks/1sp_many/Human.GO.TM.Loose.txt"
outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_"
weight = "weighted"
normalize = "type"
restarts = .7
maxiters = 50
thresh = 0.001

st2keep = 50
writepreds = 1
nfolds = 10

# queryfile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/msigdb_combo/ZHANG_BREAST_CANCER_PROGENITORS.CB.txt"
# restart = .7
# iter = 1
forcelarge = 0

property_types = c("go_curated_evidence", "go_inferred_evidence", "pfam_domain")

rwr_2stage<- function(possetfile = "C:/Users/Alex/workspace/data/msigdb/setlists/combo.txt", unifile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/hsap_universe.txt", networkfile = "C:/Users/Alex/workspace/data/msigdb/networks/1sp_many/Human.GO.TM.Loose.txt", weight = "weighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 10, st2keep = 50, property_types = c("go_curated_evidence", "go_inferred_evidence", "pfam_domain"), writepreds = 0, outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_"){

	uni = tail(unlist(strsplit(unifile, "/")),1)
	uni = gsub(".txt","",uni)
	network = tail(unlist(strsplit(networkfile, "/")),1)
	network = gsub(".edge","",network)
	network = gsub(".txt","",network)
	possetname = tail(unlist(strsplit(possetfile, "/")),1)
	possetname = gsub(".txt","",possetname)

	restable = NULL
	resfile = paste(sep="", outdir, uni, ".", network, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", st2keep, ".", paste(collapse="_", restarts), ".", possetname, ".stats")

	if(file.exists(resfile)){show(c("Already exists ", resfile)); return(1)}
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

	if(length(uniIDs)<1){ return(-1)}

	for (restart in restarts){
		# store results in table
		evaltabstart = cbind(nodenames,-1,0)
		colnames(evaltabstart) = c("node", "type", "universe")
		rownames(evaltabstart) = nodenames
		evaltabstart[uniIDs,"universe"] = 1
		
		# store node type info in table
		evaltabstart[as.character(features[,1]),"type"] = as.numeric(features[,2])
		
		# run baseline rwr network
		basefile = paste(sep="", outdir, uni, ".", network, ".", weight, ".", normalize, ".", maxiters, ".", thresh, ".", restart, ".base")
		show(basefile)
		
		blankvec = structure(rep(0,nnodes), names = nodenames)
		startvec = blankvec + 1 / nnodes
		if(!file.exists(basefile)){
			query = blankvec
			midxs = match(uniIDs, nodenames)
			query[midxs] = as.numeric(universe[uniIDs,2])
			
			rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)
			
			biter = rwr_res$iter
			evaltabstart = cbind(evaltabstart,as.numeric(rwr_res$vec))
			colnames(evaltabstart)[4] = "baseline"
			
			if(!file.exists(basefile)){
				write.table(evaltabstart, basefile, quote=F, sep="\t", row.names=T, col.names=NA)
			}
		} else { # read in base results into startvec
			evaltabstart = read.table(basefile)
		}
	
		startvec[nodenames] = as.numeric(evaltabstart[nodenames,"baseline"])  #starting distribution at baseline may speed convergence
		
		# read positive set names
		possets = as.character(read.table(possetfile)[,1])

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
			folds = sample(cut(seq(1,nquery),breaks=10,labels=FALSE))
			
			for(iter in 1:nfolds){
				# iter = 1
				## read query set
				outfile = paste(sep="", outdir, uni, ".", network, ".", weight, ".", normalize, ".",  maxiters, ".", thresh, ".", st2keep, ".", restart, ".", posset, ".", iter, ".", nfolds,".rwr")
				show(outfile)
				if(file.exists(outfile)){show(c("Already exists ", outfile)); next}

				## separate training and testing
				train_idxs = which(folds!=iter)
				test_idxs = which(folds==iter)
				train_nidxs = queryIDs[train_idxs]
				ntrain = length(train_nidxs)
				test_nidxs = queryIDs[test_idxs]
				ntest = length(test_nidxs)
				testuni = as.character(setdiff(uniIDs,train_nidxs))
				ntestuni = length(testuni)
				
				if(ntrain*ntest*1.0*ntestuni==0){ ## either no training or testing examples
					row1 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, 	iter, "baseline", -1000, -1000, -1000)
					row2 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage1", -1000, -1000, -1000)
					row3 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "diff", -1000, -1000, -1000)
					row4 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row1, row2, row3, row4)
					next
				}  
				
				# run RWR on training set
				query = blankvec
				midxs = match(train_nidxs, nodenames)
				query[midxs] = as.numeric(query_gs[train_nidxs,2])
				rwr_res = RWR(boolSparceMat, transmat, restart, query, startvec, maxiters, thresh)
			
				# store results
				qiter = rwr_res$iter
				evaltab = cbind(evaltabstart,0,0,as.numeric(rwr_res$vec[nodenames]))
				colnames(evaltab) = c("node", "type", "universe", "baseline", "train", "test", "stage1")
				diff = as.numeric(evaltab[,"stage1"]) - as.numeric(evaltab[,"baseline"])
				evaltab = cbind(evaltab,diff)
				evaltab[train_nidxs,"train"] = 1
				evaltab[test_nidxs,"test"] = 1

				# eval of baseline rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"baseline"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row	= c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "baseline", aucval, 0, length(uniIDs))
				show(row[11:14])
				restable = rbind(restable, row)

				# eval of query rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"stage1"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage1", aucval, qiter, ntrain)
				show(row[11:14])
				restable = rbind(restable, row)
				#plot(perf, lwd = 5)

				# eval of diff rwr on pos set
				model = prediction(as.numeric(evaltab[testuni,"diff"]), evaltab[testuni,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "diff", aucval, 0, ntrain)
				show(row[11:14])
				restable = rbind(restable, row)
				#plot(perf, lwd = 5)
				
				colnames(restable) = c("network", "weight", "normalize", "uni", "restart", "maxiters", "thresh", "st2keep", "posset", "nfolds", "iter", "stage1", "aucval", "niters", "ntrain")
				write.table(restable, resfile, quote=F, sep="\t", row.names=F)

				# stage 2
				# can skip in no feature nodes
				if(nfeats==0){
					row4 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row4)
					write.table(restable, resfile, quote=F, sep="\t", col.names=F, row.names=F)				
					if(iter == 1 && writepreds){
						write.table(evaltab, outfile, quote=F, sep="\t", row.names=T, col.names=NA)	
					}
					next;
				}
				
				# extract best features nodes
				keep = rep(0,nnodes)
				evaltab = cbind(evaltab, keep)
				ss = sort(as.numeric(evaltab[featnodes,"diff"]), decreasing=T, index.return = T)
				sortedfeats = featnodes[ss$ix]
				keepfeats = sortedfeats[1:nkeep]
				evaltab[keepfeats,"keep"] = 1

				# keep all non-property edges
				newedges = edges[which(edges$type %in% setdiff(all_etypes, prop_etypes)),]

				# add in edges connected to kept features
				keepidxs = which(edges[,"source"] %in% keepfeats)
				newedges = rbind(newedges, edges[keepidxs,])
				
				tmpnames = unique(c(as.character(newedges[,1]),as.character(newedges[,2])))
				train_nidxs2 = intersect(train_nidxs,tmpnames)
				test_nidxs2 = intersect(test_nidxs,tmpnames)
				testuni2 = intersect(testuni,tmpnames)

				# can skip if no overlap with the train or test set
				if(length(train_nidxs2)*length(test_nidxs2)*1.0*length(testuni2)==0){
					row4 = c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", -1000, -1000, -1000)
					restable = rbind(restable, row4)
					write.table(restable, resfile, quote=F, sep="\t", col.names=F, row.names=F)				
					if(iter == 1 && writepreds){
						write.table(evaltab, outfile, quote=F, sep="\t", row.names=T, col.names=NA)	
					}
					next;
				}
				
				# renormalized
				typetable2 = NULL
				if(normalize=="type"){
					typetable2 = aggregate(weight~type, data = newedges, FUN=sum)
					show(typetable2)
					rownames(typetable2) = as.character(typetable2[,1]	)
					newedges[,"weight"] = as.numeric(newedges[,"weight"]) / as.numeric(typetable2[as.character(newedges[,"type"]),"weight"])
					#show(newedges[1:5,])
				}

				transmat2 = NULL
				nodenames2 = NULL
				colsum2 = NULL
				# check is space matrix will work
				boolSparceMat2 = (length(tmpnames)^2 < .Machine$integer.max) * (1-forcelarge)  # will fit in sparce mat
				if(boolSparceMat2) { 
					# convert edgelist to matrix
					n2Matrix = threeCol2MaxMat( as.character(newedges[,"source"]),  as.character(newedges[,"target"]), as.numeric(newedges[,"weight"]), 1)
					nodenames2 = rownames(n2Matrix)

					# column normalize
					colsum2 = colSums(n2Matrix)
					transmat2 = t(t(n2Matrix)/colsum2)
					rm(n2Matrix)
				} else {  # must use personal list format
					ll2 = threeCol2listMat( as.character(newedges[,"source"]),  as.character(newedges[,"target"]), as.numeric(newedges[,"weight"]), 1)
					transmat2 = ll2
					nodenames2 = ll2$avals
					colsum2 = ll2$colsums
					rm(ll2)
				}
				nnodes2 = length(nodenames2)
				rm(newedges)

				# run RWR on s2 network on training set
				blankvec2 = structure(rep(0,nnodes2), names = nodenames2)
				query2 = blankvec2
				midxs = match(train_nidxs2, nodenames2)
				query2[midxs] = as.numeric(query_gs[train_nidxs2,2])
				startvec2 = blankvec2 + 1 / nnodes2
				rwr_res = RWR(boolSparceMat2, transmat2, restart, query2, startvec2, maxiters, thresh)

				# store results
				q2iter = rwr_res$iter
				stage2 = rep(0,nnodes)
				evaltab = cbind(evaltab, stage2)
				evaltab[nodenames2,"stage2"] =  rwr_res$vec[nodenames2]
				
				# eval of stage2 rwr on pos set
				model = prediction(as.numeric(evaltab[testuni2,"stage2"]), evaltab[testuni2,"test"])
				auc  = performance(model, "auc")
				perf = performance(model,"tpr","fpr")
				aucval = round(as.numeric(slot(auc, "y.values")),3)
				row	= c(network, weight, normalize, uni, restart, maxiters, thresh, st2keep, posset, nfolds, iter, "stage2", aucval, q2iter, length(train_nidxs2))
				show(row[11:14])
				restable = rbind(restable, row)
				write.table(restable, resfile, quote=F, sep="\t", row.names=F)

				if(iter == 1 && writepreds){
					write.table(evaltab, outfile, quote=F, sep="\t", row.names=T, col.names=NA)	
				}

			} #end iter
		} #end queryset
	} #end restart
} #end function


rwr_2stage(possetfile = "C:/Users/Alex/workspace/data/msigdb/setlists/combo.txt", unifile = "C:/Users/Alex/workspace/data/msigdb/gene_sets/hsap_universe.txt", networkfile = "C:/Users/Alex/workspace/data/msigdb/networks/1sp_many/Human.GO.TM.Loose.txt", weight = "weighted", normalize = "type", restarts = .7, maxiters = 50, thresh = 0.001, nfolds = 10, st2keep = 50, property_types = c("go_curated_evidence", "go_inferred_evidence", "STRING_textmining"), writepreds = 1, outdir = "C:/Users/Alex/workspace/data/msigdb/results/1st_")


