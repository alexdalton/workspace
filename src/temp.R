smoothedNetworkFile = "/workspace/R-code/Human.GO.TM.H.Loose.Smoothed.rds"
smoothedNetwork = readRDS(smoothedNetworkFile)

geneNames = rownames(smoothedNetwork)

for (geneName in geneNames) {
    sortedFeatureVector = sort(smoothedNetwork[geneName, ], decreasing = TRUE)
    zeroColumnNames = names(sortedFeatureVector[floor(0.2 * length(sortedFeatureVector)):length(sortedFeatureVector)])
    smoothedNetwork[geneName, zeroColumnNames] = 0.0
}

saveRDS(smoothedNetwork, "/workspace/R-code/Human.GO.TM.H.Loose.Smoothed.Zeroed80.rds")