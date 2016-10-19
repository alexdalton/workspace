#!/bin/sh

./embed -node ../data/node_GO_TM_PPI.txt -link ../data/link_GO_TM_PPI.txt -output GO_TM_PPI_DCA_2000.emb -binary 0 -size 2000 -negative 5 -samples 1 -iters 300 -threads 12 -model 2 -depth 5 -restart 0.7