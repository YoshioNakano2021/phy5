#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import glob
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.conversion import py2rpy,rpy2py
from rpy2.robjects import pandas2ri
from rpy2.robjects import numpy2ri

parser = argparse.ArgumentParser()

parser.add_argument("-d", "-diredtory", type=str)
parser.add_argument("-f", "-file", type=str)
parser.add_argument("-s", "-distance", type=str)
parser.add_argument("-g", "-agglomeration", type=str)

args = parser.parse_args()
print(args.d)
print(args.f)
print(args.s)
print(args.g)

if args.d is None and args.f is None:
    print('Fasta file or a folder containing fasta files is required.')
    sys.exit()

if args.d and args.f:
    print('Both file and folder are specified.  The fasta file will be used.')

distmtx = ["euclidean","manhattan", "maximum", "canberra", "binary"]
agglmtd = ["ward", "average", "single", "complete", "mcquitty", "median", "centroid"]
if args.s != None:
    distance_matrix = args.s
    if distance_matrix not in distmtx:
        print("The distance matrix is unknown. It should be chosen one out of 'euclidean','manhattan', 'maximum', 'canberra', and 'binary'.")
        distance_matrix = "manhattan"
else:
    distance_matrix = "manhattan"
    print("A distance matrix must be selected.  Manhattan distance will be used as a distance matrix.")
if args.g != None:
    aggl_method = args.g
    if aggl_method not in agglmtd:
        print("The agglomeration method is unknown. It should be chosen one out of 'ward', 'average', 'single', 'complete', 'mcquitty', 'median', and 'centroid'.")
        aggl_method = "ward.D2"
else:
    aggl_method = "ward.D2"
    print("An agglomeration method must be selected.  Ward.D2 will be used as an agglomeration method.")
if aggl_method == "ward":
    aggl_method = "ward.D2"
 
ro.r.assign("dm", distance_matrix)
ro.r.assign("am", aggl_method)

if args.f != None:
    is_file = os.path.isfile(args.f)
    if is_file:
        bstrings = importr("Biostrings")
        ro.r.assign("fl",args.f)
        ro.r('''
        dnaseq <- readDNAStringSet(fl)
        out5 <- oligonucleotideFrequency(dnaseq, width=5)
        tmp <- cbind(names(out5), out5)
        tmp.df <- data.frame(tmp)
        pentamers <- colnames(out5)
        rlist <- c()
        for (pntm in pentamers){
          nuc5 <- DNAString(pntm)
          cmp5 <- reverseComplement(nuc5)
          rpntm <- as.character(cmp5)
          if(!pntm %in% rlist){
            tmp.df[,pntm] <- tmp.df[,pntm] + tmp.df[,rpntm]
            tmp.df[,rpntm] <- NULL
            rlist[[(length(rlist) + 1)]] <- rpntm
          }
        }
        ''')
        res5 = ro.r('tmp.df')
        colnm = ro.r("colnames(tmp.df)")
        genenm = ro.r("names(dnaseq)")
        data_p = pd.DataFrame(res5)
        data_df = data_p.T
        data_df.columns = colnm 
        data_df.index = genenm
        print(data_df)
    else:
        print("The specified file does not exist.")
        sys.exit()


elif args.d != None:
    is_dir = os.path.isdir(args.d)
    if is_dir:
        path = args.d + "/*.fna"
        fastas = glob.glob(path)
        bstrings = importr("Biostrings")
        
        seqnames = []
        n = 1
        for fasta in fastas:
            with open(fasta) as f:
                seqname = f.readlines()[0]
            seqname = seqname.replace(' ','_')
            seqname = seqname.replace(',','_')
            seqname = seqname.replace('>','')
            seqname = seqname.replace('\n','')
            seqnames.append(seqname)
            ro.r.assign("fl", fasta)
            ro.r('''
            dnaseq <- readDNAStringSet(fl)
            out5 <- oligonucleotideFrequency(dnaseq, width=5)
            tmp <- cbind(names(out5), out5)
            tmp.df <- data.frame(tmp)
            pentamers <- colnames(out5)
            rlist <- c()
            for (pntm in pentamers){
              nuc5 <- DNAString(pntm)
              cmp5 <- reverseComplement(nuc5)
              rpntm <- as.character(cmp5)
              if(!pntm %in% rlist){
                tmp.df[,pntm] <- tmp.df[,pntm] + tmp.df[,rpntm]
                tmp.df[,rpntm] <- NULL
                rlist[[(length(rlist) + 1)]] <- rpntm
              }
            }
            sums5 <- colSums(tmp.df)
            ''')
            o5 = ro.r("sums5")
            onames = ro.r("names(sums5)")
            if n == 1:
                o5table = o5
                n += 1
            else:
                o5table = np.row_stack((o5table,o5))
                n += 1
        
        data_df = pd.DataFrame(o5table)
        data_df = data_df.astype("int")
        print(data_df)
        data_df.columns = onames
        data_df.index = seqnames
    else:
        print("The specified directory does not exist.")
        sys.exit()

data_df.to_csv('freq5.csv')
prop_df = data_df.apply(lambda x:x/sum(x), axis=1)
print(prop_df)
prop_df.to_csv('prop5.csv')
    
with localconverter(ro.default_converter + pandas2ri.converter):
   rprop_df = py2rpy(prop_df)
    
ro.r.assign("rprop_df", rprop_df)

ape = importr("ape")
pvclust = importr("pvclust")
ctc = importr("ctc")

ro.r('''
nr <- nrow(rprop_df)
if (nr<10) {
     plot_height <- 300
}else{
     plot_height <- nr*25
}
dist <- dist(rprop_df, method=dm)
dist.m <- as.matrix(dist)
dist_df <- data.frame(dist.m)
write.csv(dist_df, "DistTable.csv")
cluster <- hclust(dist, method=am)
cluster.nwk <- hc2Newick(cluster)
write(cluster.nwk, file = "Tree.nwk")
tcsv <- t(rprop_df)
btsv <- pvclust(tcsv, method.dist="manhattan",nboot=1000, method.hclust="ward.D2",r=seq(.5,1,by=.1))
boot <- round((btsv$edges$bp)*100)
    
pdfname <- paste("tree5_",dm,"-",am,".pdf", sep="")
pdfname <- sub("ward.D2","Ward", pdfname)
pdf(pdfname, height=(nr/2.5), width=10, pointsize=10)
plot(as.phylo(cluster))
nodelabels(boot,adj=c(1.2, -0.2), frame="n", cex=0.7)
dev.off()
''')
