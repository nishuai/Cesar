###############################################################
# File name: cnv30K.R
# detect copy number alteration from background coverage profiles
# in normal samples using Hidden Markov Model
# Author: Ni Shuai nishuai@yahoo.com
# Inspired by: Kenneth Jordan Mccallum and Ji-Ping Wang
# BGI tumor R&D department
# Input: background_coverage_file, input_coverage_file, 
# and output_directory (gwd by default)
#################################################################
setwd('E:/wrk/cnv_ctdna/')
loci_nromal_samples='inputdata/normal_samples/normal2'
loci_tumor_std='inputdata/cnv_samples/'
for (i in 1:length(list.files(loci_nromal_samples, pattern = '*freq'))){
 args=c('master_depth_std.rda', 
        file.path(loci_nromal_samples,  
                  list.files(loci_nromal_samples)[i]),
        'temp/')
 ttt(args)
}

# args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input
# processing the input
#####################
ttt=function(args){
master_depth=input_depth=output_dir=NULL
if (length(args)<2) {
  stop("Usage: Rscript cnv30k.R bg_depth_file input_mpileup output_dir", call.=FALSE)
} else if (length(args)>1) {
  # default output file
  master_depth=args[1]
  input_depth=args[2]  # Average bkgrd rate file
  if (length(args)>2) output_dir=args[3]  # 3nt bkgrd rate file
}

###check the validity of input files
out = tryCatch({
  readRDS(master_depth)
  },
  error=function(cond) {
  message(paste("file destination does not seem to exist:", master_depth))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
  },
  warning=function(cond) {
  message(paste("reading of the input background coverage file caused a warning:",
                master_depth))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
  },
  finally={master_depth=readRDS(master_depth)}
)

####read the input mpileup file
out = tryCatch({
  read.table(input_depth, header=TRUE)
  },
  error=function(cond) {
  message(paste("file destination does not seem to exist:", master_depth))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
  },
  warning=function(cond) {
  message(paste("reading of the input mpileup file caused a warning:",
                master_depth))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
  },
  finally={input_depth=read.table(input_depth, header=TRUE)}
) 


###Get target region
source('R/get_loci_in_bed.R')
ref_bed=read.table('ref/RB_lung_probe.region.original.bed')
master_depth=get_loci_in_bed(master_depth, ref_bed[44,])
input_depth=get_loci_in_bed(input_depth, ref_bed[44,])
master_depth=master_depth[complete.cases(master_depth),]
input_depth=input_depth[match( master_depth$POSITION, input_depth$POSITION),]
print(dim(master_depth));print(dim(input_depth))
####panel segmentation using coverage behavior
library(PSCBS)
segs=segmentByCBS(master_depth$DEPTH,  
                          min.width=5, verbose=0)

segs=segs$output
###remove uncessary segments
####injacent segments must differ by at least 1.15 times in mean
mergerable=abs(segs$mean[2:length(segs$mean)]/
                 segs$mean[2:length(segs$mean)-1]-1)>0.1

####injacent segments must differ by at least 500 depth in mean
mergerable=mergerable & abs(segs$mean[2:length(segs$mean)]-
                 segs$mean[2:length(segs$mean)-1]-1)>500

####segment length must span at 6 nucleotide
mergerable=c(mergerable, TRUE)

segs=segs[mergerable,]
segments=list()
segments[[1]]= 1:segs$end[1]
for (i in 2:length(segs$end)){
  segments[[i]]=ceiling(segs$end[i-1]):floor(segs$end[i])
}

####do statistics on each segment
bg_seg_depth=sapply(segments, function(I) mean(master_depth$DEPTH[I]))
bg_seg_depth= bg_seg_depth/mean(bg_seg_depth, na.rm = TRUE)
in_seg_depth=sapply(segments, function(I) mean(input_depth$DEPTH[I]))
in_seg_depth= in_seg_depth/mean(in_seg_depth, na.rm = TRUE)
# plot(log2(bg_seg_depth/in_seg_depth), type='l')
plot(log2(in_seg_depth/bg_seg_depth), type='l')
abline(h=0, col='red')

library(HMM)
hmm=initHMM(States=c('plus','minus','stable'),
            Symbols=c('P','M','S'),startProbs = c(0.2,0.2,0.6),
            transProbs=matrix(c(0.9998, 0.0001, 0.0001, 
                                0.0001, 0.9998, 0.0001, 
                                0.005, 0.005, 0.99),3, byrow=TRUE),
            emissionProbs=matrix(c(0.6, 0.05, 0.35, 0.05, 0.6, 0.35, 0.2, 0.2, 0.6),3,
                                 byrow=TRUE))

####estimating CNV using HMM with sample-specefic variance
log2FC=log2(in_seg_depth/bg_seg_depth)
standard_error=sd(log2FC, na.rm = TRUE)
print(standard_error)
log2FC[is.na(log2FC)]=0
P_seq=rep('S', length(log2FC))
P_seq[log2FC>standard_error]='P'
P_seq[log2FC< -standard_error]='M'
print(sum(P_seq=='M'))
print(sum(P_seq=='P'))
post <- posterior(hmm, P_seq)
lines(c(2,1,1.5)[apply(post, 2, function(x) which.max(x))])

# ####adding SNP information
# vafs=(input_depth[,c(7,9,11,13)] + input_depth[,c(7,9,11,13)+1])/input_depth$DEPTH
# vafs=apply(vafs, 1, function(x) max(x))
# vafs[is.na(vafs)]=0
# vafs_seg=sapply(segments, function(I) {
#   (length(vafs[I][abs(vafs[I]-0.5)<0.45])+3)/
#     (length(vafs[I][abs(vafs[I]-0.5)<=0.05 & abs(vafs[I]-1)<=0.05 ])+3)
# })
# lines(vafs_seg, col='red')
# abline(h=0, col='green')
}
