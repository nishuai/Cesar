###############################################################
# File name: cesar.R
# CESAR (CNV Estimation with Segmentation and Anchor Recalibration)
# Author: Ni Shuai nishuai@yahoo.com
# Tumor R&D department @ Beijing Genomics Institute
# Input: background_coverage_file, input_coverage_file, 
# and output_directory (gwd by default)
#################################################################
args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input
#setwd('E:/wrk/CASAR_cnv/')
#args=c(paste0('E:/wrk/CASAR_cnv/0715_test/normal/', kk), 
#       'temp/model_anchors.rda', 
#       'temp/model_parameters.rda',
#       'temp/')
# processing the input
#####################
input_pileup=model_anchor=model_parameters=output_dir=NULL
if (length(args)>4) {
  stop("Usage: Rscript casar.R input_mpileup model_anchors model_parameters output_dir", call.=FALSE)
}

if (length(args)<3) {
  stop("Usage: Rscript casar.R input_mpileup model_anchors model_parameters output_dir", call.=FALSE)
} else if (length(args)>2) {
    # default output file
  input_pileup=args[1]
  model_anchors=args[2]
  model_parameters=args[3]  # Average bkgrd rate file
  if (length(args)>3) output_dir=args[4] else output_dir='.'  # 3nt bkgrd rate file 
}
  
###checking input pileup file
out = tryCatch({
  read.table(input_pileup, header = TRUE)},
error=function(cond) {
  message(paste("file destination does not seem to exist:", input_pileup))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)},
warning=function(cond) {
  message(paste("reading of the input bed file caused a warning:",
                input_pileup))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)},
finally={input_pileup=read.table(input_pileup, header = TRUE)}
)

###check the validity of input model_anchor
out = tryCatch({
  readRDS(model_anchors)
},
error=function(cond) {
  message(paste("file destination does not seem to exist:", model_anchors))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the input background coverage file caused a warning:",
                model_anchors))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={model_anchors=readRDS(model_anchors)}
)

###check the validity of input model_parameters
out = tryCatch({
  readRDS(model_parameters)
},
error=function(cond) {
  message(paste("file destination does not seem to exist:", odel_parameters))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the input background coverage file caused a warning:",
                model_parameters))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={model_parameters=readRDS(model_parameters)}
)

###output dir existance
if(!dir.exists(output_dir)) stop('The output directory does not exist')

####sourcing scripts and libraries
source('R/get_loci_in_bed.R')
segment_bed=read.table('inputdata/segmented_bed_30k.bed')
each_mean_depth=c()
for (k in 1:nrow(segment_bed)){
  locis=get_loci_in_bed(input_pileup,segment_bed[k,])
  each_mean_depth=c(each_mean_depth, mean(locis$DEPTH, na.rm=TRUE))
}
####for each gene get is anchor and anchor ratio
cnv_seg_depth=each_mean_depth[c(128:132, 191:193, 124:127)]
anchor_depth=lapply(model_anchors, function(x) 
  mean(each_mean_depth[x], na.rm = TRUE))
anchor_ratio=as.numeric(anchor_depth)/cnv_seg_depth 
met_cnv1=dnorm(anchor_ratio[1],model_parameters$met_model1[[1]]['mean'], 
      model_parameters$met_model1[[1]]['sd'])
met_gain1=anchor_ratio[1]/model_parameters$met_model1[[1]]['mean'] 

met_cnv2=dnorm(anchor_ratio[2],model_parameters$met_model2[[1]]['mean'], 
      model_parameters$met_model2[[1]]['sd'])
met_gain2=anchor_ratio[2]/model_parameters$met_model2[[1]]['mean'] 

met_cnv3=dnorm(anchor_ratio[3],model_parameters$met_model3[[1]]['mean'], 
      model_parameters$met_model3[[1]]['sd'])
met_gain3=anchor_ratio[3]/model_parameters$met_model3[[1]]['mean'] 

met_cnv4=dnorm(anchor_ratio[4],model_parameters$met_model4[[1]]['mean'], 
      model_parameters$met_model4[[1]]['sd'])
met_gain4=anchor_ratio[4]/model_parameters$met_model4[[1]]['mean'] 

met_cnv5=dnorm(anchor_ratio[5],model_parameters$met_model5[[1]]['mean'], 
      model_parameters$met_model5[[1]]['sd'])
met_gain5=anchor_ratio[5]/model_parameters$met_model5[[1]]['mean'] 

erbb2_cnv1=dnorm(anchor_ratio[6],model_parameters$erbb2_model1[[1]]['mean'], 
      model_parameters$erbb2_model1[[1]]['sd'])
erbb2_gain1=anchor_ratio[6]/model_parameters$erbb2_model1[[1]]['mean'] 

erbb2_cnv2=dnorm(anchor_ratio[7],model_parameters$erbb2_model2[[1]]['mean'], 
      model_parameters$erbb2_model2[[1]]['sd'])
erbb2_gain2=anchor_ratio[7]/model_parameters$erbb2_model2[[1]]['mean'] 

erbb2_cnv3=dnorm(anchor_ratio[8],model_parameters$erbb2_model3[[1]]['mean'], 
      model_parameters$erbb2_model3[[1]]['sd'])
erbb2_gain3=anchor_ratio[8]/model_parameters$erbb2_model3[[1]]['mean'] 

egfr1_cnv=dnorm(anchor_ratio[9],model_parameters$egfr1[[1]]['mean'], 
      model_parameters$egfr1[[1]]['sd'])
egfr1_gain=anchor_ratio[9]/model_parameters$egfr1_model[[1]]['mean'] 

egfr2_cnv=dnorm(anchor_ratio[10],model_parameters$egfr2[[1]]['mean'], 
      model_parameters$egfr2[[1]]['sd'])
egfr2_gain=anchor_ratio[10]/model_parameters$egfr2_model[[1]]['mean'] 

egfr3_cnv=dnorm(anchor_ratio[11],model_parameters$egfr3[[1]]['mean'], 
      model_parameters$egfr3[[1]]['sd'])
egfr3_gain=anchor_ratio[11]/model_parameters$egfr3_model[[1]]['mean'] 

egfr4_cnv=dnorm(anchor_ratio[12],model_parameters$egfr4[[1]]['mean'], 
      model_parameters$egfr4[[1]]['sd'])
egfr4_gain=anchor_ratio[12]/model_parameters$egfr4[[1]]['mean'] 

cnv_overview=c(met_cnv1,met_cnv2, met_cnv3, met_cnv4, met_cnv5, 
               erbb2_cnv1, erbb2_cnv2, erbb2_cnv3, egfr1_cnv,
               egfr2_cnv, egfr3_cnv, egfr4_cnv)
cnv_overview=-log10(cnv_overview)
cnv_gainorloss=c(met_gain1, met_gain2, met_gain3, met_gain4, met_gain5, 
                 erbb2_gain1, erbb2_gain2, erbb2_gain3, egfr1_gain,
               egfr2_gain, egfr3_gain, egfr4_gain)
cnv_gainorloss=1/cnv_gainorloss
cnv_gainorloss=round(cnv_gainorloss, 3)
names(cnv_overview)=cnv_gainorloss
####CNV detection in MET
###if other 3 segments says aye or no opinion, and 
###have the same cnv status(gain or loss), then calculate
###the average cnv, vioce form the last segment can be ignored 
###because its not always accurate 
###defining significance level
Sig_level=1.6

opinion_leader=cnv_overview[1]
opinions=ifelse(cnv_overview[1:5]>Sig_level, 1, 0)*ifelse(names(cnv_overview)[1:5]>=1, 1, -1)
opinion_agreed=max(opinions[1:4])-min(opinions[1:4])!=2
###only 5 all agreed or the first 4 agreed, allowing for small perturbations
cnv_agreed=length(unique(as.numeric(names(cnv_overview))[1:5]>=1))==1 | 
  length(unique(names(cnv_overview)[1:4]>=1.1))==1 | 
  length(unique(names(cnv_overview)[1:4]>=0.9))==1 
  
MET_status=NA
MET_status=mean(as.numeric(names(cnv_overview))[1:4])
if (cnv_agreed & opinion_agreed){
  names(MET_status)=max(cnv_overview[1:4])
} else names(MET_status)=0.0

print(cnv_overview)
####CNV detection in ERBB2
opinions=ifelse(cnv_overview[6:8]>Sig_level, 1, 0)*ifelse(names(cnv_overview)[6:8]>=1, 1, -1)
print(opinions)
opinion_agreed=max(opinions)-min(opinions)!=2
###only all agreed
cnv_agreed=length(unique(as.numeric(names(cnv_overview))[6:8]>=1))==1 
ERBB2_status=NA
ERBB2_status=mean(as.numeric(names(cnv_overview))[6:8])

if (cnv_agreed & opinion_agreed){
  names(ERBB2_status)=max(as.numeric(cnv_overview[6:8]))
} else names(ERBB2_status)=0.0

####CNV detection in EGFR
opinion_leader=cnv_overview[10]
opinions=ifelse(cnv_overview[9:12]>Sig_level, 1, 0)*ifelse(names(cnv_overview)[9:12]>=1, 1, -1)
opinion_agreed=max(opinions)-min(opinions)!=2
###only all agreed
cnv_agreed=length(unique(as.numeric(names(cnv_overview))[9:12]>=1))==1 
EGFR_status=NA
EGFR_status=mean(as.numeric(names(cnv_overview))[9:12])
if (cnv_agreed & opinion_agreed){
  names(EGFR_status)=max(cnv_overview[9:12])
} else names(EGFR_status)=0.0

report=data.frame(Gene=c('Significant', 'Pvalue','CNV_status', 'CNV_level'), 
                  EGFR=c(ifelse(names(EGFR_status)>2, 'Yes','No'),
                  ifelse(as.numeric(names(EGFR_status))>Sig_level, round(as.numeric(names(EGFR_status)), 3), 0),
                  ifelse(as.numeric(names(EGFR_status))>Sig_level, ifelse(EGFR_status>=1, 'Gain', 'Loss'), 'No_CNV'),
                  ifelse(as.numeric(names(EGFR_status))>Sig_level, round(EGFR_status, 3), 'No_CNV')),
                  ERBB2=c(ifelse(as.numeric(names(ERBB2_status))>Sig_level, 'Yes','No'),
                  ifelse(as.numeric(names(ERBB2_status))>Sig_level, round(as.numeric(names(ERBB2_status)), 3), 0),
                  ifelse(as.numeric(names(ERBB2_status))>Sig_level, ifelse(ERBB2_status>=1, 'Gain', 'Loss'), 'No_CNV'),
                  ifelse(as.numeric(names(ERBB2_status))>Sig_level, round(ERBB2_status, 3), 'No_CNV')),
                  MET=c(ifelse(as.numeric(names(MET_status))>Sig_level, 'Yes','No'),
                  ifelse(as.numeric(names(MET_status))>Sig_level, round(as.numeric(names(MET_status)), 3), 0),
                  ifelse(as.numeric(names(MET_status))>Sig_level, ifelse(MET_status>=1, 'Gain', 'Loss'), 'No_CNV'),
                  ifelse(as.numeric(names(MET_status))>Sig_level, round(MET_status, 3), 'No_CNV')))

print(report)

write.table(report, file.path(output_dir, paste0(args[1],'.cnv')), quote=FALSE, row.names=FALSE)

