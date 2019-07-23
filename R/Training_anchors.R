###############################################################
# File name: Training_anchors.R
# This script is part of CASAR for cnv detection in MET, EGFR and ERBB2
# in 30K panel
# This script trains a model from autumated segmentations to defind
# a set of anchor segments for recalibrating coverage ratio
# This script is only suitable for 30K panel, to use this script in 
# another panel, one should first generate a segment file using CBS
#
# Author: Ni Shuai nishuai@yahoo.com
# BGI tumor R&D department
# Input: a segment file in bed format, a directory of a set of input 
# mpileup files (negative CNV samples) and output_directory (gwd by default)
#################################################################

####Getting arguments
args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input

setwd('E:/wrk/CASAR_cnv/')
args=c('E:/wrk/CASAR_cnv/inputdata/segmented_bed_30k.bed', 
       'E:/wrk/CASAR_cnv/0715_std_test/met2.1875ERBB2.625/',
       'E:/wrk/CASAR_cnv/temp/')

list.files(args[2])

####assigning character arguments to variables
segment_bed=input_mpileup_dir=output_dir=NULL
if (length(args)<2) {
  stop("Usage: Rscript Training_anchors.R segment_bed input_mpileup_dir output_dir", call.=FALSE)
} else if (length(args)>1) {
  # default output file
  segment_bed=args[1]
  input_mpileup_dir=args[2]  # Average bkgrd rate file
  if (length(args)>2) output_dir=args[3] else output_dir='.'  # 3nt bkgrd rate file
}


####processing input

###output dir exist
if(!dir.exists(output_dir)) stop('The output directory does not exist')

###checking segetment file
out = tryCatch({
  read.table(segment_bed)
},
error=function(cond) {
  message(paste("file destination does not seem to exist:", segment_bed))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the input bed file caused a warning:",
                segment_bed))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={segment_bed=read.table(segment_bed)}
)

###checking segetment file
out = tryCatch({
  list.files(input_mpileup_dir, pattern = '*freq')
},
error=function(cond) {
  message(paste("pileup file directory does not seem to exist:", input_mpileup_dir))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the input bed file caused a warning:",
                input_mpileup_dir))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={input_mpileup=list.files(input_mpileup_dir, pattern = '*freq')}
)

if (length(input_mpileup)<5) stop('At least 5 samples are required to determine
                                  anchor segments')

####sourcing scripts and libraries
library(MASS)
source('R/get_loci_in_bed.R')
mean_coorelations=list()

###reading mean depth for each segment in each sample
for (i in 1:length(input_mpileup)){
  file=read.table(file.path(input_mpileup_dir, input_mpileup[i]), 
                  header=TRUE)
  each_mean_depth=c()
  for (k in 1:nrow(segment_bed)){
    locis=get_loci_in_bed(file,segment_bed[k,])
    each_mean_depth=c(each_mean_depth, mean(locis$DEPTH, na.rm=TRUE))
  }
  mean_coorelations[[i]]=each_mean_depth
}

mean_coorelations=do.call(rbind, mean_coorelations)
# bb=mean_coorelations/rowMeans(mean_coorelations)
# plot(bb[1,c(124:128, 187)], type='l', ylim=c(0, 3), col='red', main='before CNV recalibration',
#      xlab='EGFR1, EGFR2, EGFR3, EGFR4, MET, HER2')
# for (i in 2:28){lines(bb[i,c(124:128, 187)], type='l')}
# for (i in c(12, 22:25)){lines(bb[i,c(124:128, 187)], type='l', col='red')}
# for (i in c(2:4, 26:28)){lines(bb[i,c(124:128, 187)], type='l', col='green')}
# for (i in c(5:10)){lines(bb[i,c(124:128, 187)], type='l', col='gold')}

###Met gene anchor
n_anchors=16
met_top_I=tail(order(cor(mean_coorelations)[128,]), n_anchors)
met_top_I=met_top_I[!is.na(cor(mean_coorelations)[128,][met_top_I])]
met_top_I=met_top_I[-length(met_top_I)]

met_top_I2=tail(order(cor(mean_coorelations)[129,]), n_anchors)
met_top_I2=met_top_I2[!is.na(cor(mean_coorelations)[129,][met_top_I2])]
met_top_I2=met_top_I2[-length(met_top_I2)]

met_top_I3=tail(order(cor(mean_coorelations)[130,]), n_anchors)
met_top_I3=met_top_I3[!is.na(cor(mean_coorelations)[130,][met_top_I3])]
met_top_I3=met_top_I3[-length(met_top_I3)]

met_top_I4=tail(order(cor(mean_coorelations)[131,]), n_anchors)
met_top_I4=met_top_I4[!is.na(cor(mean_coorelations)[131,][met_top_I4])]
met_top_I4=met_top_I4[-length(met_top_I4)]

met_top_I5=tail(order(cor(mean_coorelations)[132,]), 90)
met_top_I5=met_top_I5[!is.na(cor(mean_coorelations)[132,][met_top_I5])]
met_top_I5=met_top_I5[-length(met_top_I5)]

###ERBB2 anchor
erbb_top_I=tail(order(cor(mean_coorelations)[191,]), n_anchors)
erbb_top_I=erbb_top_I[!is.na(cor(mean_coorelations)[191,][erbb_top_I])]
erbb_top_I=erbb_top_I[-length(erbb_top_I)]

erbb_top_I2=tail(order(cor(mean_coorelations)[192,]), n_anchors)
erbb_top_I2=erbb_top_I2[!is.na(cor(mean_coorelations)[192,][erbb_top_I2])]
erbb_top_I2=erbb_top_I2[-length(erbb_top_I2)]

erbb_top_I3=tail(order(cor(mean_coorelations)[193,]), n_anchors)
erbb_top_I3=erbb_top_I3[!is.na(cor(mean_coorelations)[193,][erbb_top_I3])]
erbb_top_I3=erbb_top_I3[-length(erbb_top_I3)]
###EGFR anchor

egfr1_top_I=tail(order(cor(mean_coorelations)[124,]), n_anchors)
egfr1_top_I=egfr1_top_I[!is.na(cor(mean_coorelations)[124,][egfr1_top_I])]
egfr1_top_I=egfr1_top_I[!egfr1_top_I %in% c(124:127)]

egfr2_top_I=tail(order(cor(mean_coorelations)[125,]), n_anchors)
egfr2_top_I=egfr2_top_I[!is.na(cor(mean_coorelations)[125,][egfr2_top_I])]
egfr2_top_I=egfr2_top_I[!egfr2_top_I %in% c(124:127)]

egfr3_top_I=tail(order(cor(mean_coorelations)[126,]), 90)
egfr3_top_I=egfr3_top_I[!is.na(cor(mean_coorelations)[126,][egfr3_top_I])]
egfr3_top_I=egfr3_top_I[!egfr3_top_I %in% c(124:127)]

egfr4_top_I=tail(order(cor(mean_coorelations)[127,]), 90)
egfr4_top_I=egfr4_top_I[!is.na(cor(mean_coorelations)[127,][egfr4_top_I])]
egfr4_top_I=egfr4_top_I[!egfr4_top_I %in% c(124:127)]

####generating 6 models for 6 possible snv segments

##MET
met_normal=fitdistr(
  rowMeans(mean_coorelations[,met_top_I])/mean_coorelations[,128],'normal')
met_normal=met_normal$estimate
met_cauchy=fitdistr(
  rowMeans(mean_coorelations[,met_top_I])/mean_coorelations[,128],'cauchy')
met_cauchy=met_cauchy$estimate
met_model=list(met_normal, met_cauchy)

met2_normal=fitdistr(
  rowMeans(mean_coorelations[,met_top_I2])/mean_coorelations[,129],'normal')
met2_normal=met2_normal$estimate
met2_cauchy=fitdistr(
  rowMeans(mean_coorelations[,met_top_I2])/mean_coorelations[,129],'cauchy')
met2_cauchy=met2_cauchy$estimate
met2_model=list(met2_normal, met2_cauchy)

met3_normal=fitdistr(
  rowMeans(mean_coorelations[,met_top_I3])/mean_coorelations[,130],'normal')
met3_normal=met3_normal$estimate
met3_cauchy=fitdistr(
  rowMeans(mean_coorelations[,met_top_I3])/mean_coorelations[,130],'cauchy')
met3_cauchy=met3_cauchy$estimate
met3_model=list(met3_normal, met3_cauchy)

met4_normal=fitdistr(
  rowMeans(mean_coorelations[,met_top_I4])/mean_coorelations[,131],'normal')
met4_normal=met4_normal$estimate
met4_cauchy=fitdistr(
  rowMeans(mean_coorelations[,met_top_I4])/mean_coorelations[,131],'cauchy')
met4_cauchy=met4_cauchy$estimate
met4_model=list(met4_normal, met4_cauchy)

met5_normal=fitdistr(
  rowMeans(mean_coorelations[,met_top_I5])/mean_coorelations[,132],'normal')
met5_normal=met5_normal$estimate
met5_cauchy=fitdistr(
  rowMeans(mean_coorelations[,met_top_I5])/mean_coorelations[,132],'cauchy')
met5_cauchy=met5_cauchy$estimate
met5_model=list(met5_normal, met5_cauchy)

##ERBB2
erbb2_normal=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I])/mean_coorelations[,191],'normal')
erbb2_normal=erbb2_normal$estimate
erbb2_cauchy=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I])/mean_coorelations[,191],'cauchy')
erbb2_cauchy=erbb2_cauchy$estimate
erbb2_model=list(erbb2_normal, erbb2_cauchy)

erbb2_2_normal=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I2])/mean_coorelations[,192],'normal')
erbb2_2_normal=erbb2_2_normal$estimate
erbb2_2_cauchy=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I2])/mean_coorelations[,192],'cauchy')
erbb2_2_cauchy=erbb2_2_cauchy$estimate
erbb2_2_model=list(erbb2_2_normal, erbb2_2_cauchy)

erbb2_3_normal=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I3])/mean_coorelations[,193],'normal')
erbb2_3_normal=erbb2_3_normal$estimate
erbb2_3_cauchy=fitdistr(
  rowMeans(mean_coorelations[,erbb_top_I3])/mean_coorelations[,193],'cauchy')
erbb2_3_cauchy=erbb2_3_cauchy$estimate
erbb2_3_model=list(erbb2_3_normal, erbb2_3_cauchy)
##EGFR1234
egfr1_normal=fitdistr(
  rowMeans(mean_coorelations[,egfr1_top_I])/mean_coorelations[,124],'normal')
egfr1_normal=egfr1_normal$estimate
egfr1_cauchy=fitdistr(
  rowMeans(mean_coorelations[,egfr1_top_I])/mean_coorelations[,124],'cauchy')
egfr1_cauchy=egfr1_cauchy$estimate
egfr1_model=list(egfr1_normal, egfr1_cauchy)

egfr2_normal=fitdistr(
  rowMeans(mean_coorelations[,egfr2_top_I])/mean_coorelations[,125],'normal')
egfr2_normal=egfr2_normal$estimate
egfr2_cauchy=fitdistr(
  rowMeans(mean_coorelations[,egfr2_top_I])/mean_coorelations[,125],'cauchy')
egfr2_cauchy=egfr2_cauchy$estimate
egfr2_model=list(egfr2_normal, egfr2_cauchy)

egfr3_normal=fitdistr(
  rowMeans(mean_coorelations[,egfr3_top_I])/mean_coorelations[,126],'normal')
egfr3_normal=egfr3_normal$estimate
egfr3_cauchy=fitdistr(
  rowMeans(mean_coorelations[,egfr3_top_I])/mean_coorelations[,126],'cauchy')
egfr3_cauchy=egfr3_cauchy$estimate
egfr3_model=list(egfr3_normal, egfr3_cauchy)

egfr4_normal=fitdistr(
  rowMeans(mean_coorelations[,egfr4_top_I])/mean_coorelations[,127],'normal')
egfr4_normal=egfr4_normal$estimate
egfr4_cauchy=fitdistr(
  rowMeans(mean_coorelations[,egfr4_top_I])/mean_coorelations[,127],'cauchy')
egfr4_cauchy=egfr4_cauchy$estimate
egfr4_model=list(egfr4_normal, egfr4_cauchy)

models=list(met_model, met2_model, met3_model, met4_model, 
            met5_model,erbb2_model,erbb2_2_model,erbb2_3_model,egfr1_model, 
            egfr2_model, egfr3_model, egfr4_model)
names(models)=c('met_model1','met_model2','met_model3','met_model4','met_model5',
                'erbb2_model1','erbb2_model2','erbb2_model3','egfr1_model', 
              'egfr2_model','egfr3_model','egfr4_model')

saveRDS(models, file.path(output_dir, 'model_parameters.rda'))

anchors=list(met_top_I, met_top_I2, met_top_I3, met_top_I4, met_top_I5, 
             erbb_top_I, erbb_top_I2, erbb_top_I3, egfr1_top_I,
             egfr2_top_I, egfr3_top_I, egfr4_top_I)

names(anchors)=c('met1','met2','met3','met4','met5','erbb2_1','erbb2_2','erbb2_3','egfr1','egfr2','egfr3','egfr4')
saveRDS(anchors,file.path(output_dir, 'model_anchors.rda'))

