###############################################################
# File name: cesar.R
# CESAR (CNV Estimation with Segmentation and Anchor Recalibration)
# Author: Ni Shuai nishuai@yahoo.com
# Tumor R&D department @ Beijing Genomics Institute
# Input: background_coverage_file, input_coverage_file, 
# segmentation bed file and output_directory (gwd by default)
#################################################################
args = commandArgs(trailingOnly=TRUE)  # 1-3 arguments from input
# setwd('E:/wrk/cesar_508/')
# args=c(paste0('E:/wrk/cesar_508/samples/4-1//LoD2_LOD2-1-4GW-OGTM001-6-6_LOD2-1-4GW-OGTM001-6_HD752-13.case.reg.depth'),
#        'temp/model_anchors.rda', 
#        'temp/model_parameters.rda',
#        'temp/')

####################################################################
# processing the input
#####################
cat('Cesar is a software for efficient calling gene-level CNVs from NGS sequencing data')
cat('\n')
cat('\n')

input_pileup=model_anchor=model_parameters=output_dir=NULL
if (length(args)>5) {
  stop("Cesar Usage: Rscript cesar.R input_mpileup model_anchors model_parameters bed_file output_dir", call.=FALSE)
}

if (length(args)<4) {
  stop("Cesar Usage: Rscript cesar.R input_mpileup model_anchors model_parameters bed_file output_dir", call.=FALSE)
} else if (length(args)>3) {
    # default output file
  input_pileup=args[1]
  model_anchors=args[2]
  model_parameters=args[3]  # Average bkgrd rate file
  bed_file=args[4]  # Average bkgrd rate file
  if (length(args)>4) output_dir=args[5] else output_dir=dirname(input_pileup)  # 3nt bkgrd rate file 
}
  
pileup_filename=input_pileup
####################################################################
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

####################################################################
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

####################################################################
###check the validity of input model_parameters
out = tryCatch({
  readRDS(model_parameters)
},
error=function(cond) {
  message(paste("file destination does not seem to exist:", model_parameters))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the model_parameter file caused a warning:",
                model_parameters))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={model_parameters=readRDS(model_parameters)}
)

####################################################################
###output dir existance
if(!dir.exists(output_dir)) stop('The output directory does not exist')

####################################################################
###check the validity of input bed file
out = tryCatch({
  read.table(bed_file)
},
error=function(cond) {
  message(paste("file destination does not seem to exist:", bed_file))
  message("Here's the original error message:")
  message(cond)
  # Choose a return value in case of error
  return(NA)
},
warning=function(cond) {
  message(paste("reading of the input bed file caused a warning:",
                bed_file))
  message("Here's the original warning message:")
  message(cond)
  # Choose a return value in case of warning
  return(NULL)
},
finally={segment_bed=read.table(bed_file)}
)


cat('Cesar: Sourcing scripts and libraries')
cat('\n')
####################################################################
####sourcing scripts and libraries
source('R/bed_depth_in_pileup.R')
locis=bed_depth_in_pileup(input_pileup,segment_bed)

cat('Cesar: Computing CNVs')
cat('\n')
####################################################################
####for each gene get is anchor and anchor ratio
p_values=c()
gain_loss=c()
anchor_ratio=c()
for (i in 1:length(locis)){
  ##get anchor relative ratio
  anchor_ratio=c(anchor_ratio, mean(locis[model_anchors[[i]]])/locis[i])
  
  gain_loss=c(gain_loss, model_parameters[[i]][1]/anchor_ratio[i])
  ##test anchor relative ratio
  ###tried cauchy but not as good as normal
  #p_values=c(p_values, dcauchy(anchor_ratio[i], model_parameters[[i]][1], model_parameters[[i]][2]))
  p_values=c(p_values, dnorm(anchor_ratio[i], model_parameters[[i]][1], model_parameters[[i]][2]))
}

####################################################################
###output results
pdf(file.path(output_dir, 
  paste0(basename(pileup_filename), '_CNV_plot.pdf')),
  width = 16, height = 8)
p_values=p_values+1e-100
K=1:length(p_values) 
aa=loess(-log10(p_values)~ K, span=.01)
plot(predict(aa), type='l', main=pileup_filename)
abline(v=grep('MET',segment_bed$V4)[1], col='green', lty=2)
abline(v=grep('ERBB2',segment_bed$V4)[1], col='green', lty=2)
dev.off() -> msg_supressor

####################################################################
###test result
genes=strsplit(as.character(segment_bed$V4), '\\|NM_|\\|EN')
genes=unlist(lapply(genes, function(x) x[[1]]))
df=data.frame(genes, gain_loss, confidence=-log10(p_values))
result=aggregate(df[,2:3], by=list(df$genes), FUN=mean, na.rm=TRUE)
result=result[order(result$confidence, decreasing = TRUE), ]
names(result)=c('GeneName', 'CNV', 'Confidence')

write.table(head(result, 10), 
            file.path(output_dir, 
                      paste0(basename(pileup_filename), '_CNV_result.txt')),
            quote = FALSE,
            row.names = FALSE,
            sep = '\t')
cat('Cesar: CNV results successfully stored')
cat('\n')
