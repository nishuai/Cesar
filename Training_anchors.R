  ###############################################################
  # File name: Training_anchors.R
  # This script is part of CESAR for cnv detection in MET, EGFR and ERBB2
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
  
  # setwd('E:/wrk/cesar_508//')
  # args=c('E:/wrk/cesar_508/revised_cnv_anno_508.bed', 
  #        'E:/wrk/cesar_508/samples/m002_hd734/',
  #        'E:/wrk/cesar_508/temp/')
  
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
  
  cat('*** Cesar is learning CNV status for each segment in a group of control samples *** \n')
  cat('\n')
  
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
  finally={input_mpileup=list.files(input_mpileup_dir, pattern = '*depth$')}
  )
  
  if (length(input_mpileup)<5) stop('At least 5 samples are required to determine
                                    anchor segments')
  
  ####sourcing scripts and libraries
  library(MASS)
  library(data.table)
  source('R/bed_depth_in_pileup.R')
  
  mean_coorelations=list()
  ###reading mean depth for each segment in each sample
  cat('Calculating covariance matrix')
  cat('\n')
  ###reading mean depth for each segment in each sample
  for (i in 1:length(input_mpileup)){
    mfile=as.data.frame(fread(file.path(input_mpileup_dir, input_mpileup[i]), 
                        sep='\t', header=TRUE))
    locis=bed_depth_in_pileup(mfile,segment_bed)
    
    mean_coorelations[[i]]=locis
  }
  
  ###remove bed regions with extremely low coverage
  ##  average min_depth in each region to be 200
  mean_coorelations=do.call(rbind, mean_coorelations)
  mean_coorelations[,colMeans(mean_coorelations)<200]=0
  #################################################
  ####temporarily added for testing
  I_met=grep('MET', segment_bed[,4])
  mean_coorelations[grep('1-8GW-OGTM001', input_mpileup), I_met]=
    mean_coorelations[grep('1-8GW-OGTM001', input_mpileup), I_met]/1.09375
  I_ERBB2=grep('ERBB2', segment_bed[,4])
  mean_coorelations[grep('1-8GW-OGTM001', input_mpileup), I_ERBB2]=
    mean_coorelations[grep('1-8GW-OGTM001', input_mpileup), I_ERBB2]/1.31
  ##################################################
  ####temporarily added for testing
  I_met=grep('MET', segment_bed[,4])
  mean_coorelations[grep('1-4GW-OGTM001', input_mpileup), I_met]=
    mean_coorelations[grep('1-4GW-OGTM001', input_mpileup), I_met]/1.1875
  I_ERBB2=grep('ERBB2', segment_bed[,4])
  mean_coorelations[grep('1-4GW-OGTM001', input_mpileup), I_ERBB2]=
    mean_coorelations[grep('1-4GW-OGTM001', input_mpileup), I_ERBB2]/1.625

  
  ###Find the best anchor number for each segment with minimal sd
  cat('Finding the best anchor number for each segment to  minimize variance withing training samples')
  cat('\n')
  cov_matrix=cor(mean_coorelations)
  max_cov=apply(cov_matrix, 2, function(x) head(order(x, decreasing = TRUE), 100))
  
  ###calculate normal sd
  gene_anchor_sd=list()
  min_ac=80
  for (n_anchor in (min_ac+1):100){
    sds=apply(max_cov, 2, function(x)
        rowMeans(mean_coorelations[,x[2:n_anchor]])/mean_coorelations[,x[1]])
    sds=apply(sds, 2, function(x) sd(x)/mean(x))
    gene_anchor_sd[[n_anchor]]=sds
  }
  gene_anchor_sd=do.call(rbind, gene_anchor_sd)
  
  anchors=apply(gene_anchor_sd, 2, function(x) which.min(x))
  anchors=sapply(anchors, function(x) ifelse(length(x)==0, NA, x))
  
  #### get mean and sd for each distribution
  model_parameters=list()
  model_anchors=list()
  anchors_ratio=list()
  for (i in 1:length(anchors)){
    if (is.na(anchors[i])) {
      model_anchors[[i]]=NA
      model_parameters[[i]]=c(NA, NA) 
      } else{
      model_anchors[[i]]=max_cov[2:(anchors[i]+min_ac),i]
      anchor_mean=rowMeans(mean_coorelations[,model_anchors[[i]]]) 
      loci_depth=mean_coorelations[,i]
      anchors_ratio[[i]]=anchor_mean/loci_depth
      #model_parameters[[i]]=mgedist(anchors_ratio[[i]], 'cauchy', lower = 0.001, upper = 10)$estimate
      model_parameters[[i]]=fitdistr(anchors_ratio[[i]], 'normal')$estimate
      }
  }
#   RD007=36:42
#   M002=grep('M002', input_mpileup)
#   M001=grep('M001', input_mpileup)
#   HD734=grep('HD734', input_mpileup)
#   
#   input_mpileup
#   for (I in 1:20){
# plot(density(anchors_ratio[[I_met[I]]]), ylim=c(0, 20))
# lines(density(anchors_ratio[[I_met[I]]][RD007]), col='red')
# lines(density(anchors_ratio[[I_met[I]]][M002]), col='green')
# lines(density(anchors_ratio[[I_met[I]]][HD734]), col='gold')
# lines(density(anchors_ratio[[I_met[I]]][M001]), col='blue')
#   }
  

  saveRDS(model_parameters, file.path(output_dir, 'model_parameters.rda'))
  
  saveRDS(model_anchors,file.path(output_dir, 'model_anchors.rda'))
cat(paste0('Model files successfully stored'))
cat('\n')
  
