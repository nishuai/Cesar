setwd('E:/wrk/cnv_ctdna/')
K=128  
met_top_I=c(2,16,23,29,46,51,54,55,61,77,126,132,157,158,175,176,190)

anchor_normalized_met_neg=list()
anchor_normalized_met_pos=list()
for (i in 1:length(list.files('inputdata/cnv_samples/inhouse/', 
                                pattern = '*freq'))){
    args=c('inputdata/Master_depth.rda', 
           file.path('inputdata/cnv_samples/inhouse//',  
                     list.files('inputdata/cnv_samples/inhouse/')[i]),
           'temp/')
    source('R/get_loci_in_bed.R')
    ref_bed=read.table('inputdata/segmented_bed_30k.bed')
    ####read the input mpileup file
    input_depth=args[2]
    input_depth=read.table(input_depth, header=TRUE)
    
    total_depth=get_loci_in_bed(input_depth, ref_bed)
    total_depth_mean=mean(input_depth$DEPTH, na.rm=TRUE)
    
    anchor_depth=get_loci_in_bed(input_depth, ref_bed[sort(met_top_I),])
    anchor_depth_mean=mean(anchor_depth$DEPTH, na.rm=TRUE)
    
    
    input_depth=get_loci_in_bed(input_depth, ref_bed[K,])
    total_nor_DEPTH=input_depth$DEPTH/total_depth_mean
    anchor_nor_DEPTH=input_depth$DEPTH/anchor_depth_mean
    anchor_normalized_met_neg[[i]]=total_nor_DEPTH
}
  
  
  ###adding amplified samples
for (i in 1:length(list.files('inputdata/cnv_samples/std///', 
                                pattern = '*freq'))){
    
    args=c('inputdata/Master_depth.rda', 
           file.path('inputdata/cnv_samples/std/',  
                     list.files('inputdata/cnv_samples/std/')[i]),
           'temp/')
    source('R/get_loci_in_bed.R')
    ref_bed=read.table('inputdata/segmented_bed_30k.bed')
    ####read the input mpileup file
    input_depth=args[2]
    input_depth=read.table(input_depth, header=TRUE)
    
    total_depth=get_loci_in_bed(input_depth, ref_bed)
    total_depth_mean=mean(input_depth$DEPTH, na.rm=TRUE)
    
    anchor_depth=get_loci_in_bed(input_depth, ref_bed[c(sort(met_top_I)),])
    anchor_depth_mean=mean(anchor_depth$DEPTH, na.rm=TRUE)
    
    
    input_depth=get_loci_in_bed(input_depth, ref_bed[K,])
    
    total_nor_DEPTH=input_depth$DEPTH/total_depth_mean
    anchor_nor_DEPTH=input_depth$DEPTH/anchor_depth_mean
    anchor_normalized_met_pos[[i]]=total_nor_DEPTH
}


library(MASS)
neg_mean=sapply(anchor_normalized_met_neg, function(x) mean(x))
pos_mean=sapply(anchor_normalized_met_pos, function(x) mean(x))

fit_neg <- fitdistr(neg_mean, 'normal')
fit_pos <- fitdistr(pos_mean, 'normal')


x=seq(0.55, 0.82, 0.001)
plot(dnorm(x,mean=fit_neg$estimate['mean'],sd=fit_neg$estimate['sd']), type='l')
lines(dnorm(x,mean=fit_pos$estimate['mean'],sd=fit_pos$estimate['sd']), type='l')
