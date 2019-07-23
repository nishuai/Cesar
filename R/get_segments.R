setwd('E:/wrk/cnv_ctdna/')
for (K in 1:190){
  
  pdf(paste0('region_',K,'.pdf'))
for (i in 1:length(list.files('inputdata/normal_samples/normal1', 
                              pattern = '*freq'))){
  args=c('inputdata/Master_depth.rda', 
         file.path('inputdata/normal_samples/normal1//',  
                   list.files('inputdata/normal_samples/normal1//')[i]),
         'temp/')
  source('R/get_loci_in_bed.R')
  ref_bed=read.table('inputdata/segmented_bed_30k.bed')
  ####read the input mpileup file
  input_depth=args[2]
  input_depth=read.table(input_depth, header=TRUE)
  
  total_depth=get_loci_in_bed(input_depth, ref_bed)
  total_depth_mean=mean(input_depth$DEPTH, na.rm=TRUE)
  
  anchor_depth=get_loci_in_bed(input_depth, ref_bed[c(3, 21,54, 59, 62, 129),])
  anchor_depth_mean=mean(anchor_depth$DEPTH, na.rm=TRUE)
  
  
  input_depth=get_loci_in_bed(input_depth, ref_bed[K,])
  # ####find segmentation for regions longer than 500
  # library(PSCBS)
  # long_bed=ref_bed[ref_bed$V3-ref_bed$V2>600,]
  # ref_bed=ref_bed[ref_bed$V3-ref_bed$V2<=600,]
  # 
  # for (i in 1:nrow(long_bed)){
  #   depth_vector=get_loci_in_bed(input_depth, long_bed[i,])
  #   changepoints=find_changepoints(depth_vector$DEPTH)
  #   n=length(changepoints)-1
  #   bed_coor=long_bed[i,2]+changepoints-1
  #   seg_bed=data.frame(V1=rep(long_bed[i,1], n), 
  #                      V2=bed_coor[-length(bed_coor)], 
  #                      V3=bed_coor[-1]+1)
  #   ref_bed=rbind(ref_bed, seg_bed)
  # }
  # ref_bed=ref_bed[order(as.numeric(substr(as.character(ref_bed$V1), 4,7)), ref_bed$V2),]
  # write.table(ref_bed, 'segmented_bed_30k.bed', quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # segments=list()
  # segments[[1]]= 1:segs$end[1]
  # for (i in 2:length(segs$end)){
  #   segments[[i]]=ceiling(segs$end[i-1]):floor(segs$end[i])
  # }
  
  total_nor_DEPTH=input_depth$DEPTH/total_depth_mean
  anchor_nor_DEPTH=input_depth$DEPTH/anchor_depth_mean
  
  if(i==1) plot(total_nor_DEPTH, type='l', ylim=c(0, 3), 
                main=paste('total region normalized: bed region ;',K,  
                           'total length:', 
                           length(input_depth$DEPTH)))
  else lines(total_nor_DEPTH, type='l')
}


###adding amplified samples
for (i in 1:length(list.files('inputdata/cnv_samples//', 
                              pattern = '*freq'))){
  
  args=c('inputdata/Master_depth.rda', 
         file.path('inputdata/cnv_samples//',  
                   list.files('inputdata/cnv_samples/')[i]),
         'temp/')
  source('R/get_loci_in_bed.R')
  ref_bed=read.table('inputdata/segmented_bed_30k.bed')
  ####read the input mpileup file
  input_depth=args[2]
  input_depth=read.table(input_depth, header=TRUE)
  
  total_depth=get_loci_in_bed(input_depth, ref_bed)
  total_depth_mean=mean(input_depth$DEPTH, na.rm=TRUE)
  
  anchor_depth=get_loci_in_bed(input_depth, ref_bed[c(3, 21,54, 59, 62, 129),])
  anchor_depth_mean=mean(anchor_depth$DEPTH, na.rm=TRUE)
  
  
  input_depth=get_loci_in_bed(input_depth, ref_bed[K,])

  total_nor_DEPTH=input_depth$DEPTH/total_depth_mean
  anchor_nor_DEPTH=input_depth$DEPTH/anchor_depth_mean
  
  lines(total_nor_DEPTH, type='l', col='red')
}
  
  dev.off()
}
