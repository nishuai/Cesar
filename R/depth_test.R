setwd('E:/wrk/cnv_ctdna/')
source('R/get_loci_in_bed.R')
mybed=read.table('inputdata/segmented_bed_30k.bed')
mean_coorelations=list()

file_list=list.files(path = 'inputdata/cnv_samples/inhouse/',
                    pattern='freq')
for (i in 1:length(file_list)){
  file=read.table(file.path('inputdata/cnv_samples/inhouse/', file_list[i]), 
                  header=TRUE)
  # total_depth=get_loci_in_bed(file,mybed)
  # total_mean_depth=median(total_depth$DEPTH, na.rm = TRUE)
  each_mean_depth=c()
  for (k in 1:nrow(mybed)){
    locis=get_loci_in_bed(file,mybed[k,])
    each_mean_depth=c(each_mean_depth, mean(locis$DEPTH, na.rm=TRUE))
  }
  mean_coorelations[[i]]=each_mean_depth
}


####selecting anchor sites for met gene
mean_coorelations=do.call(rbind, mean_coorelations)
write.table(mean_coorelations, 'mean_coorelations_in_9_no_cnv_stds.txt',
            quote = FALSE, row.names = FALSE, col.names = FALSE)

met_top_I=tail(order(cor(mean_coorelations)[128,]), 20)
met_top_I=met_top_I[!is.na(cor(mean_coorelations)[128,][met_top_I])]
met_top_I=met_top_I[!met_top_I=='128']

plot(rowMeans(mean_coorelations, na.rm = TRUE), mean_coorelations[,128])
plot(rowMeans(mean_coorelations[,met_top_I], na.rm = TRUE),  mean_coorelations[,128])

