find_changepoints=function(depth_vector){
  library(PSCBS)
  
  segs=segmentByCBS(depth_vector,
                    min.width=5, verbose=0)
  loci=segs$segRows
  segs=segs$output
  ###remove uncessary segments
  ####injacent segments must differ by at least 1.15 times in mean
  mergerable=abs(segs$mean[2:length(segs$mean)]/
                   segs$mean[2:length(segs$mean)-1]-1)>0.1
  
  ####injacent segments must differ by at least 500 depth in mean
  
  mergerable=mergerable & abs(segs$mean[2:length(segs$mean)]-
                                segs$mean[2:length(segs$mean)-1]-1)>total_depth_mean/20
  
  mergerable=mergerable & segs$nbrOfLoci[-1]>5
  ####segment length must span at 6 nucleotide
  mergerable=c(mergerable, TRUE)
  loci=loci[mergerable,]
  changepoints=c(loci$startRow[1], loci$endRow)
  changepoints=changepoints[c(TRUE, diff(changepoints[-length(changepoints)])>50, TRUE)]
  changepoints
}

