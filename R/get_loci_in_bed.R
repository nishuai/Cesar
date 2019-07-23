get_loci_in_bed=function(loci, bed){
  nu_chr=substr(bed[,1], 4, 7)
  nu_chr=gsub( 'X', 23, nu_chr)
  nu_chr=gsub( 'Y', 23, nu_chr)
  nu_chr=as.numeric(nu_chr)*1e10
  bed_pos=c(nu_chr+bed[,2]-1, nu_chr+bed[,3]+1)
  
  ###make sure bed do not overmap, merge if overlaps
  aa=cbind(nu_chr+bed[,2], nu_chr+bed[,3])
  if (any(aa[-1,1]+4<aa[-nrow(aa), 2])) stop('please make sure the input bed file is ordered, and do not overlap by over 5 nucleotide')
  
  nu_chr=substr(loci[,1], 4, 7)
  nu_chr=gsub( 'X', 23, nu_chr)
  nu_chr=gsub( 'Y', 23, nu_chr)
  nu_chr=as.numeric(nu_chr)*1e10
  
  loci_pos=c(nu_chr+loci[,2])
  
  total_pos=c(loci_pos, bed_pos)
  total_order=order(total_pos)
  loci_start=total_order[which(total_order>length(loci_pos))[c(TRUE,FALSE)]+1]
  loci_end=total_order[which(total_order>length(loci_pos))[c(FALSE,TRUE)]-1]
  
  dummy=c()
  for (i in 1:length(loci_start)){
    dummy=c(dummy,loci_start[i]:loci_end[i])
  } 
  
 return(loci[dummy,])
}
