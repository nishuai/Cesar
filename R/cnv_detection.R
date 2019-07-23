setwd('E:/wrk/algos/TNER-test/mpileup//')

####generate a average depth profile
master_panel=list()
setwd('E:/wrk/cnv_ctdna/inputdata/normal_samples/normal1//')
for (i in list.files()){
  aa=read.table(i, header=TRUE)
 aa=aa[,1:3]
  master_panel[[i]]=aa
}

master_panel=do.call(cbind, master_panel)
master_panel=cbind(master_panel[,1:2],
rowMeans(master_panel[,grep('DEPTH', names(master_panel))]))
names(master_panel)=c('CHR','POSITION','DEPTH')



###generate for each sample
for (i in list.files()){
  aa=read.table(i, header=TRUE)
  vafs=(aa[,c(7,9,11,13)] + aa[,c(7,9,11,13)+1])/aa$DEPTH
  ###remove invalid rows in aa
  aa=aa[!rowSums(vafs)>1,]
  vafs=vafs[!rowSums(vafs)>1,]
  ###
  # plot(density(vafs[vafs>0.2 & vafs<0.8], na.rm = TRUE, bw = 0.001))
  SNPs=vafs[rowSums(vafs>0.2 & vafs<1)>=1,]
  plot(aa$POSITION, type='l', axes =FALSE, xlab='', ylab='', col='red')
  par(new=TRUE)
  plot(rownames(SNPs), apply(SNPs, 1, function(x) max(max(x), 1-max(x))), pch=19, col='blue')
  abline(h=0.5, col='blue')
  ###coverage plot
  window_size=20
  seq_depth1= sapply(seq(1,length(aa$DEPTH), window_size), 
                    function(x) mean(aa$DEPTH[x:x+window_size]))
  # plot(aa$POSITION, type='l', axes =FALSE, xlab='', ylab='', col='red')
  seq_depth2= sapply(seq(1,length(master_panel$DEPTH), window_size), 
                    function(x) mean(master_panel$DEPTH[x:x+window_size]))
  par(new=TRUE)
  cov_ratio=log2(seq_depth1/seq_depth2)
  cov_ratio=cov_ratio[!is.na(cov_ratio)]
  cov_ratio=cov_ratio-mean(cov_ratio)+1
  plot(cov_ratio, type='l', col='green')
  abline(h=1, col='green')
  
}






