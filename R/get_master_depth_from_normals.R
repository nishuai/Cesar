get_master_depth=function(location){
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
  master_panel
}