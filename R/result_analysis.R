setwd('E:/Software/30k_cnv_Casar_v2.0/mpileups/')
cnv2.18=list()
for (i in list.files('met2.1875ERBB2.625', pattern = '.cnv')){
  aa=read.table(file.path('met2.1875ERBB2.625', i), header=TRUE)
  cnv2.18[[i]]=aa[4,]
}

cnvnormal=list()
for (i in list.files('normal', pattern = '.cnv')){
  aa=read.table(file.path('normal', i), header=TRUE)
  cnvnormal[[i]]=as.vector(aa[4,])
}
cnvnormal=do.call(rbind, cnvnormal)

###3.25
cnv2.35=list()
for (i in list.files('met2.375ERBB3.25/', pattern = '.cnv')){
  aa=read.table(file.path('met2.375ERBB3.25/', i), header=TRUE)
  cnv2.35[[i]]=as.vector(aa[4,])
}
cnv2.35=do.call(rbind, cnv2.35)

##2.75
cnv2.75=list()
for (i in list.files('met2.75ERBB4.75/', pattern = '.cnv')){
  aa=read.table(file.path('met2.75ERBB4.75/', i), header=TRUE)
  cnv2.75[[i]]=as.vector(aa[4,])
}
cnv2.75=do.call(rbind, cnv2.75)
cnv_all=rbind(cnvnormal, cnv2.35, cnv2.75)

cnv_all$CNV=c(rep('Loss', 10), rep('Gain_1', 6), rep('Gain_2', 6))
cnv_all$ERBB2=as.numeric(as.character(cnv_all$ERBB2))
cnv_all$MET=as.numeric(as.character(cnv_all$MET))
cnv_all$EGFR=as.numeric(as.character(cnv_all$EGFR))
cnv_all$EGFR[is.na(cnv_all$EGFR)]=1
cnv_all$MET[is.na(cnv_all$MET)]=1
cnv_all$ERBB2[is.na(cnv_all$ERBB2)]=1

df=data.frame(value=cbind(c(cnv_all$EGFR, cnv_all$ERBB2, cnv_all$MET)), 
           CNV=rep(cnv_all$CNV,3), Gene=rep(c('EGFR','ERBB2','MET'), each=22))

library(ggplot2)
library(reshape2)
head(df)
ggplot(df, aes(CNV, value))+geom_boxplot(aes(color=Gene))


###clinical data
setwd('E:/Software/30k_cnv_Casar_v2.0//')
cnv=list()
for (i in list.files('clinical_data/', pattern = '.cnv')){
  aa=read.table(file.path('clinical_data/', i), header=TRUE)
  cnv[[i]]=aa[4,]
}
cnv=do.call(rbind, cnv)
dim(cnv)

cnv$ERBB2=as.numeric(as.character(cnv$ERBB2))
cnv$MET=as.numeric(as.character(cnv$MET))
cnv$EGFR=as.numeric(as.character(cnv$EGFR))
cnv$EGFR[is.na(cnv$EGFR)]=1
cnv$MET[is.na(cnv$MET)]=1
cnv$ERBB2[is.na(cnv$ERBB2)]=1

df=data.frame(value=cbind(c(cnv$EGFR, cnv$ERBB2, cnv$MET)), 
              Gene=rep(c('EGFR','ERBB2','MET'), each=35))

library(ggplot2)
library(reshape2)
head(df)
ggplot(df, aes(Gene, value, color=Gene))+geom_jitter(aes(size=abs(value-1)), width=0.5)
