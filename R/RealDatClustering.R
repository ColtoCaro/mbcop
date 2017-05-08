#Code to apply CSENS and CSPEC to a pair of real data sets

setwd('C:/Users/Jonathon/Dropbox/Qaqish/Clustering/QsClustering/')
cancerdat<-read.table('Multi_A.txt',header=T)

library(cluster)
library(miscFuncs)

#Define functions

clus1 = function(x,k)
  #implements kmeans clustering for k groups on sample x
  #returns a vector of integers from [1,..,k]
{
  clust<-kmeans(x,k,nstart=10)
  return(clust$cluster)
}

clus2 = function(x,k)
  #implements HC with average linkange for k groups on sample x
  #returns a vector of integers from [1,..,k]
{
  clust<-cutree(agnes(x),k=k)
  return(clust)
}

clus3 = function(x,k)
  #implements HC with single linkage for k groups on sample x
  #returns a vector of integers from [1,..,k]
{
  clust<-cutree(agnes(x,method="single"),k=k)
  return(clust)
}

clus4 = function(x,k)
  #implements HC with complete linkange for k groups on sample x
  #returns a vector of integers from [1,..,k]
{
  clust<-cutree(agnes(x,method="complete"),k=k)
  return(clust)
}


ch2 = function(n) {n * (n - 1) / 2}

smry2x2 = function(a, g) {
  # a : true cluster id
  # g : assigned cluster id
  # return the 2x2 table of pairs (as a vector)
  n   = table(a, g)
  mdd = ch2(sum(n))  # total #pairs
  m22 = sum(ch2(n))  # #links by both
  m2d = sum(ch2(apply(n, 1, sum)))  #true links
  md2 = sum(ch2(apply(n, 2, sum)))  #assigned links
  m12 = md2 - m22
  m21 = m2d - m22
  m11 = mdd - m2d - md2 + m22
  return(c(m11, m12, m21, m22))
}

get_indices = function(freq)
{
  if(is.null(dim(freq))){
    sens = freq[4] / (freq[3] + freq[4])
    spec = freq[1] / (freq[1] + freq[2])
  }
  if(!(is.null(dim(freq)))){
    sens = freq[,4] / (freq[,3] + freq[,4])
    spec = freq[,1] / (freq[,1] + freq[,2])
  }
  return(list(sens=sens, spec=spec))
}


#Analyze the cancer dat
truelab<-substr(colnames(cancerdat[,-c(1:2)]),1,2)


labelfreq<-table(truelab)/length(truelab)
prev = labelfreq %*% labelfreq  # "known" prevalence

sens = matrix(NA,nrow=4,ncol=7)    # sensitivity
spec = matrix(NA,nrow=4,ncol=7)    # specificity
sums = matrix(NA,nrow=4,ncol=7)
ppv  = matrix(NA,nrow=4,ncol=7)
npv  = matrix(NA,nrow=4,ncol=7)
sumpv = matrix(NA,nrow=4,ncol=7)

for(k in 2:8){
g1 = clus1(t(cancerdat[,-c(1:2)]),k)
g2 = clus2(t(cancerdat[,-c(1:2)]),k)
g3 = clus3(t(cancerdat[,-c(1:2)]),k)
g4 = clus4(t(cancerdat[,-c(1:2)]),k)

freq1 = smry2x2(truelab, g1)
freq2 = smry2x2(truelab, g2)
freq3 = smry2x2(truelab, g3)
freq4 = smry2x2(truelab, g4)  
  
indi1 = get_indices(freq1)
indi2 = get_indices(freq2)
indi3 = get_indices(freq3)
indi4 = get_indices(freq4)

sens[1,k-1] = indi1$sens    # sensitivity
spec[1,k-1] = indi1$spec    # specificity
sums[1,k-1] = sens[1,k-1]+spec[1,k-1]
ppv[1,k-1]  = prev*sens[1,k-1] / (prev*sens[1,k-1]+(1-prev)*(1-spec[1,k-1]))
npv[1,k-1]  = (1-prev)*spec[1,k-1] / ((1-prev)*spec[1,k-1] + prev*(1-sens[1,k-1]))
sumpv[1,k-1] = ppv[1,k-1]+npv[1,k-1]
   
  sens[2,k-1] = indi2$sens    # sensitivity
  spec[2,k-1] = indi2$spec    # specificity
  sums[2,k-1] = sens[2,k-1]+spec[2,k-1]
  ppv[2,k-1]  = prev*sens[2,k-1] / (prev*sens[2,k-1]+(1-prev)*(1-spec[2,k-1]))
  npv[2,k-1]  = (1-prev)*spec[2,k-1] / ((1-prev)*spec[2,k-1] + prev*(1-sens[2,k-1]))
  sumpv[2,k-1] = ppv[2,k-1]+npv[2,k-1]
  
  sens[3,k-1] = indi3$sens    # sensitivity
  spec[3,k-1] = indi3$spec    # specificity
  sums[3,k-1] = sens[3,k-1]+spec[3,k-1]
  ppv[3,k-1]  = prev*sens[3,k-1] / (prev*sens[3,k-1]+(1-prev)*(1-spec[3,k-1]))
  npv[3,k-1]  = (1-prev)*spec[3,k-1] / ((1-prev)*spec[3,k-1] + prev*(1-sens[3,k-1]))
  sumpv[3,k-1] = ppv[3,k-1]+npv[3,k-1]

  sens[4,k-1] = indi4$sens    # sensitivity
  spec[4,k-1] = indi4$spec    # specificity
  sums[4,k-1] = sens[4,k-1]+spec[4,k-1]
  ppv[4,k-1]  = prev*sens[4,k-1] / (prev*sens[4,k-1]+(1-prev)*(1-spec[4,k-1]))
  npv[4,k-1]  = (1-prev)*spec[4,k-1] / ((1-prev)*spec[4,k-1] + prev*(1-sens[4,k-1]))
  sumpv[4,k-1] = ppv[4,k-1]+npv[4,k-1]
  
  
}#end for loop for k  


plot(c(2:8),sums[1,],main="Sensitivity + Specificity for Mixed Tumor Data \n True Groups = 4",
  ylab="CSENS+CSPEC",xlab="Number of Clusters",type='o',col=1,ylim=c(1,2.5))
points(c(2:8),sums[2,],type='o',col=2,lty=2)
points(c(2:8),sums[3,],type='o',col=3,lty=3)
points(c(2:8),sums[4,],type='o',col=4,lty=4)
legend("topright",c("K-Means","HC:Average","HC:Single","HC:Complete"),fill=c(1,2,3,4),lty=c(1,2,3,4))

rnames<-c("K-Means","Average Linkage HC","Single Linkage HC","Complete Linkage HC")
latextable(sens,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(spec,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sums,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(ppv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(npv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sumpv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)


#Cross tabulation to explain plot

HCavg = clus2(t(cancerdat[,-c(1:2)]),6)
crosstab<-table(truelab,HCavg)
latextable(crosstab,colnames=c('','1','2','3','4','5','6'),rownames=c("Breast","Colon","Lung","Prostate"))


HCavg = clus2(t(cancerdat[,-c(1:2)]),7)
crosstab<-table(truelab,HCavg)
latextable(crosstab,colnames=c('','1','2','3','4','5','6','7'),rownames=c("Breast","Colon","Lung","Prostate"))


HCavg = clus2(t(cancerdat[,-c(1:2)]),8)
crosstab<-table(truelab,HCavg)
latextable(crosstab,colnames=c('','1','2','3','4','5','6','7','8'),rownames=c("Breast","Colon","Lung","Prostate"))



############Now do the lung cancer data############
cancerdat<-read.csv('DatasetA.csv')

truelab<-substr(colnames(cancerdat[,-c(1:2)]),1,2)


labelfreq<-table(truelab)/length(truelab)
prev = labelfreq %*% labelfreq  # "known" prevalence

sens = matrix(NA,nrow=4,ncol=7)    # sensitivity
spec = matrix(NA,nrow=4,ncol=7)    # specificity
sums = matrix(NA,nrow=4,ncol=7)
ppv  = matrix(NA,nrow=4,ncol=7)
npv  = matrix(NA,nrow=4,ncol=7)
sumpv = matrix(NA,nrow=4,ncol=7)

for(k in 2:8){
  g1 = clus1(t(cancerdat[,-c(1:2)]),k)
  g2 = clus2(t(cancerdat[,-c(1:2)]),k)
  g3 = clus3(t(cancerdat[,-c(1:2)]),k)
  g4 = clus4(t(cancerdat[,-c(1:2)]),k)
  
  freq1 = smry2x2(truelab, g1)
  freq2 = smry2x2(truelab, g2)
  freq3 = smry2x2(truelab, g3)
  freq4 = smry2x2(truelab, g4)  
  
  indi1 = get_indices(freq1)
  indi2 = get_indices(freq2)
  indi3 = get_indices(freq3)
  indi4 = get_indices(freq4)
  
  sens[1,k-1] = indi1$sens    # sensitivity
  spec[1,k-1] = indi1$spec    # specificity
  sums[1,k-1] = sens[1,k-1]+spec[1,k-1]
  ppv[1,k-1]  = prev*sens[1,k-1] / (prev*sens[1,k-1]+(1-prev)*(1-spec[1,k-1]))
  npv[1,k-1]  = (1-prev)*spec[1,k-1] / ((1-prev)*spec[1,k-1] + prev*(1-sens[1,k-1]))
  sumpv[1,k-1] = ppv[1,k-1]+npv[1,k-1]
  
  sens[2,k-1] = indi2$sens    # sensitivity
  spec[2,k-1] = indi2$spec    # specificity
  sums[2,k-1] = sens[2,k-1]+spec[2,k-1]
  ppv[2,k-1]  = prev*sens[2,k-1] / (prev*sens[2,k-1]+(1-prev)*(1-spec[2,k-1]))
  npv[2,k-1]  = (1-prev)*spec[2,k-1] / ((1-prev)*spec[2,k-1] + prev*(1-sens[2,k-1]))
  sumpv[2,k-1] = ppv[2,k-1]+npv[2,k-1]
  
  sens[3,k-1] = indi3$sens    # sensitivity
  spec[3,k-1] = indi3$spec    # specificity
  sums[3,k-1] = sens[3,k-1]+spec[3,k-1]
  ppv[3,k-1]  = prev*sens[3,k-1] / (prev*sens[3,k-1]+(1-prev)*(1-spec[3,k-1]))
  npv[3,k-1]  = (1-prev)*spec[3,k-1] / ((1-prev)*spec[3,k-1] + prev*(1-sens[3,k-1]))
  sumpv[3,k-1] = ppv[3,k-1]+npv[3,k-1]
  
  sens[4,k-1] = indi4$sens    # sensitivity
  spec[4,k-1] = indi4$spec    # specificity
  sums[4,k-1] = sens[4,k-1]+spec[4,k-1]
  ppv[4,k-1]  = prev*sens[4,k-1] / (prev*sens[4,k-1]+(1-prev)*(1-spec[4,k-1]))
  npv[4,k-1]  = (1-prev)*spec[4,k-1] / ((1-prev)*spec[4,k-1] + prev*(1-sens[4,k-1]))
  sumpv[4,k-1] = ppv[4,k-1]+npv[4,k-1]
  
  
}#end for loop for k  


plot(c(2:8),sums[1,],main="Sensitivity + Specificity for Lung Cancer Data \n Histological Groups = 5",
  ylab="CSENS+CSPEC",xlab="Number of Clusters",type='o',col=1,ylim=c(1,2))
points(c(2:8),sums[2,],type='o',col=2,lty=2)
points(c(2:8),sums[3,],type='o',col=3,lty=3)
points(c(2:8),sums[4,],type='o',col=4,lty=4)
legend("topright",c("K-Means","HC:Average","HC:Single","HC:Complete"),fill=c(1,2,3,4),lty=c(1,2,3,4))

rnames<-c("K-Means","Average Linkage HC","Single Linkage HC","Complete Linkage HC")
latextable(sens,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(spec,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sums,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(ppv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(npv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sumpv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)




#########################Repeat above analysis with reduced data#####################
cancerdat<-read.table('Multi_A.txt',header=T)
mads<-apply(as.matrix(cancerdat[,-c(1:2)]),1,mad)
index<-order(mads,decreasing=TRUE)
index<-index[1:500]
cancerdat<-cancerdat[index,]

truelab<-substr(colnames(cancerdat[,-c(1:2)]),1,2)

labelfreq<-table(truelab)/length(truelab)
prev = labelfreq %*% labelfreq  # "known" prevalence

sens = matrix(NA,nrow=4,ncol=7)    # sensitivity
spec = matrix(NA,nrow=4,ncol=7)    # specificity
sums = matrix(NA,nrow=4,ncol=7)
ppv  = matrix(NA,nrow=4,ncol=7)
npv  = matrix(NA,nrow=4,ncol=7)
sumpv = matrix(NA,nrow=4,ncol=7)

for(k in 2:8){
  g1 = clus1(t(cancerdat[,-c(1:2)]),k)
  g2 = clus2(t(cancerdat[,-c(1:2)]),k)
  g3 = clus3(t(cancerdat[,-c(1:2)]),k)
  g4 = clus4(t(cancerdat[,-c(1:2)]),k)
  
  freq1 = smry2x2(truelab, g1)
  freq2 = smry2x2(truelab, g2)
  freq3 = smry2x2(truelab, g3)
  freq4 = smry2x2(truelab, g4)  
  
  indi1 = get_indices(freq1)
  indi2 = get_indices(freq2)
  indi3 = get_indices(freq3)
  indi4 = get_indices(freq4)
  
  sens[1,k-1] = indi1$sens    # sensitivity
  spec[1,k-1] = indi1$spec    # specificity
  sums[1,k-1] = sens[1,k-1]+spec[1,k-1]
  ppv[1,k-1]  = prev*sens[1,k-1] / (prev*sens[1,k-1]+(1-prev)*(1-spec[1,k-1]))
  npv[1,k-1]  = (1-prev)*spec[1,k-1] / ((1-prev)*spec[1,k-1] + prev*(1-sens[1,k-1]))
  sumpv[1,k-1] = ppv[1,k-1]+npv[1,k-1]
  
  sens[2,k-1] = indi2$sens    # sensitivity
  spec[2,k-1] = indi2$spec    # specificity
  sums[2,k-1] = sens[2,k-1]+spec[2,k-1]
  ppv[2,k-1]  = prev*sens[2,k-1] / (prev*sens[2,k-1]+(1-prev)*(1-spec[2,k-1]))
  npv[2,k-1]  = (1-prev)*spec[2,k-1] / ((1-prev)*spec[2,k-1] + prev*(1-sens[2,k-1]))
  sumpv[2,k-1] = ppv[2,k-1]+npv[2,k-1]
  
  sens[3,k-1] = indi3$sens    # sensitivity
  spec[3,k-1] = indi3$spec    # specificity
  sums[3,k-1] = sens[3,k-1]+spec[3,k-1]
  ppv[3,k-1]  = prev*sens[3,k-1] / (prev*sens[3,k-1]+(1-prev)*(1-spec[3,k-1]))
  npv[3,k-1]  = (1-prev)*spec[3,k-1] / ((1-prev)*spec[3,k-1] + prev*(1-sens[3,k-1]))
  sumpv[3,k-1] = ppv[3,k-1]+npv[3,k-1]
  
  sens[4,k-1] = indi4$sens    # sensitivity
  spec[4,k-1] = indi4$spec    # specificity
  sums[4,k-1] = sens[4,k-1]+spec[4,k-1]
  ppv[4,k-1]  = prev*sens[4,k-1] / (prev*sens[4,k-1]+(1-prev)*(1-spec[4,k-1]))
  npv[4,k-1]  = (1-prev)*spec[4,k-1] / ((1-prev)*spec[4,k-1] + prev*(1-sens[4,k-1]))
  sumpv[4,k-1] = ppv[4,k-1]+npv[4,k-1]
  
  
}#end for loop for k  


plot(c(2:8),sums[1,],main="Sensitivity + Specificity for Mixed Tumor Data \n 500 top MADs,True Groups = 4",
  ylab="CSENS+CSPEC",xlab="Number of Clusters",type='o',col=1,ylim=c(1,2))
points(c(2:8),sums[2,],type='o',col=2)
points(c(2:8),sums[3,],type='o',col=3)
points(c(2:8),sums[4,],type='o',col=4)
legend("topright",c("K-Means","HC:Average","HC:Single","HC:Complete"),fill=c(1,2,3,4))

rnames<-c("K-Means","Average Linkage HC","Single Linkage HC","Complete Linkage HC")
latextable(sens,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(spec,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sums,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(ppv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(npv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sumpv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)

#lung data
cancerdat<-read.csv('DatasetA.csv')
mads<-apply(as.matrix(cancerdat[,-c(1:2)]),1,mad)
index<-order(mads,decreasing=TRUE)
index<-index[1:500]
cancerdat<-cancerdat[index,]

truelab<-substr(colnames(cancerdat[,-c(1:2)]),1,2)

labelfreq<-table(truelab)/length(truelab)
prev = labelfreq %*% labelfreq  # "known" prevalence

sens = matrix(NA,nrow=4,ncol=7)    # sensitivity
spec = matrix(NA,nrow=4,ncol=7)    # specificity
sums = matrix(NA,nrow=4,ncol=7)
ppv  = matrix(NA,nrow=4,ncol=7)
npv  = matrix(NA,nrow=4,ncol=7)
sumpv = matrix(NA,nrow=4,ncol=7)

for(k in 2:8){
  g1 = clus1(t(cancerdat[,-c(1:2)]),k)
  g2 = clus2(t(cancerdat[,-c(1:2)]),k)
  g3 = clus3(t(cancerdat[,-c(1:2)]),k)
  g4 = clus4(t(cancerdat[,-c(1:2)]),k)
  
  freq1 = smry2x2(truelab, g1)
  freq2 = smry2x2(truelab, g2)
  freq3 = smry2x2(truelab, g3)
  freq4 = smry2x2(truelab, g4)  
  
  indi1 = get_indices(freq1)
  indi2 = get_indices(freq2)
  indi3 = get_indices(freq3)
  indi4 = get_indices(freq4)
  
  sens[1,k-1] = indi1$sens    # sensitivity
  spec[1,k-1] = indi1$spec    # specificity
  sums[1,k-1] = sens[1,k-1]+spec[1,k-1]
  ppv[1,k-1]  = prev*sens[1,k-1] / (prev*sens[1,k-1]+(1-prev)*(1-spec[1,k-1]))
  npv[1,k-1]  = (1-prev)*spec[1,k-1] / ((1-prev)*spec[1,k-1] + prev*(1-sens[1,k-1]))
  sumpv[1,k-1] = ppv[1,k-1]+npv[1,k-1]
  
  sens[2,k-1] = indi2$sens    # sensitivity
  spec[2,k-1] = indi2$spec    # specificity
  sums[2,k-1] = sens[2,k-1]+spec[2,k-1]
  ppv[2,k-1]  = prev*sens[2,k-1] / (prev*sens[2,k-1]+(1-prev)*(1-spec[2,k-1]))
  npv[2,k-1]  = (1-prev)*spec[2,k-1] / ((1-prev)*spec[2,k-1] + prev*(1-sens[2,k-1]))
  sumpv[2,k-1] = ppv[2,k-1]+npv[2,k-1]
  
  sens[3,k-1] = indi3$sens    # sensitivity
  spec[3,k-1] = indi3$spec    # specificity
  sums[3,k-1] = sens[3,k-1]+spec[3,k-1]
  ppv[3,k-1]  = prev*sens[3,k-1] / (prev*sens[3,k-1]+(1-prev)*(1-spec[3,k-1]))
  npv[3,k-1]  = (1-prev)*spec[3,k-1] / ((1-prev)*spec[3,k-1] + prev*(1-sens[3,k-1]))
  sumpv[3,k-1] = ppv[3,k-1]+npv[3,k-1]
  
  sens[4,k-1] = indi4$sens    # sensitivity
  spec[4,k-1] = indi4$spec    # specificity
  sums[4,k-1] = sens[4,k-1]+spec[4,k-1]
  ppv[4,k-1]  = prev*sens[4,k-1] / (prev*sens[4,k-1]+(1-prev)*(1-spec[4,k-1]))
  npv[4,k-1]  = (1-prev)*spec[4,k-1] / ((1-prev)*spec[4,k-1] + prev*(1-sens[4,k-1]))
  sumpv[4,k-1] = ppv[4,k-1]+npv[4,k-1]
  
  
}#end for loop for k  


plot(c(2:8),sums[1,],main="Sensitivity + Specificity for Lung Cancer Data \n 500 top MADs, True Groups = 5",
  ylab="CSENS+CSPEC",xlab="Number of Clusters",type='o',col=1,ylim=c(1,2))
points(c(2:8),sums[2,],type='o',col=2)
points(c(2:8),sums[3,],type='o',col=3)
points(c(2:8),sums[4,],type='o',col=4)
legend("topright",c("K-Means","HC:Average","HC:Single","HC:Complete"),fill=c(1,2,3,4))

rnames<-c("K-Means","Average Linkage HC","Single Linkage HC","Complete Linkage HC")
latextable(sens,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(spec,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sums,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(ppv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(npv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)
latextable(sumpv,colnames=c('','K=2','K=3','K=4','K=5','K=6','K=7','k=8'),rownames=rnames)


