git.dir='~/Dropbox/mutation_rate/20180705miseq/'
out.dir=paste0(git.dir, 'out_downsample/') 
sample.alignments=readRDS(paste0(out.dir, 'sample.alignments.RDS'))
##############################
#record mutation type
variant=c("31451_A_G","31663_A_G","31867_T_C","33002_A_C","33344_C_T","33551_A_G","33613_-_A","33782_C_A","33851_C_T","34040_A_G","34050_G_A","34052_T_A","34324_-_A","34493_C_T","34622_T_C")   
mutation=list()
for(libr in names(sample.alignments) ){
  print(libr)
  for (samp in names(sample.alignments[[libr]])){
    print(samp)
    mutation[[libr]][[samp]]=mutation.spect.table(sample.alignments,libr,samp,variant)
  }
}
save(mutation,file=paste0(out.dir,'mutation.Rdata'))

########################
#expected reads number 
##colnum is the number of segregants of each type
colnum=list(281,230,252,277,281,230,252,277) 
names(colnum)=names(sample.alignments)

exp.read=function(sample.alignments,libr,samp){
  onesample.read=sample.alignments[[libr]][[samp]]
  total.read=sum(onesample.read$rr.lengths)
  exp.read=total.read/colnum[[libr]]
  return(c(total.read,exp.read))
}

exp.read2=function(sample.alignments,libr,samp){
  onesample.read=sample.alignments[[libr]][[samp]]
  total.read=sum(onesample.read$rr.lengths)
  return(c(total.read))
}

read=function(sample.alignments,libr,samp){
  onesample.read=sample.alignments[[libr]][[samp]]
  total.read=sum(onesample.read$rr.lengths)
  #exp.read=total.read/colnum[[libr]]
  return(total.read)
}

#########################
#just see how many reads match from start how many not
#count total read count and the cutoff read count (but the exp is underestimated; because some spots may not have colony)
match.data=list()
readcount=list()
readout=list()
readcount2=list()

##colnum is the number of segregants of each type
colnum=list(281,230,252,277,281,230,252,277)
names(colnum)=names(sample.alignments)

for(libr in names(sample.alignments) ){
  print(libr)
  for (samp in names(sample.alignments[[libr]])){
    print(samp)
    #match.data[[libr]][[samp]]=match.con(sample.alignments,libr,samp)
    readcount[[libr]][[samp]]=exp.read(sample.alignments,libr,samp)
    readcount2[[libr]][[samp]]=exp.read2(sample.alignments,libr,samp)
    #readout[[libr]][[samp]]=read(sample.alignments,libr,samp)
  }
}
save(match.data,file=paste0(out.dir,'match.data.Rdata'))
save(readcount,file=paste0(out.dir,'readcount.Rdata'))
save(readout,file=paste0(out.dir,'readout.Rdata'))    ###counts of whole library
load(file=paste0(out.dir,'readcount.Rdata'))

########count of the library
lib.count=c()
for (libr in names(sample.alignments)){
lib.count[[libr]]=sum(unlist(readcount2[[libr]]))
}
#names(lib.count)=c("RAD5:R;MKT1:B","RAD5:R;MKT1:Br","RAD5:R;MKT1:R","RAD5:R;MKT1:Rr","RAD5:B;MKT1:B","RAD5:B;MKT1:Br","RAD5:B;MKT1:R","RAD5:B;MKT1:Rr")
names(lib.count)=c("RAD5:R;MKT1:B","RAD5:R;MKT1:R","RAD5:B;MKT1:B","RAD5:B;MKT1:R")
barplot(lib.count,col="light yellow",xlab="sample",ylab="read.counts")

###print the readcounts
#names=c("RAD5:R;MKT1:R","RAD5:R;MKT1:B","RAD5:B;MKT1:R","RAD5:B;MKT1:B")
names=c("Group1","Group2","Group3","Group4")
lib.count=data.frame(count = as.vector(lib.count), ID = names)
df = melt(lib.count)
p=ggplot(data=lib.count, aes(ID,count)) + 
  geom_bar(stat="identity")
print(p)

######read count sum
total.lib.read=unique(unlist(readcount)[seq(1,length(unlist(type.sum)),2)])
c(sum(total.lib.read[1:12]),sum(total.lib.read[12:23]),sum(total.lib.read[24:35]),sum(total.lib.read[36:47]),sum(total.lib.read[48:59]),sum(total.lib.read[60:70]),sum(total.lib.read[71:82]),sum(total.lib.read[83:94]))
###the distribution of the reads
hist(sample.alignments$LG1_S1$`LG1_S1:CAN1_11_12:chrV:31636:31915`$rr.lengths,breaks = 10000,xlim=c(0,100),col="light yellow",main="sample1_region11",xlab = "read counts for a certain mutation type")
abline(v=20,col="red")

########################
#remove reads shorter than the expected one and re_cal the mutation
#filter BY RM variants ref_alt
variant=c("31451_A_G","31663_A_G","31867_T_C","33002_A_C","33344_C_T","33551_A_G","33613_-_A","33782_C_A","33851_C_T","34040_A_G","34050_G_A","34052_T_A","34324_-_A","34493_C_T","34622_T_C")   
mutation.fil=list()

for(libr in names(sample.alignments) ){
  print(libr)
  for (samp in names(sample.alignments[[libr]])){
    print(samp)
    #filter out with very little reads
    if (readcount[[libr]][[samp]][1] < 200) next
    mutation.fil[[libr]][[samp]]=mutation.spect.table.read(sample.alignments,libr,samp,variant)
  }
}
save(mutation.fil,file=paste0(out.dir,'mutation.fil.Rdata'))
load(file=paste0(out.dir,'mutation.fil.Rdata'))

#########################
#generate the mutation distribution status; make sure mutation not just at two ends
mut.dis=list()
for(libr in names(sample.alignments)){
  print(libr)
  for (samp in names(mutation[[libr]])){
    print(samp)
    mut.dis[[libr]][[samp]]=mut.pos(libr,samp,mutation)
  }
}
save(mut.dis,file=paste0(out.dir,'mut.dis.Rdata'))
load(file=paste0(out.dir,'mut.dis.Rdata'))

###
mut.fil.dis=list()
for(libr in names(sample.alignments)){
  print(libr)
  for (samp in names(mutation.fil[[libr]])){
    print(samp)
    mut.fil.dis[[libr]][[samp]]=mut.pos(libr,samp,mutation.fil)
  }
}
save(mut.fil.dis,file=paste0(out.dir,'mut.fil.dis.Rdata'))
load(file=paste0(out.dir,'mut.fil.dis.Rdata'))

#########################
#replication reproducibility
#the ranking kendall correlation of the intercept of two replicates 
#Kendall's coefficient of concordance
cor=matrix(NA,nrow=4,ncol=10)
for (i in 1:4) {
  for (r in 1:10){
    repone=mutation.fil[[(i*2-1)]][[r]]
    reptwo=mutation.fil[[(i*2)]][[r]]

    inter=intersect(repone[1,],reptwo[1,])
    reponei=repone[,which(repone[1,] %in% inter)]
    reptwoi=reptwo[,which(reptwo[1,] %in% inter)]

    cor[i,r]=cor(as.numeric(reponei[2,]),as.numeric(reptwoi[2,]), method = "kendall")
  }
}
rownames(cor)=c("LG1","LG2","LG3","LG4")
#colnames(cor)=names(mutation.fil$LG1_S1)
colnames(cor)=c("1","10","11","12","3","4","5","6","7","8")
save(cor,file=paste0(out.dir,'cor.Rdata'))
load(file=paste0(out.dir,'cor.Rdata'))

#########################
#summary of mutation types
type.sum=list()
for (i in names(mutation.fil)) {
  for (j in names(mutation.fil[[i]])){    
    types=c()
    for (r in 1:length(mutation.fil[[i]][[j]][1,])){
    twoletter=strsplit(mutation.fil[[i]][[j]][1,],'_')[[r]][2:3]
    type=paste0(as.character(twoletter),collapse='_')
    types=c(types,type)
    }
    types=unique(types)
  
    type.count=c()
    for (k in 1:length(types)){
    type.index=grep(types[k],mutation.fil[[i]][[j]][1,])
    type.count[k]=sum(as.numeric(mutation.fil[[i]][[j]][2,type.index]))
    }
    type.sum[[i]][[j]]=rbind(types,type.count)
  }
}
save(type.sum,file=paste0(out.dir,'type.sum.Rdata'))
load(file=paste0(out.dir,'type.sum.Rdata'))

#########################
#combine mutation types of a group
all_types=unique(unlist(type.sum)[seq(1,length(unlist(type.sum)),2)])
all_sam = names(type.sum)
all_types_rec=matrix(0, nrow=length(all_sam), ncol=length(all_types))
colnames(all_types_rec)=all_types
rownames(all_types_rec)=all_sam

unlist_type_sum = unlist(type.sum)
for(i in seq(1,length(unlist_type_sum),2))
{
  sam=strsplit(names(unlist_type_sum[i]),"[.]")[[1]][1]
  typ=unlist_type_sum[i]
  count=as.numeric(unlist_type_sum[i+1])
  a=which(all_sam==sam)
  b=which(all_types==typ)
  all_types_rec[a,b] = all_types_rec[a,b]+count
}
save(all_types_rec,file=paste0(out.dir,'all_types_rec.Rdata'))
load(file=paste0(out.dir,'all_types_rec.Rdata'))

#all_types_rec_por=round((all_types_rec/apply(all_types_rec,1,sum)),digit=3)*100
all_types_rec_por=round((all_types_rec/lib.count),digit=3)*1000

png(file=paste0(out.dir,'musum_por.png'), width=1000, height=500)
heatmap.2(all_types_rec, Rowv=F, Colv=F, dendrogram='none', trace='none', scale='row')
#heatmap.2(all_types_rec_por, Rowv=F, Colv=F, dendrogram='none', trace='none', scale='row')
dev.off()

#########################
#the unique overlapping proportion between replicates
#
overlapp=list()
for (i in 1:4) {
  overlapp[[i]]=list()
  for (r in 1:10){

    repone=mutation.fil[[(i*2-1)]][[r]]
    reptwo=mutation.fil[[(i*2)]][[r]]
    
    uni=union(repone[1,],reptwo[1,])
    uni_count=c()
    uni_count[which(uni %in% repone[1,])]=1
    uni_count[which(uni %in% reptwo[1,])]=uni_count[which(uni %in% reptwo[1,])]+1
   
    overlapp[[i]][[r]]=rbind(uni,uni_count)
  }
}

names(overlapp)=c("LG1","LG2","LG3","LG4")
for (libr in names(overlapp)){
names(overlapp[[libr]])=c("1","10","11","12","3","4","5","6","7","8")
}

save(overlapp,file=paste0(out.dir,'overlapp.Rdata'))
load(file=paste0(out.dir,'overlapp.Rdata'))

#########################
#the unique overlapping proportion between replicates
chi=matrix(NA,nrow=4,ncol=10)
fisher=matrix(NA,nrow=4,ncol=10)
df.out=matrix(NA,nrow=8,ncol=20)

for (i in 1:4) {
  for (r in 1:10){
    repone=mutation.fil[[(i*2-1)]][[r]]
    reptwo=mutation.fil[[(i*2)]][[r]]
    
    inter=intersect(repone[1,],reptwo[1,])
    reponei=repone[,which(repone[1,] %in% inter)]
    reptwoi=reptwo[,which(reptwo[1,] %in% inter)]
    
    one.non=length(repone[1,])-length(inter)
    two.non=length(reptwo[1,])-length(inter)
    
    one=c(length(inter),one.non)
    two=c(length(inter),two.non)
    
    df.out[(i*2-1),(r*2-1):(r*2)]=one
    df.out[(i*2),(r*2-1):(r*2)]=two
    
    ###not sure the chi square here
    df=as.data.frame(rbind(one,two))
    names(df)=c('overlap','nonoverlap')
    #chi[i,r]=chisq.test(df)$p.value
    fisher[i,r]=fisher.test(df)$p.value
  }
}

rownames(chi)=c("LG1","LG2","LG3","LG4")
colnames(chi)=c("1","10","11","12","3","4","5","6","7","8")

rownames(df.out)=c("RAD5:R;MKT1:B","RAD5:R;MKT1:Br","RAD5:R;MKT1:R","RAD5:R;MKT1:Rr","RAD5:B;MKT1:B","RAD5:B;MKT1:Br","RAD5:B;MKT1:R","RAD5:B;MKT1:Rr")
#colnames(df.out)=rep(c("overlap","non"),10)
colnames(df.out)=c("ovlap1","non1","ovlap10","non10","ovlap11","non11","ovlap12","non12","ovlap3","non3","ovlap4","non4","ovlap5","non5","ovlap6","non6","ovlap7","non7","ovlap8","non8")

save(df.out,file=paste0(out.dir,'df.out.Rdata'))
load(file=paste0(out.dir,'df.out.Rdata'))

df.out.sum=transform(df.out,sum=rowSums(df.out))
#########################
#find where the 'A_-" are distributed
type.sum=list()
loc=c()
loc.tem=c()

for (i in c("LG1_S1","LG1rep_S5")) {
  for (j in names(mutation.fil[[i]])){    
    for (r in 1:length(mutation.fil[[i]][[j]][1,])){
      twoletter=strsplit(mutation.fil[[i]][[j]][1,],'_')[[r]][2:3]
      type=paste0(as.character(twoletter),collapse='_')
      if (type=='A_-'){
        loc.tem=strsplit(mutation.fil[[i]][[j]][1,],'_')[[r]][1]
        loc=c(loc,loc.tem)
        } else {next}
      #types=c(types,type)
    }
    #types=unique(types)
    
    #type.count=c()
    #for (k in 1:length(types)){
    #  type.index=grep(types[k],mutation.fil[[i]][[j]][1,])
    #  type.count[k]=sum(as.numeric(mutation.fil[[i]][[j]][2,type.index]))
    }
    #type.sum[[i]][[j]]=rbind(types,type.count)
}


save(type.sum,file=paste0(out.dir,'type.sum.Rdata'))
load(file=paste0(out.dir,'type.sum.Rdata'))









