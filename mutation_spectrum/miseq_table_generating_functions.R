##extract difference characters in two DNA sequences (must of same lengh)
list.string.diff<-function(query,ref,samp){
  seq.query<-unlist(strsplit(query,split=""))
  seq.ref<-unlist(strsplit(ref,split=""))
  diff.d<-rbind(seq.query,seq.ref)
  only.diff<-diff.d[,diff.d[1,]!=diff.d[2,]]
  
  start.pos=as.numeric(strsplit(samp,":")[[1]][4])
  pos<-which(diff.d[1,]!=diff.d[2,])+start.pos
  #only.diff<-rbind(pos,only.diff)
  if ((length(pos)>0)==TRUE){
    diff=paste(pos,only.diff[[1]],only.diff[[2]],sep='_')} else {diff = NA}
  return(diff)
}
#####################
#for not removing reads not match from start
###not sure where did I use the diff2
list.string.diff2<-function(query,ref,samp,startindex){
  seq.query<-unlist(strsplit(query,split=""))
  seq.ref<-unlist(strsplit(ref,split=""))
  diff.d<-rbind(seq.query,seq.ref)
  only.diff<-diff.d[,diff.d[1,]!=diff.d[2,]]
  
  start.pos=as.numeric(strsplit(samp,":")[[1]][4])
  pos<-which(diff.d[1,]!=diff.d[2,])+start.pos+startindex-2
  #only.diff<-rbind(pos,only.diff)
  if ((length(pos)>0)==TRUE){
    diff=paste(pos,only.diff[[1]],only.diff[[2]],sep='_')} else {diff = NA}
  return(diff)
}

##generate a table with all mutation/variants types and count information
mutation.spect.table=function(sample.alignments,libr,samp,variant){
  id_vec=c()
  count_vec=c()
  ref.start=c()
  
  onesample.read=sample.alignments[[libr]][[samp]]
  startsite=regexpr(onesample.read$ref[1],onesample.read$ref)
  #filter out reads not match from the start
  index=which(startsite==1)
  
  for (r in index){
    query=onesample.read$query[r]
    ref=onesample.read$ref[r]
    count=onesample.read$rr.lengths[r]
    diff=list.string.diff(query,ref,samp)
    
    if (is.na(diff[1]==TRUE)) next
    
    for (i in 1:length(diff)){
      if(diff[i] %in% id_vec){
        count_vec[which(id_vec==diff[i])]=count_vec[which(id_vec==diff[i])]+count
      } else {
        id_vec=c(id_vec,diff[i])
        count_vec[which(id_vec==diff[i])]=count
      }
    }
  }
  
  if (length(which(id_vec %in% variant))>0) {
    count_vec = count_vec[-which(id_vec %in% variant)]
    id_vec=id_vec[-which(id_vec %in% variant)]
  }
  
  spect=rbind(id_vec,count_vec)
  return(spect)
}


#################not moving reads not from the start
###this is the one I used
mutation.spect.table.read=function(sample.alignments,libr,samp,variant){
  id_vec=c()
  count_vec=c()
  ref.start=c()
  
  onesample.read=sample.alignments[[libr]][[samp]]
  #startsite=regexpr(onesample.read$ref[1],onesample.read$ref)
  
  #startsite=regexpr(onesample.read$ref,onesample.read$ref[1])
  #filter out reads not match from the start
  #index=which(startsite==1)
  
  #filter readout with less number of reads 
  #few.exp=which(onesample.read$rr.lengths > readcount[[libr]][[samp]][2])
  #few.exp=which(onesample.read$rr.lengths > 20)
  few.exp=which(onesample.read$rr.lengths > 8)
  #index=intersect(few.exp,index)
  
  for (r in few.exp){
  #for (r in 1:length(onesample.read$query)){
    query=onesample.read$query[r]
    ref=onesample.read$ref[r]
    count=onesample.read$rr.lengths[r]
    ref.index=gsub('-','',ref)
    startindex=regexpr(ref.index,onesample.read$ref[1])[1]
    diff=list.string.diff2(query,ref,samp,startindex)
    
    if (is.na(diff[1]==TRUE)) next
    
    for (i in 1:length(diff)){
      if(diff[i] %in% id_vec){
        count_vec[which(id_vec==diff[i])]=count_vec[which(id_vec==diff[i])]+count
      } else {
        id_vec=c(id_vec,diff[i])
        count_vec[which(id_vec==diff[i])]=count
      }
    }
  }
  
  #remove the BY and RM variants
  if (length(which(id_vec %in% variant))>0) {
    count_vec=count_vec[-which(id_vec %in% variant)]
    id_vec=id_vec[-which(id_vec %in% variant)]
  }
  
  spect=rbind(id_vec,count_vec)
  return(spect)
}

######just see how many reads match from start how many not
match.con=function(sample.alignments,libr,samp){
  onesample.read=sample.alignments[[libr]][[samp]]
  startsite=regexpr(onesample.read$ref[1],onesample.read$ref)
  #filter out reads not match from the start
  index=which(startsite==1)
  notindex=which(startsite!=1)
  match.start=sum(onesample.read$rr.lengths[index])
  not.match.start=sum(onesample.read$rr.lengths[notindex])
  return(c(match.start,not.match.start))
}

######mutation filtered by variants and readcount
mutation.spect.table.read=function(sample.alignments,libr,samp,variant){
  id_vec=c()
  count_vec=c()
  ref.start=c()
  
  onesample.read=sample.alignments[[libr]][[samp]]
  startsite=regexpr(onesample.read$ref[1],onesample.read$ref)
  #filter out reads not match from the start
  index=which(startsite==1)
  
  #filter readout with less number of reads 
  #few.exp=which(onesample.read$rr.lengths > readcount[[libr]][[samp]][2])
  few.exp=which(onesample.read$rr.lengths > 20)
  index=intersect(few.exp,index)
  
  for (r in index){
    query=onesample.read$query[r]
    ref=onesample.read$ref[r]
    count=onesample.read$rr.lengths[r]
    diff=list.string.diff(query,ref,samp)
    
    if (is.na(diff[1]==TRUE)) next
    
    for (i in 1:length(diff)){
      if(diff[i] %in% id_vec){
        count_vec[which(id_vec==diff[i])]=count_vec[which(id_vec==diff[i])]+count
      } else {
        id_vec=c(id_vec,diff[i])
        count_vec[which(id_vec==diff[i])]=count
      }
    }
  }
  
  #remove the BY and RM variants
  if (length(which(id_vec %in% variant))>0) {
    count_vec=count_vec[-which(id_vec %in% variant)]
    id_vec=id_vec[-which(id_vec %in% variant)]
  }
  
  spect=rbind(id_vec,count_vec)
  return(spect)
}

########################
##distribution of the mutations 
mut.pos=function(libr,samp,mutation){
  mu.count=mutation[[libr]][[samp]][2,]
  mu.pos=c()
  for (r in 1:length(mu.count)){
    tem.pos=strsplit(mutation[[libr]][[samp]][1,],'_')[[r]][1]
    mu.pos=c(mu.pos,tem.pos)
  }
  res=rbind(as.numeric(mu.pos),as.numeric(mu.count))
  return(res)
}

