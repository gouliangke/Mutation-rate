library(reshape2)
library(ggplot2)

all_types = unlist(mutation.fil)
all_types = all_types[seq(1,length(all_types),2)]
all_types = unique(
  unlist(lapply(all_types, function(x) paste0(strsplit(x,"_")[[1]][2:3],collapse="_")
  )
  ))
all_types = gsub("-", "N", all_types)

count_mat = matrix(nrow=0, ncol=length(all_types)+1)
colnames(count_mat) = c('ID', all_types)
count_mat = as.data.frame(count_mat, stringsAsFactors=F)

#the code is for replicates, here I did not have replicates, so 
#just rep the mutation.fil to make it compatible for rest codes
mutation.fil = rep(mutation.fil, each = 2)
replicate_mat = matrix(names(mutation.fil), ncol=2, byrow = T)
rownames(replicate_mat) = c('S1', 'S2', 'S3', 'S4')
#rownames(replicate_mat) = c('RAD5:R;MKT1:R', 'RAD5:R;MKT1:B', 'RAD5:B;MKT1:R', 'RAD5:B;MKT1:B')

for(i in 1:nrow(replicate_mat))
{
  this_id = rownames(replicate_mat)[i]
  this = get_replicated_type(replicate_mat[i,1], replicate_mat[i,2])
  count_mat = expand_count(this_id, this$S_common, count_mat)
}

save(count_mat,file=paste0(out.dir,'count_mat.Rdata'))

#selected_mutation_type = apply(count_mat[,-1],2,sum)<10

#keep the one base pair indels
#selected_mutation_type = grep("N", colnames(count_mat), invert=T)

df = melt(count_mat[,selected_mutation_type])
p=ggplot(data=df, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") +
    scale_fill_brewer(palette = "Set3") + theme_bw()
print(p)
#change position="dodge" to position="stack" to make stack graph

##make a combine dataframe (just 6 spectrum) for count_mat
#count_mat_tmp=count_mat[,selected_mutation_type]
count_mat_tmp=count_mat
count_mat_new=data.frame(
                         CG_GC = count_mat_tmp$C_G+count_mat_tmp$G_C,
                         CA_GT = count_mat_tmp$C_A+count_mat_tmp$G_T,
                         CT_GA = count_mat_tmp$C_T+count_mat_tmp$G_A,
                         AT_TA = count_mat_tmp$A_T+count_mat_tmp$T_A,
                         AC_TG = count_mat_tmp$A_C+count_mat_tmp$T_G,
                         AG_TC = count_mat_tmp$A_G+count_mat_tmp$T_C,
                         insertion = count_mat_tmp$A_N+count_mat_tmp$G_N+count_mat_tmp$T_N+count_mat_tmp$C_N,
                         deletion = count_mat_tmp$N_A+count_mat_tmp$N_G+count_mat_tmp$N_T+count_mat_tmp$N_C
)


count_mat_new3 = as.matrix(count_mat_new)
group_total = apply(count_mat_new3, 1, sum)

count_effectRAD5 = data.frame(RM = as.vector(count_mat_new3[1,]+count_mat_new3[2,]),
                              BY = as.vector(count_mat_new3[3,]+count_mat_new3[4,]),
                              ID = c("CG_GC", "CA_GT", "CT_GA", "AT_TA","AC_TG","AG_TC","insertion","deletion") )

## "CG_GC" "CA_GT" "CT_GA" "AT_TA" "AC_TG" "AG_TC" "insertion" "deletion"
group_total = as.numeric(colSums(count_effectRAD5[,1:2]))
rad5.p.val = c()
for (i in 1:8){
  rad5.p.val[i] = prop.test(as.numeric(count_effectRAD5[i,1:2]), group_total)$p.value
}

#MKT1 effect
count_effectMKT1 = data.frame(RM = as.vector(count_mat_new3[1,]+count_mat_new3[3,]),
                              BY = as.vector(count_mat_new3[2,]+count_mat_new3[4,]),
                              ID = c("CG_GC", "CA_GT", "CT_GA", "AT_TA","AC_TG","AG_TC","insertion","deletion") )

## "CG_GC" "CA_GT" "CT_GA" "AT_TA" "AC_TG" "AG_TC"
group_total = as.numeric(colSums(count_effectMKT1[,1:2]))

mkt1.p.val = c()
for (i in 1:8){
  mkt1.p.val[i] = prop.test(as.numeric(count_effectMKT1[i,1:2]), group_total)$p.value
}

### Record FDR 
all.p.val = c(rad5.p.val, mkt1.p.val)
p.adjust(all.p.val, method='fdr')


count_mat_new2 = count_mat_new/group_total
count_mat_new2$ID = count_mat$ID
df = melt(count_mat_new2)

#order x axis by mutation rate from low to high
p=ggplot(data=df, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") +
  scale_fill_brewer(palette = "Set3") + scale_x_discrete(limits=c('S3', 'S1', 'S4', 'S2') ) +
  geom_text(aes(label=round(value,2)), size =3, position = position_stack(vjust = 0.5)) +
  theme_bw()
print(p)


part = count_mat_new2[,1:6]/rowSums(count_mat_new2[,1:6])
part$ID = count_mat$ID
df2=melt(part)
p2=ggplot(data=df2, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") +
  scale_fill_brewer(palette = "Set3") + scale_x_discrete(limits=c('S3', 'S1', 'S4', 'S2') ) +
  geom_text(aes(label=round(value,2)), size =3, position = position_stack(vjust = 0.5)) +
  theme_bw()
print(p2)






##if want to make the plot to be ordered by value
a$ID=with(a, factor(ID,levels=a$ID[order(a$value)]))

#calculate the percentage (ID column cannot be added as numeric)
sam_sum=apply(count_mat[,2:ncol(count_mat)], 1, sum)
count_mat2 = count_mat[,2:ncol(count_mat)]/sam_sum
count_mat2=cbind(ID=count_mat$ID, count_mat2)

count_mat2=count_mat2[,selected_mutation_type]
###mutation.fil is pos-seq-ref  need to reverse the *_* order
type=names(count_mat2)
type_change=type
type_change[2:length(type)]=sapply(type[2:length(type)],function(x) paste(strsplit(x,'_')[[1]][2],strsplit(x,'_')[[1]][1],sep ="_"))

names(count_mat2)=type_change

#combine two types into one
count_mat3=data.frame(ID=count_mat$ID,
                      CG_GC = count_mat2$C_G+count_mat2$G_C,
                      CA_GT = count_mat2$C_A+count_mat2$G_T,
                      CT_GA = count_mat2$C_T+count_mat2$G_A,
                      AT_TA = count_mat2$A_T+count_mat2$T_A,
                      AC_TG = count_mat2$A_C+count_mat2$T_G,
                      AG_TC = count_mat2$A_G+count_mat2$T_C
                      )
#count_mat3 = cbind(count_mat3, CG_GC = count_mat2$C_G+count_mat2$G_C)


normmat3 = count_mat3[2:7]/rowSums(count_mat3[2:7])
normmat3$ID = count_mat3$ID

selected_mutation_type = grep("N", colnames(normmat3), invert=T)
df = melt(normmat3[,selected_mutation_type])

p=ggplot(data=df, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") +
  scale_fill_brewer(palette = "Set3") + 
  geom_text(aes(label=round(value,2)), size =3, position = position_stack(vjust = 0.5)) +
  theme_bw()
print(p)

#order x axis by mutation rate from low to high
p=ggplot(data=df, aes(x=ID, y=value, fill=variable)) + geom_bar(stat="identity", position="stack") +
  scale_fill_brewer(palette = "Set3") + scale_x_discrete(limits=c('S3', 'S1', 'S4', 'S2') ) +
  geom_text(aes(label=round(value,2)), size =3, position = position_stack(vjust = 0.5)) +
  theme_bw()
print(p)


#function
get_replicated_type  =function(s1, s2)
{
  S1 = unlist(mutation.fil[[s1]])
  S1 = S1[seq(1, length(S1),2)]
  S1_rep = unlist(mutation.fil[[s2]])
  S1_rep = S1_rep[seq(1, length(S1_rep),2)]
  S1_common = intersect(S1, S1_rep)
  #S1_common = union(S1, S1_rep)
  
  all_types = unique(
    unlist(lapply(S1_common, function(x) paste0(strsplit(x,"_")[[1]][2:3],collapse="_")
    )
    ))
  all_types = gsub("-", "N", all_types)
  return(list(this_type=all_types, S_common=S1_common))
}

expand_count = function(ID, S_common, count_mat)
{
this_vec = NULL
this_vec$ID = ID
for(j in 2:ncol(count_mat))
  this_vec[[colnames(count_mat)[j]]] = 0

for(i in 1:length(S_common))
{
  this_type = paste0(strsplit(S_common[i],"_")[[1]][2:3], collapse="_")
  this_type = gsub("-", "N", this_type)
  if(length(this_vec[[this_type]])<1)
  {
    this_vec[[this_type]]=1
  } else {
    this_vec[[this_type]] = this_vec[[this_type]] + 1
  }
}

count_mat =rbind.data.frame(count_mat, this_vec, stringsAsFactors=F)
return(count_mat)
}                   