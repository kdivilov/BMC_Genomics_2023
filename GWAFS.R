#Function to obtain Fisher's exact test p-values from pre-planting and post-mortality samples
fishers_exact_test = function(pre,post){
  pval=c()
  for(i in 1:ncol(pre)){
    pval = c(pval,fisher.test(rbind(table(pre[,i]),table(post[i])))$p.value)
  }
  return(pval)
}


#Family 30.004

phase = read.table("30.004.phase")
snpnames = read.table("snpnames_30.004.txt")

phase = phase[-c(1:4),]
phase[phase==9] = NA

paternal = phase[seq(1,nrow(phase),2),]
colnames(paternal) = snpnames[,1]
rownames(paternal) = paternal[,1]
paternal = paternal[,-1]
paternal = paternal[,-which(colMeans(paternal,na.rm=T)==0 | colMeans(paternal,na.rm=T)==1 | is.na(colMeans(paternal,na.rm=T)))]
paternal = paternal[,-which((colSums(is.na(paternal))/nrow(paternal))>0.1)]

paternal_pre = paternal[grep("pre",rownames(paternal)),]
paternal_post = paternal[grep("post",rownames(paternal)),]
paternal_post_rep1 = paternal[grep("post_rep1",rownames(paternal)),]
paternal_post_rep2 = paternal[grep("post_rep2",rownames(paternal)),]
paternal_post_rep3 = paternal[grep("post_rep3",rownames(paternal)),]

#GWAFS for paternal haplotype
paternal_pre_post = fishers_exact_test(paternal_pre,paternal_post)
paternal_pre_post_rep1 = fishers_exact_test(paternal_pre,paternal_post_rep1)
paternal_pre_post_rep2 = fishers_exact_test(paternal_pre,paternal_post_rep2)
paternal_pre_post_rep3 = fishers_exact_test(paternal_pre,paternal_post_rep3)


maternal = phase[seq(2,nrow(phase),2),]
colnames(maternal) = snpnames[,1]
rownames(maternal) = maternal[,1]
maternal = maternal[,-1]
maternal = maternal[,-which(colMeans(maternal,na.rm=T)==0 | colMeans(maternal,na.rm=T)==1 | is.na(colMeans(maternal,na.rm=T)))]
maternal = maternal[,-which((colSums(is.na(maternal))/nrow(maternal))>0.1)]

maternal_pre = maternal[grep("pre",rownames(maternal)),]
maternal_post = maternal[grep("post",rownames(maternal)),]
maternal_post_rep1 = maternal[grep("post_rep1",rownames(maternal)),]
maternal_post_rep2 = maternal[grep("post_rep2",rownames(maternal)),]
maternal_post_rep3 = maternal[grep("post_rep3",rownames(maternal)),]

#GWAFS for maternal haplotype
maternal_pre_post = fishers_exact_test(maternal_pre,maternal_post)
maternal_pre_post_rep1 = fishers_exact_test(maternal_pre,maternal_post_rep1)
maternal_pre_post_rep2 = fishers_exact_test(maternal_pre,maternal_post_rep2)
maternal_pre_post_rep3 = fishers_exact_test(maternal_pre,maternal_post_rep3)


#Family 30.058

phase = read.table("30.058.phase")
snpnames = read.table("snpnames_30.058.txt")

phase = phase[-c(1:4),]
phase[phase==9] = NA

paternal = phase[seq(1,nrow(phase),2),]
colnames(paternal) = snpnames[,1]
rownames(paternal) = paternal[,1]
paternal = paternal[,-1]
paternal = paternal[,-which(colMeans(paternal,na.rm=T)==0 | colMeans(paternal,na.rm=T)==1 | is.na(colMeans(paternal,na.rm=T)))]
paternal = paternal[,-which((colSums(is.na(paternal))/nrow(paternal))>0.1)]

paternal_pre = paternal[grep("pre",rownames(paternal)),]
paternal_post = paternal[grep("post",rownames(paternal)),]
paternal_post_rep1 = paternal[grep("post_rep1",rownames(paternal)),]
paternal_post_rep2 = paternal[grep("post_rep2",rownames(paternal)),]
paternal_post_rep3 = paternal[grep("post_rep3",rownames(paternal)),]

#GWAFS for paternal haplotype
paternal_pre_post = fishers_exact_test(paternal_pre,paternal_post)
paternal_pre_post_rep1 = fishers_exact_test(paternal_pre,paternal_post_rep1)
paternal_pre_post_rep2 = fishers_exact_test(paternal_pre,paternal_post_rep2)
paternal_pre_post_rep3 = fishers_exact_test(paternal_pre,paternal_post_rep3)


maternal = phase[seq(2,nrow(phase),2),]
colnames(maternal) = snpnames[,1]
rownames(maternal) = maternal[,1]
maternal = maternal[,-1]
maternal = maternal[,-which(colMeans(maternal,na.rm=T)==0 | colMeans(maternal,na.rm=T)==1 | is.na(colMeans(maternal,na.rm=T)))]
maternal = maternal[,-which((colSums(is.na(maternal))/nrow(maternal))>0.1)]

maternal_pre = maternal[grep("pre",rownames(maternal)),]
maternal_post = maternal[grep("post",rownames(maternal)),]
maternal_post_rep1 = maternal[grep("post_rep1",rownames(maternal)),]
maternal_post_rep2 = maternal[grep("post_rep2",rownames(maternal)),]
maternal_post_rep3 = maternal[grep("post_rep3",rownames(maternal)),]

#GWAFS for maternal haplotype
maternal_pre_post = fishers_exact_test(maternal_pre,maternal_post)
maternal_pre_post_rep1 = fishers_exact_test(maternal_pre,maternal_post_rep1)
maternal_pre_post_rep2 = fishers_exact_test(maternal_pre,maternal_post_rep2)
maternal_pre_post_rep3 = fishers_exact_test(maternal_pre,maternal_post_rep3)


#Family 30.062

phase = read.table("30.062.phase")
snpnames = read.table("snpnames_30.062.txt")

phase = phase[-c(1:4),]
phase[phase==9] = NA

paternal = phase[seq(1,nrow(phase),2),]
colnames(paternal) = snpnames[,1]
rownames(paternal) = paternal[,1]
paternal = paternal[,-1]
paternal = paternal[,-which(colMeans(paternal,na.rm=T)==0 | colMeans(paternal,na.rm=T)==1 | is.na(colMeans(paternal,na.rm=T)))]
paternal = paternal[,-which((colSums(is.na(paternal))/nrow(paternal))>0.1)]

paternal_pre = paternal[grep("pre",rownames(paternal)),]
paternal_post = paternal[grep("post",rownames(paternal)),]
paternal_post_rep1 = paternal[grep("post_rep1",rownames(paternal)),]
paternal_post_rep2 = paternal[grep("post_rep2",rownames(paternal)),]
paternal_post_rep3 = paternal[grep("post_rep3",rownames(paternal)),]

#GWAFS for paternal haplotype
paternal_pre_post = fishers_exact_test(paternal_pre,paternal_post)
paternal_pre_post_rep1 = fishers_exact_test(paternal_pre,paternal_post_rep1)
paternal_pre_post_rep2 = fishers_exact_test(paternal_pre,paternal_post_rep2)
paternal_pre_post_rep3 = fishers_exact_test(paternal_pre,paternal_post_rep3)


maternal = phase[seq(2,nrow(phase),2),]
colnames(maternal) = snpnames[,1]
rownames(maternal) = maternal[,1]
maternal = maternal[,-1]
maternal = maternal[,-which(colMeans(maternal,na.rm=T)==0 | colMeans(maternal,na.rm=T)==1 | is.na(colMeans(maternal,na.rm=T)))]
maternal = maternal[,-which((colSums(is.na(maternal))/nrow(maternal))>0.1)]

maternal_pre = maternal[grep("pre",rownames(maternal)),]
maternal_post = maternal[grep("post",rownames(maternal)),]
maternal_post_rep1 = maternal[grep("post_rep1",rownames(maternal)),]
maternal_post_rep2 = maternal[grep("post_rep2",rownames(maternal)),]
maternal_post_rep3 = maternal[grep("post_rep3",rownames(maternal)),]

#GWAFS for maternal haplotype
maternal_pre_post = fishers_exact_test(maternal_pre,maternal_post)
maternal_pre_post_rep1 = fishers_exact_test(maternal_pre,maternal_post_rep1)
maternal_pre_post_rep2 = fishers_exact_test(maternal_pre,maternal_post_rep2)
maternal_pre_post_rep3 = fishers_exact_test(maternal_pre,maternal_post_rep3)


#Family 30.065

phase = read.table("30.065.phase")
snpnames = read.table("snpnames_30.065.txt")

phase = phase[-c(1:4),]
phase[phase==9] = NA

paternal = phase[seq(1,nrow(phase),2),]
colnames(paternal) = snpnames[,1]
rownames(paternal) = paternal[,1]
paternal = paternal[,-1]
paternal = paternal[,-which(colMeans(paternal,na.rm=T)==0 | colMeans(paternal,na.rm=T)==1 | is.na(colMeans(paternal,na.rm=T)))]
paternal = paternal[,-which((colSums(is.na(paternal))/nrow(paternal))>0.1)]

paternal_pre = paternal[grep("pre",rownames(paternal)),]
paternal_post = paternal[grep("post",rownames(paternal)),]
paternal_post_rep1 = paternal[grep("post_rep1",rownames(paternal)),]
paternal_post_rep2 = paternal[grep("post_rep2",rownames(paternal)),]
paternal_post_rep3 = paternal[grep("post_rep3",rownames(paternal)),]

#GWAFS for paternal haplotype
paternal_pre_post = fishers_exact_test(paternal_pre,paternal_post)
paternal_pre_post_rep1 = fishers_exact_test(paternal_pre,paternal_post_rep1)
paternal_pre_post_rep2 = fishers_exact_test(paternal_pre,paternal_post_rep2)
paternal_pre_post_rep3 = fishers_exact_test(paternal_pre,paternal_post_rep3)


maternal = phase[seq(2,nrow(phase),2),]
colnames(maternal) = snpnames[,1]
rownames(maternal) = maternal[,1]
maternal = maternal[,-1]
maternal = maternal[,-which(colMeans(maternal,na.rm=T)==0 | colMeans(maternal,na.rm=T)==1 | is.na(colMeans(maternal,na.rm=T)))]
maternal = maternal[,-which((colSums(is.na(maternal))/nrow(maternal))>0.1)]

maternal_pre = maternal[grep("pre",rownames(maternal)),]
maternal_post = maternal[grep("post",rownames(maternal)),]
maternal_post_rep1 = maternal[grep("post_rep1",rownames(maternal)),]
maternal_post_rep2 = maternal[grep("post_rep2",rownames(maternal)),]
maternal_post_rep3 = maternal[grep("post_rep3",rownames(maternal)),]

#GWAFS for maternal haplotype
maternal_pre_post = fishers_exact_test(maternal_pre,maternal_post)
maternal_pre_post_rep1 = fishers_exact_test(maternal_pre,maternal_post_rep1)
maternal_pre_post_rep2 = fishers_exact_test(maternal_pre,maternal_post_rep2)
maternal_pre_post_rep3 = fishers_exact_test(maternal_pre,maternal_post_rep3)