library(vcfR)

#Family 30.004

vcf = read.vcfR("30.004_snps.vcf.gz")

ad = extract.gt(vcf, element = "AD")
colnames(ad) = gsub("\\.bam*","",basename(colnames(ad)))
ad = t(ad)

gt = extract.gt(vcf, element = "GT")
colnames(gt) = gsub("\\.bam*","",basename(colnames(gt)))

gt[gt=="0/0"] = 0
gt[gt=="0/1"] = 1
gt[gt=="1/1"] = 2
gt = t(gt)
mode(gt) = "numeric"

dp = extract.gt(vcf, element = "DP")
colnames(dp) = gsub("\\.bam*","",basename(colnames(dp)))
dp = t(dp)
mode(dp) = "numeric"

gt[which(dp<5 | dp>100)] = NA
ad[which(dp<5 | dp>100)] = "0,0"

MAF = colMeans(gt,na.rm=T)/2
gt = gt[,-which(MAF<0.15 | MAF>0.85)]
ad = ad[,-which(MAF<0.15 | MAF>0.85)]

ad = ad[,which(colSums(is.na(gt))<nrow(ad)*0.8)]

seqfile = matrix(NA,nrow(ad)*2,ncol(ad))
for(i in 1:nrow(ad)){
  seqfile[seq(1,nrow(ad)*2,2)[i],] = sapply(strsplit(ad[i,],","), `[`, 1)
  seqfile[seq(1,nrow(ad)*2,2)[i]+1,] = sapply(strsplit(ad[i,],","), `[`, 2)
}
mode(seqfile) = "numeric"
rownames(seqfile) = rep(rownames(ad),each=2)

map = cbind(gsub("Chr","",sapply(strsplit(colnames(ad),"_"), `[`, 1)),sapply(strsplit(colnames(ad),"_"), `[`, 2))

pedigree = rbind(cbind(rownames(ad)[2],0,0),
                 cbind(rownames(ad)[1],0,0),
                 cbind(rownames(ad)[-c(1:2)],rownames(ad)[2],rownames(ad)[1]))

write.table(seqfile,"seqfile_30.004.txt",col.names = F,quote = F)
write.table(map,"map_30.004.txt",row.names = F,col.names = F, quote = F)
write.table(pedigree,"pedigree_30.004.txt",row.names = F,col.names = F, quote = F)
write.table(colnames(ad),"snpnames_30.004.txt",row.names = F,col.names = F,quote = F)


#Family 30.058

vcf = read.vcfR("30.058_snps.vcf.gz")

ad = extract.gt(vcf, element = "AD")
colnames(ad) = gsub("\\.bam*","",basename(colnames(ad)))
ad = t(ad)

gt = extract.gt(vcf, element = "GT")
colnames(gt) = gsub("\\.bam*","",basename(colnames(gt)))

gt[gt=="0/0"] = 0
gt[gt=="0/1"] = 1
gt[gt=="1/1"] = 2
gt = t(gt)
mode(gt) = "numeric"

dp = extract.gt(vcf, element = "DP")
colnames(dp) = gsub("\\.bam*","",basename(colnames(dp)))
dp = t(dp)
mode(dp) = "numeric"

gt[which(dp<5 | dp>100)] = NA
ad[which(dp<5 | dp>100)] = "0,0"

MAF = colMeans(gt,na.rm=T)/2
gt = gt[,-which(MAF<0.15 | MAF>0.85)]
ad = ad[,-which(MAF<0.15 | MAF>0.85)]

ad = ad[,which(colSums(is.na(gt))<nrow(ad)*0.8)]

seqfile = matrix(NA,nrow(ad)*2,ncol(ad))
for(i in 1:nrow(ad)){
  seqfile[seq(1,nrow(ad)*2,2)[i],] = sapply(strsplit(ad[i,],","), `[`, 1)
  seqfile[seq(1,nrow(ad)*2,2)[i]+1,] = sapply(strsplit(ad[i,],","), `[`, 2)
}
mode(seqfile) = "numeric"
rownames(seqfile) = rep(rownames(ad),each=2)

map = cbind(gsub("Chr","",sapply(strsplit(colnames(ad),"_"), `[`, 1)),sapply(strsplit(colnames(ad),"_"), `[`, 2))

pedigree = rbind(cbind(rownames(ad)[2],0,0),
                 cbind(rownames(ad)[1],0,0),
                 cbind(rownames(ad)[-c(1:2)],rownames(ad)[2],rownames(ad)[1]))

write.table(seqfile,"seqfile_30.058.txt",col.names = F,quote = F)
write.table(map,"map_30.058.txt",row.names = F,col.names = F, quote = F)
write.table(pedigree,"pedigree_30.058.txt",row.names = F,col.names = F, quote = F)
write.table(colnames(ad),"snpnames_30.058.txt",row.names = F,col.names = F,quote = F)


#Family 30.062

vcf = read.vcfR("30.062_snps.vcf.gz")

ad = extract.gt(vcf, element = "AD")
colnames(ad) = gsub("\\.bam*","",basename(colnames(ad)))
ad = t(ad)

gt = extract.gt(vcf, element = "GT")
colnames(gt) = gsub("\\.bam*","",basename(colnames(gt)))

gt[gt=="0/0"] = 0
gt[gt=="0/1"] = 1
gt[gt=="1/1"] = 2
gt = t(gt)
mode(gt) = "numeric"

dp = extract.gt(vcf, element = "DP")
colnames(dp) = gsub("\\.bam*","",basename(colnames(dp)))
dp = t(dp)
mode(dp) = "numeric"

gt[which(dp<5 | dp>100)] = NA
ad[which(dp<5 | dp>100)] = "0,0"

MAF = colMeans(gt,na.rm=T)/2
gt = gt[,-which(MAF<0.15 | MAF>0.85)]
ad = ad[,-which(MAF<0.15 | MAF>0.85)]

ad = ad[,which(colSums(is.na(gt))<nrow(ad)*0.8)]

seqfile = matrix(NA,nrow(ad)*2,ncol(ad))
for(i in 1:nrow(ad)){
  seqfile[seq(1,nrow(ad)*2,2)[i],] = sapply(strsplit(ad[i,],","), `[`, 1)
  seqfile[seq(1,nrow(ad)*2,2)[i]+1,] = sapply(strsplit(ad[i,],","), `[`, 2)
}
mode(seqfile) = "numeric"
rownames(seqfile) = rep(rownames(ad),each=2)

map = cbind(gsub("Chr","",sapply(strsplit(colnames(ad),"_"), `[`, 1)),sapply(strsplit(colnames(ad),"_"), `[`, 2))

pedigree = rbind(cbind(rownames(ad)[2],0,0),
                 cbind(rownames(ad)[1],0,0),
                 cbind(rownames(ad)[-c(1:2)],rownames(ad)[2],rownames(ad)[1]))

write.table(seqfile,"seqfile_30.062.txt",col.names = F,quote = F)
write.table(map,"map_30.062.txt",row.names = F,col.names = F, quote = F)
write.table(pedigree,"pedigree_30.062.txt",row.names = F,col.names = F, quote = F)
write.table(colnames(ad),"snpnames_30.062.txt",row.names = F,col.names = F,quote = F)


#Family 30.065

vcf = read.vcfR("30.065_snps.vcf.gz")

ad = extract.gt(vcf, element = "AD")
colnames(ad) = gsub("\\.bam*","",basename(colnames(ad)))
ad = t(ad)

gt = extract.gt(vcf, element = "GT")
colnames(gt) = gsub("\\.bam*","",basename(colnames(gt)))

gt[gt=="0/0"] = 0
gt[gt=="0/1"] = 1
gt[gt=="1/1"] = 2
gt = t(gt)
mode(gt) = "numeric"

dp = extract.gt(vcf, element = "DP")
colnames(dp) = gsub("\\.bam*","",basename(colnames(dp)))
dp = t(dp)
mode(dp) = "numeric"

gt[which(dp<5 | dp>100)] = NA
ad[which(dp<5 | dp>100)] = "0,0"

MAF = colMeans(gt,na.rm=T)/2
gt = gt[,-which(MAF<0.15 | MAF>0.85)]
ad = ad[,-which(MAF<0.15 | MAF>0.85)]

ad = ad[,which(colSums(is.na(gt))<nrow(ad)*0.8)]

seqfile = matrix(NA,nrow(ad)*2,ncol(ad))
for(i in 1:nrow(ad)){
  seqfile[seq(1,nrow(ad)*2,2)[i],] = sapply(strsplit(ad[i,],","), `[`, 1)
  seqfile[seq(1,nrow(ad)*2,2)[i]+1,] = sapply(strsplit(ad[i,],","), `[`, 2)
}
mode(seqfile) = "numeric"
rownames(seqfile) = rep(rownames(ad),each=2)

map = cbind(gsub("Chr","",sapply(strsplit(colnames(ad),"_"), `[`, 1)),sapply(strsplit(colnames(ad),"_"), `[`, 2))

pedigree = rbind(cbind(rownames(ad)[2],0,0),
                 cbind(rownames(ad)[1],0,0),
                 cbind(rownames(ad)[-c(1:2)],rownames(ad)[2],rownames(ad)[1]))

write.table(seqfile,"seqfile_30.065.txt",col.names = F,quote = F)
write.table(map,"map_30.065.txt",row.names = F,col.names = F, quote = F)
write.table(pedigree,"pedigree_30.065.txt",row.names = F,col.names = F, quote = F)
write.table(colnames(ad),"snpnames_30.065.txt",row.names = F,col.names = F,quote = F)