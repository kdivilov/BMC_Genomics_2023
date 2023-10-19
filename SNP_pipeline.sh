mkdir families
mkdir families/30.004
mkdir families/30.004/raw
mkdir families/30.004/trim
mkdir families/30.004/map
mkdir families/30.058
mkdir families/30.058/raw
mkdir families/30.058/trim
mkdir families/30.058/map
mkdir families/30.062
mkdir families/30.062/raw
mkdir families/30.062/trim
mkdir families/30.062/map
mkdir families/30.065
mkdir families/30.065/raw
mkdir families/30.065/trim
mkdir families/30.065/map


#Cgigas_RefSeq.fasta is GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fasta (RefSeq GCF_902806645.1) from NCBI with linkage groups renamed to chromosomes according to Peñaloza et al. (2021) and unplaced scaffolds removed
bwa index Cgigas_RefSeq.fasta

#Download and place raw reads from NCBI BioProject PRJNA873124 into their respective families/30.0XX/raw folders

#Merge the reads of the parents of the four families, which were sequenced twice (except the male parent of family 30.004)
cat families/30.004/raw/30.004_F_1.fastq.gz families/30.004/raw/30.004_F_2.fastq.gz > families/30.004/raw/30.004_F.fastq.gz
mv families/30.004/raw/30.004_M_1.fastq.gz families/30.004/raw/30.004_M.fastq.gz
cat families/30.058/raw/30.058_F_1.fastq.gz families/30.058/raw/30.058_F_2.fastq.gz > families/30.058/raw/30.058_F.fastq.gz
cat families/30.058/raw/30.058_M_1.fastq.gz families/30.058/raw/30.058_M_2.fastq.gz > families/30.058/raw/30.058_M.fastq.gz
cat families/30.062/raw/30.062_F_1.fastq.gz families/30.062/raw/30.062_F_2.fastq.gz > families/30.062/raw/30.062_F.fastq.gz
cat families/30.062/raw/30.062_M_1.fastq.gz families/30.062/raw/30.062_M_2.fastq.gz > families/30.062/raw/30.062_M.fastq.gz
cat families/30.065/raw/30.065_F_1.fastq.gz families/30.065/raw/30.065_F_2.fastq.gz > families/30.065/raw/30.065_F.fastq.gz
cat families/30.065/raw/30.065_M_1.fastq.gz families/30.065/raw/30.065_M_2.fastq.gz > families/30.065/raw/30.065_M.fastq.gz

#Trim raw reads
sh families_trim.sh

#Map trimmed reads to the reference genome
sh families_map.sh

#Call biallelic SNPs
bcftools mpileup -f Cgigas_RefSeq.fasta families/30.004/map/*.bam -C 50 -Q 30 -q 40 -a "AD,DP" -I -Ou | bcftools call -vm -Ou | bcftools view -m2 -M2 -Oz -o 30.004_snps.vcf.gz
bcftools mpileup -f Cgigas_RefSeq.fasta families/30.058/map/*.bam -C 50 -Q 30 -q 40 -a "AD,DP" -I -Ou | bcftools call -vm -Ou | bcftools view -m2 -M2 -Oz -o 30.058_snps.vcf.gz
bcftools mpileup -f Cgigas_RefSeq.fasta families/30.062/map/*.bam -C 50 -Q 30 -q 40 -a "AD,DP" -I -Ou | bcftools call -vm -Ou | bcftools view -m2 -M2 -Oz -o 30.062_snps.vcf.gz
bcftools mpileup -f Cgigas_RefSeq.fasta families/30.065/map/*.bam -C 50 -Q 30 -q 40 -a "AD,DP" -I -Ou | bcftools call -vm -Ou | bcftools view -m2 -M2 -Oz -o 30.065_snps.vcf.gz

#Run the AlphaFamImpute_prep.R script to filter SNPs and generate the three input files per family necessary for AlphaFamImpute

#Phase and impute SNPs
AlphaFamImpute -seqfile seqfile_30.004.txt -pedigree pedigree_30.004.txt -map map_30.004.txt -gbs -calling_threshold 0.9 -out 30.004
AlphaFamImpute -seqfile seqfile_30.058.txt -pedigree pedigree_30.058.txt -map map_30.058.txt -gbs -calling_threshold 0.9 -out 30.058
AlphaFamImpute -seqfile seqfile_30.062.txt -pedigree pedigree_30.062.txt -map map_30.062.txt -gbs -calling_threshold 0.9 -out 30.062
AlphaFamImpute -seqfile seqfile_30.065.txt -pedigree pedigree_30.065.txt -map map_30.065.txt -gbs -calling_threshold 0.9 -out 30.065


