# Main code of the pipeline for NIPT data GWAS

# fq2bam
```
#fastp -q 20 -u 10 -n 1 --in1 fq.gz -o clean.fq.gz
SOAPnuke filter -n 0.1 -q 0.1 -l 5 -G --seqType 1 -Q 2 --fq1 fq.gz -o ./ -C clean.fq.gz
#35bp
bwa aln GRCh38.fa clean.fq.gz|bwa samse GRCh38.fa - clean.fq.gz|samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20|samtools view -Sb -F 2820 -q 1 -|sambamba sort -t 4 -m 2G --tmpdir=./tmp -o sample.bam /dev/stdin
#135bp
bwa mem GRCh38.fa clean.fq.gz|samblaster --excludeDups --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20|samtools view -Sb -F 2820 -q 1 -|sambamba sort -t 4 -m 2G --tmpdir=./tmp -o sample.bam /dev/stdin
mosdepth qc sample.bam -f GRCh38.fa --fast-mode --no-per-base --by 1000000 --thresholds 1,2,4
java -Xmx3g -jar gatk.jar BaseRecalibrator -R GRCh38.fa -I sample.bam -O bqsr.file --known-sites dbsnp_146.hg38.vcf.gz --known-sites 1000G_omni2.5.hg38.vcf.gz --known-sites hapmap_3.3.hg38.vcf.gz
java -Xmx3g -jar gatk.jar ApplyBQSR -R GRCh38.fa -I sample.bam -O sample.bqsr.bam --bqsr-recal-file bqsr.file
```

# bam2SNV
```
ls path2bam/*.bqsr.bam > bam.list
BaseVar basetype -t 10 -b 100 -i bam.list -r GRCh38.fa -o sample
zcat sample.vcf.gz|cut -f '1-8'|bcftools norm - -m -|bcftools view - -i "QUAL > 100 && FORMAT/CM_DP > 100 && FORMAT/CM_AF > 0.01 && FORMAT/CM_AF < 0.99" -T 35bp.bed -O z -o sample.filter.vcf.gz
java -Xmx30g -jar gatk.jar VariantRecalibrator --tmp-dir ./tmp -R GRCh38.fa -V sample.filter.vcf.gz \
             -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
             -resource:omini,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
             -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
             -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz \
             --max-attempts 2 --maximum-training-variants 5000000 \
             -an BaseQRankSum -an ReadPosRankSum -an FS -an QD -an SOR -an MQRankSum -an CM_DP \
             -mode SNP \
             -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 -tranche 89.0 -tranche 88.0 -tranche 87.0 -tranche 86.0 -tranche 85.0 -tranche 80.0 \
             --rscript-file vqsr.R \
             --tranches-file vqsr.tranches \
            -O vqsr.recal
java -Xmx30g -jar gatk.jar ApplyVQSR --tmp-dir ./tmp -R GRCh38.fa -V sample.filter.vcf.gz \
            -ts-filter-level 99.0 \
            --tranches-file vqsr.tranches \
            --recal-file vqsr.recal \
            -mode SNP \
            -O sample.filter.vqsr.vcf.gz
zcat sample.filter.vqsr.vcf.gz|vawk --header '{if($7=="PASS")print}'|bgzip -f > sample.pass.vcf.gz
bcftools stats sample.pass.vcf.gz | grep "TSTV"
```

# SNP imputation
```
#reference panel
cat bam.list|awk -F '/' '{print $NF}'|awk -F '.' '{print $1}' > sample.id
tabix Han.vcf.gz chr1:1-10000000 -h|bcftools convert - --haplegendsample ref
tabix sample.pass.vcf.gz chr1:1-10000000|cut -f '1,2,4,5' > pos.txt
#QUILT for sample size < 10000
Rscript QUILT.R  --outputdir=.  --tempdir=./tmp  --sampleNames_file=sample.id  --chr=chr1  --regionStart=1  --regionEnd=10000000  --buffer=250000  --bamlist=bam.list  --posfile=pos.txt  --reference_haplotype_file=ref.hap.gz  --reference_legend_file=ref.legend.gz  --nGen=100  --nCores=20
zcat quilt.chr1.1.10000000.vcf.gz|vawk '{if(I$INFO_SCORE > 0.8 && I$EAF > 0.025 && I$EAF < 0.975 && I$HWE > 0.0005){print}}' --header|bgzip -f > chr1.1.10000000.SNP.vcf.gz
#STITCH for sample size > 10000
Rscript STITCH.R  --outputdir=.  --tempdir=./tmp  --sampleNames_file=sample.id  --chr=chr1  --regionStart=1  --regionEnd=10000000  --buffer=250000  --bamlist=bam.list  --posfile=pos.txt  --reference_haplotype_file=ref.hap.gz  --reference_legend_file=ref.legend.gz  --nGen=100  --nCores=20
zcat stitch.chr1.1.10000000.vcf.gz|vawk '{if(I$INFO_SCORE > 0.8 && I$EAF > 0.025 && I$EAF < 0.975 && I$HWE > 0.0005 && I$AC/I$AN > 0.8){print}}' --header|bgzip -f > chr1.1.10000000.SNP.vcf.gz
```

# VCF concat
```
ls path2vcf/*.SNP.vcf.gz|awk -F '/' '{print $NF}'|tr '.' ' '|sort -k1,1 -k2,2n|tr ' ' '.'|awk '{print "path2vcf/"$0}' > vcf.list
bcftools concat -f vcf.list -O z -o SNP.vcf.gz
tabix -s1 -b2 -e2 -f SNP.vcf.gz
```

# PCA
```
zcat SNP.vcf.gz |vawk --header '{if(I$EAF>0.05 && I$EAF<0.95 && NR % 10 == 1)print $1,$2,$3,$4,$5,$6,$7,$8,"GT",S$*$GT}'|bgzip -f > pca.vcf.gz
plink --vcf pca.vcf.gz --make-bed --out pca
flashpca --bfile pca --outpc pca.pc --batch 10
```

# GWAS
```
zcat SNP.vcf.gz|vawk '{print $1"_"$2"_"$4"_"$5,S$*$GT}'|bgzip -f > gwas.data.gz
Rscript gwas.r
```

# GWAS R code [gwas.r]
```
gdf<-read.table(gzfile('gwas.data.gz'))
fam<-read.table('gwas.fam',head=T)
pc<-read.table('pca.pc',head=T,head=T)
#names(fam)
# [1] "FID"            "IID"            "Age"            "Height"         "BMI"            "Birthweight"           "PPH"
#names(pc)
# [1] "FID"            "IID"            "PC1"            "PC2"         "PC3"            "PC4"           "PC5"
y<-fam$PPH

sdf<-as.data.frame(t(apply(gdf[,-1],1,function(x){
  x<-as.numeric(as.character(factor(x,levels=c('0|0','0|1','1|0','1|1'),labels=c(0,1,1,2))))
  coef<-summary(lm(y~fam$Age+fam$Height+fam$BMI+fam$Birthweight+pc$PC1+pc$PC2+pc$PC3+pc$PC4+pc$PC5+x))$coefficients
  return(coef[nrow(coef),c(1,2,4)])
})))
names(sdf)<-c('Beta','se','Pvalue')
sdf$SNP<-gdf[,1]
write.table(sdf[,c(4,1,2,3)],file=gzfile('gwas.summary.gz'),sep=' ',quote=F,row.names=F,col.names=T)
```

# Heritability
```
zcat gwas.summary.gz|awk '$NF<0.01{print $1}'|awk -F '_' '{print $1,$2-1,$2}'|tr ' ' '\t' > s.bed
tabix SNP.vcf.gz -R s.bed -h|plink --vcf - --make-bed --out SNP
gcta --bfile SNP --ld-score-region 500 --out SNP
gcta --bfile SNP --make-grm --out sample
Rscript split.ld.r
gcta --bfile sample --extract SNP_group1.txt --make-grm --out sample_group1
gcta --bfile sample --extract SNP_group2.txt --make-grm --out sample_group2
gcta --bfile sample --extract SNP_group3.txt --make-grm --out sample_group3
gcta --bfile sample --extract SNP_group4.txt --make-grm --out sample_group4
#echo "sample_group1 sample_group2 sample_group3 sample_group4"|tr ' ' '\n' > multi_GRMs.txt
gcta --reml --mgrm multi_GRMs.txt --pheno phe.txt --qcovar cov.txt --mpheno 1 --out H2
```
# R code [split.ld.r]
```
lds_seg = read.table("SNP.score.ld",header=T,colClasses=c("character",rep("numeric",8)))
quartiles=summary(lds_seg$ldscore_SNP)
lb1 = which(lds_seg$ldscore_SNP <= quartiles[2])
lb2 = which(lds_seg$ldscore_SNP > quartiles[2] & lds_seg$ldscore_SNP <= quartiles[3])
lb3 = which(lds_seg$ldscore_SNP > quartiles[3] & lds_seg$ldscore_SNP <= quartiles[5])
lb4 = which(lds_seg$ldscore_SNP > quartiles[5])
lb1_snp = lds_seg$SNP[lb1]
lb2_snp = lds_seg$SNP[lb2]
lb3_snp = lds_seg$SNP[lb3]
lb4_snp = lds_seg$SNP[lb4]
write.table(lb1_snp, "SNP_group1.txt", row.names=F, quote=F, col.names=F)
write.table(lb2_snp, "SNP_group2.txt", row.names=F, quote=F, col.names=F)
write.table(lb3_snp, "SNP_group3.txt", row.names=F, quote=F, col.names=F)
write.table(lb4_snp, "SNP_group4.txt", row.names=F, quote=F, col.names=F)'
```

# Enrichment R code
```
library(EWCE)
PPH_gene<-c("ADAMTS20","ERVFC1","GATA6","HTR1B","KIF26A","NEGR1","RERG","SEMA6A","TMC1")
#load('CellTypeData_muscle.rda')
load('CellTypeData_uterus.rda')

bg  = attr(ctd[[1]]$specificity,'dimnames')[[1]]
hits=unique(PPH_gene[PPH_gene %in% bg])
reps=10000
set.seed(2022)
full_results = bootstrap_enrichment_test(sct_data=ctd,
                                         hits=hits,
                                         bg=bg,
                                         reps=reps,
                                         annotLevel=1)
full_results$results$padjust<-p.adjust(full_results$results$p,method='fdr')
```

# GGM Network r code
```
library(huge)
library(GGally)
library(network)
library(ggplot2)

df<-read.table('gene.exp',head=T)
df<-df[,colSums(df)>0]

mbModel <- huge(as.matrix(df), method = "mb")
mbOptRIC = huge.select(mbModel)
mbOptRICGraph = mbOptRIC$refit
net<-network(mbOptRICGraph)
network.vertex.names(net)<-names(df)

tmp<-as.matrix(sparseMatrix(i=mbOptRIC$refit@i+1,p=mbOptRIC$refit@p,x=mbOptRIC$refit@x,dims=c(32,32)),diag=T)
network.edgecount(net)

net %v% "color" = 'black'
cor.matrix<- cor(df,use = "pairwise.complete.obs")
net %e% "weights"<-abs(cor.matrix[which(tmp==1 & lower.tri(tmp),arr.ind=TRUE)])*20

set.seed(1)
ggnet2(net,
       size = 10,
       hjust=-0.4,
       label = TRUE,
       label.size=4,
       edge.size = "weights",
       color= "color",
       layout.exp=0.8)
```

# Choose reference panel
```
# 20130606_g1k_3202_samples_ped_population.txt
# 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
# cram files
# find these files here
# https://42basepairs.com/browse/s3/1000genomes/1000G_2504_high_coverage

awk '$6=="CDX" && $3==0 && $4==0{print $2}' 20130606_g1k_3202_samples_ped_population.txt > Dai.list
awk '($6=="CHS" || $6=="CHB") && $3==0 && $4==0{print $2}' 20130606_g1k_3202_samples_ped_population.txt > Han.list
awk '$6!="CDX" && $3==0 && $4==0{print $2}' 20130606_g1k_3202_samples_ped_population.txt > global.list
bcftools view 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -S Dai.list -r chr21:10000001-20000000|vawk --header '{if(length($4)==1 && length($5)==1 && I$AC/I$AN > 0.01 && I$AC/I$AN < 0.99)print}'|bgzip -f > Dai.vcf.gz
bcftools view 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -S Han.list -r chr21:10000001-20000000|vawk --header '{if(length($4)==1 && length($5)==1 && I$AC/I$AN > 0.01 && I$AC/I$AN < 0.99)print}'|bgzip -f > Han.vcf.gz
bcftools view 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -S global.list -r chr21:10000001-20000000|vawk --header '{if(length($4)==1 && length($5)==1 && I$AC/I$AN > 0.01 && I$AC/I$AN < 0.99)print}'|bgzip -f > global.vcf.gz
for i in `cat Dai.list`
do
samtools view path2cram/$i.cram -T GRCh38.fa chr21:10000001-20000000 --subsample 0.015 --cram -o $i.downsample.cram
samtools index $i.downsample.cram
echo "`pwd`/$i.downsample.cram"
done > cram.list
zcat global.vcf.gz|vawk '{print $1,$2,$4,$5}' > pos.txt
`date`
bcftools convert Han.vcf.gz --haplegendsample ref
zcat Han.vcf.gz|vawk '{print $1,$2,$4,$5}'|cut -f '1,2,4,5' > ref.pos
Rscript QUILT.R  --outputdir=./Han_ref --tempdir=./Han_ref/tmp  --sampleNames_file=Dai.list  --chr=chr21  --regionStart=10250001  --regionEnd=19750000  --buffer=250000 --bamlist=cram.list --posfile=pos.txt  --reference_haplotype_file=ref.hap.gz  --reference_legend_file=ref.legend.gz  --nGen=100  --nCores=20 --reference=GRCh38.fa
`date`
bcftools convert global.vcf.gz --haplegendsample ref
zcat global.vcf.gz|vawk '{print $1,$2,$4,$5}'|cut -f '1,2,4,5' > ref.pos
Rscript QUILT.R  --outputdir=./global_ref --tempdir=./global_ref/tmp  --sampleNames_file=Dai.list  --chr=chr21  --regionStart=10250001  --regionEnd=19750000  --buffer=250000 --bamlist=cram.list --posfile=pos.txt  --reference_haplotype_file=ref.hap.gz  --reference_legend_file=ref.legend.gz  --nGen=100  --nCores=20 --reference=GRCh38.fa
`date`
zcat ./Han_ref/quilt.chr21.10250001.19750000.vcf.gz|vawk '{if($2>10250001 && $2<19750000 && I$INFO_SCORE> 0.8 && I$HWE >0.0005 && I$EAF>0.025 && I$EAF<0.975){print $1"_"$2"_"$4"_"$5,I$EAF,S$*$GT}}'|sort -k1,1|gzip -f > a.gz
zcat ./global_ref/quilt.chr21.10250001.19750000.vcf.gz|vawk '{if($2>10250001 && $2<19750000 && I$INFO_SCORE> 0.8 && I$HWE >0.0005 && I$EAF>0.025 && I$EAF<0.975){print $1"_"$2"_"$4"_"$5,I$EAF,S$*$GT}}'|sort -k1,1|gzip -f > b.gz
bcftools view Dai.vcf.gz -r chr21:10250001-19750000 |vawk '{print $1"_"$2"_"$4"_"$5,I$AC/I$AN,S$*$GT}'|sort -k1,1|gzip -f > c.gz
zcat a.gz|cut -f 1 > a.pos
zcat b.gz|cut -f 1 > b.pos
wc -l a.pos
wc -l b.pos
zcat c.gz|cut -f 1 > c.pos
join a.pos c.pos > ac.pos
join b.pos c.pos > bc.pos
zcat a.gz|join - ac.pos|gzip -f > ac.gz
zcat b.gz|join - bc.pos|gzip -f > bc.gz
zcat c.gz|join - ac.pos|gzip -f > ca.gz
zcat c.gz|join - bc.pos|gzip -f > cb.gz
Rscript accuracy.r
```

# R code [accuracy.r]
```
ac<-read.table(gzfile('ac.gz'))
bc<-read.table(gzfile('bc.gz'))
ca<-read.table(gzfile('ca.gz'))
cb<-read.table(gzfile('cb.gz'))
for(i in 3:ncol(ac)){
  ac[,i]<-as.numeric(as.character(factor(ac[,i],levels=c('0|0','0|1','1|0','1|1'),labels=c(0,1,1,2))))
  bc[,i]<-as.numeric(as.character(factor(bc[,i],levels=c('0|0','0|1','1|0','1|1'),labels=c(0,1,1,2))))
  ca[,i]<-as.numeric(as.character(factor(ca[,i],levels=c('0|0','0|1','1|0','1|1'),labels=c(0,1,1,2))))
  cb[,i]<-as.numeric(as.character(factor(cb[,i],levels=c('0|0','0|1','1|0','1|1'),labels=c(0,1,1,2))))
}
accuracy_Han<-sum(ac[,3:nn]==ca[,3:nn])/(nrow(ac)*(ncol(ac)-2))
accuracy_global<-sum(bc[,3:nn]==cb[,3:nn])/(nrow(bc)*(ncol(bc)-2))
print(c(accuracy_Han,accuracy_globa))
```
