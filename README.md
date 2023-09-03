## WGS 流程
## step 1 去接头比对

```
ls *_1.fq.gz|while read id;
do
echo "fastp -q 20 -u 50 -l 35 -n 15 -i ${id} ${id%%_*}_2.fq.gz -o ../01.clean_data/${id%%_*}_1.fastq.gz   ../01.clean_data/${id%%_*}_2.fastq.gz --cut_tail --n_base_limit 3 --length_required 60 --correction --thread 4 " >>setp1.sh;
done

cd 01.clean_data
ls *_clean_1.fastq.gz|while read id;
do
echo "bwa mem -t 12  -k 32 -M -R '@RG\tID:foo\tPL:illumina\tSM:${id%%_*}' \
/public/agis/liuyuwen_group/leibowen/ref/camel_new/camel.fa ${id} ${id%%_*}_clean_2.fastq.gz | samtools view -Sb - > ../02.bwa/${id%%_*}.unsorted.bam" >> setp2_bwa.sh;
done
```
## step 2 排序，标记重复
```
ls *.sorted.bam >last
#ls *.sorted.bam|cut -d "." -f 1 >last 
cat last | while read id
do
name=${id%%.*}
echo "samtools sort -@ 6 -m 4G -O bam -o ${name}.sorted.bam ${id}" >> bam.sort.sh
done

cat last | while read id
do
name=${id%%.*}
echo "java -jar /public/home/leibowen/software/picard/build/libs/picard.jar MarkDuplicates I=${id} O=${id%%.*}.markdup.bam M=${id%%.*}.markdup.metrc.csv" >> markdup.commandlines.sh
done
```
## step 3 按染色体callg vcf
```
1每个个体call gvcf
ls ../02.bwa/*.markdup.bam|while read id ;do
echo "gatk --java-options "-Xmx15g" HaplotypeCaller -R $ref -I ../02.bwa/${id} -O ${id%%.*}.g.vcf.gz -ERC GVCF" >>gvcf.sh
done

2按染色体合并GCVF 
ls *.g.vcf >GVCF.list
cat chr|while read id;
do
echo "gatk --java-options -Xmx20g CombineGVCFs -R  $ref --variant GVCF.list -L ${id} -O ./VCF/${id}.g.vcf.gz" >>merge_gVCF.sh;
done

3 每个染色体基因分型
ls *.g.vcf.gz|while read id;do 
echo "gatk --java-options -Xmx20g GenotypeGVCFs  -R $ref -V ${id} -O ../COMVCF/${id%%.*}.vcf.gz" >>GenotypeGVCFs.sh;
done

```

### 合并cvf
```
 gatk --java-options "-Xmx20g " MergeVcfs -I raw_vcf.list -O camel.merge_raw.vcf.gz
```


## step 4 过滤VCF
```
#../softwore/gatk-4.2.6.0/gatk --java-options "-Xmx4g" VariantFiltration -R camel.fa -V  all.raw.snp.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 \
#|| SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'SNP_filter'  -O all.filter.snp.vcf


../softwore/gatk-4.2.6.0/gatk --java-options "-Xmx4g" VariantFiltration -R camel.fa -V  all.raw.indel.vcf --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 \
|| SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'indel_filter'  -O all.filter.indel.vcf
```
## step 5 群体结构
```
基因流分析
vcftools --gzvcf /data/dongwenjun/new_camel/quality/camel.vcf.gz --plink --out camel_non_missing
plink -real-ref-alleles --chr-set 36 --file camel_non_missing --indep-pairwise 50 5 0.5 --out camel_non_missing_filterLD
plink -real-ref-alleles --chr-set 36 --file camel_non_missing --extract camel_non_missing_filterLD.prune.in --make-bed --out camel_non_missing_filterLD
plink -real-ref-alleles --chr-set 36 --bfile camel_non_missing_filterLD --recode --out camel_non_missing_filterLD

plink -real-ref-alleles --chr-set 36 --bfile camel_non_missing_filterLD --freq --missing --within pop.txt --out camel_non_missing_filterLD
gzip camel_non_missing_filterLD.frq.strat
python2 /data/dongwenjun/software/plink2treemix.py camel_non_missing_filterLD.frq.strat.gz camel_non_missing_filterLD.treemix.frq.gz
gunzip camel_non_missing_filterLD.frq.strat
gunzip camel_non_missing_filterLD.treemix.frq.gz
awk 'BEGIN{print "scaffold_pos\tscaffold\tpos"}{split($2,pos,":");print $2"\t"pos[1]"\t"pos[2]}' camel_non_missing_filterLD.map > camel_non_missing_filterLD.positions
paste camel_non_missing_filterLD.positions camel_non_missing_filterLD.treemix.frq > camel_non_missing_filterLD.frequencies
awk '{printf $0; for(i = 4; i <= NF; i++){
                split($i,values,",")
                if((values[1]+values[2])>0) freq=values[1]/(values[1]+values[2])
                else freq=0
                printf freq"\t"};printf "\n"}' camel_non_missing_filterLD.frequencies > camel_non_missing_filterLD.frequencies2
mv camel_non_missing_filterLD.frequencies2 camel_non_missing_filterLD.frequencies
awk 'BEGIN{scaffold="";pos=0;newpos=0};{if($2==scaffold){newpos=pos+$3}else{scaffold=$2;pos=newpos};chpos=pos+$3;print $0,chpos}' camel_non_missing_filterLD.frequencies > camel_non_missing_filterLD.frequencies.newpos
gzip camel_non_missing_filterLD.treemix.frq



for i in {0..10}
do
 for j in {1..5}
 do
  treemix -i camel_non_missing_filterLD.treemix.frq.gz -m $i -o m10_rep5_500_outgroup/migration_${i}_dup${j} -k 500 -root Wild
 done
done

for i in {0..10}
do
 treemix -i camel_non_missing_filterLD.treemix.frq.gz -m $i -o m4/migration_${i} -k 500 -root Wild
done

for (i in 0:10){
+ pdf(paste('migration_',i,'_tree','.pdf', sep = ""), width = 8, height = 10)
+ plot_tree(paste('migration_',i,sep=""))
+ dev.off()
+ }

for (i in 0:10){
+ pdf(paste('migration_',i,'_resid','.pdf', sep = ""), width = 8, height = 10)
+ plot_resid(paste('migration_',i,sep=""), 'breed')
+ dev.off()
+ }

祖先成分分析
#!/bin/bash
cat K |while read id
do
 arr=(${id})
 K=${arr[0]}
 admixture --cv -j10 camel_filterLD.bed $K | tee CV/log${K}.out
done
```
## step6 选择信号的
```
方法1XPEHH：
1.拆分染色体
bcftools view --threads 20 -r 1 domestic.vcf.gz -Oz -o chr/1.vcf.gz
2.定相（拆分染色体会比较快）
java -Xmx2000m -jar /data/dongwenjun/miniconda3/share/beagle-5.2_21Apr21.304-0/beagle.jar\
 gt=/data/dongwenjun/camel/couple_pop_selection/domestic/chr/1.vcf.gz \
out=/data/dongwenjun/camel/couple_pop_selection/domestic/phased/1.vcf.gz ne=100 nthreads=2
两个群体都需做定相
3.计算群体间xpehh值（循环计算合并一般耗时半天到一天，看文件大小，单个染色体计算再合并比较快）
library(vcfR)
library(rehh)
library(R.utils)
n <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,38)
for(i in n) {
hh1 <- data2haplohh(hap_file = "Sunite_phased.vcf.gz",chr.name = i,polarize_vcf = FALSE,vcf_reader = "vcfR")
hh2 <- data2haplohh(hap_file = "domestic_phased.vcf.gz",chr.name = i,polarize_vcf = FALSE,vcf_reader = "vcfR")
scan1 <- scan_hh(hh1,threads=4)
scan2 <- scan_hh(hh2,threads=4)
if (i == 1) { wgscan1 <- scan1 } else { wgscan1 <- rbind(wgscan1, scan1) }
if (i == 1) { wgscan2 <- scan2 } else { wgscan2 <- rbind(wgscan2, scan2) }
}
xpehh.cgu_eut <- ies2xpehh(scan_pop1 = wgscan1,scan_pop2 = wgscan2,popname1 = "Sunite",popname2 = "domestic",standardize=FALSE)
#按窗口计算
cr.cgu <- calc_candidate_regions(xpehh.cgu_eut,threshold = -10000,window_size = 50000,overlap = 25000,min_n_mrk = 10,
min_n_extr_mrk = 1, min_perc_extr_mrk = 0,join_neighbors = FALSE)
write.table(cr.cgu,file="/data/dongwenjun/camel/couple_pop_selection/xpehh/s_d/Sunite_domestic_window_result.txt",
sep="\t",row.names=F,col.names=T,quote=F)


方法2 iHS：
1.2步同xpehh
3.计算iHS
library(vcfR)
library(rehh)
library(R.utils)
n <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,38)
for(i in n) {
hh1 <- data2haplohh(hap_file = "Sunite_phased.vcf.gz",chr.name = i,polarize_vcf = FALSE,vcf_reader = "vcfR")
scan1 <- scan_hh(hh1,threads=4)
if (i == 1) { wgscan1 <- scan1 } else { wgscan1 <- rbind(wgscan1, scan1) }
}
ihs <- ihh2ihs(wgscan1)
#按窗口计算
result_windows <- calc_candidate_regions(ihs,threshold = -1000,window_size = 50000,overlap = 25000,min_n_mrk = 10,
min_n_extr_mrk = 1, min_perc_extr_mrk = 0,join_neighbors =FALSE,pval = TRUE)
write.table(result_windows,file="/data/dongwenjun/camel/single_pop_selection/Sunite_ihs_window.txt",sep="\t",
row.names=F,col.names=T,quote=F)


方法3 Fst：
1.两两群体间计算
vcftools --vcf camel_control_remove11.vcf --weir-fst-pop Sunite.txt --weir-fst-pop Wild.txt --out\
 Sunite_Wild --fst-window-size 50000 --fst-window-step 25000
2.负值为计算有误，变为0
cat Sunite_Wild.windowed.weir.fst | awk '{if ($6>0){print $6}else{print 0}}' > tt
paste Sunite_Wild.windowed.weir.fst tt |cut -f 1-5,7 > Sunite_Wild1.windowed.weir.fst_none-.txt


方法4 PI：
1.单个群体计算
vcftools --vcf /data/dongwenjun/camel/quality/control_breed_vcf/Sunite.vcf --window-pi 50000\
 --window-pi-step 25000 --out Sunite
vcftools --vcf /data/dongwenjun/camel/quality/control_breed_vcf/Wild.vcf --window-pi 50000\
 --window-pi-step 25000 --out Wild
2.取-log10(s/w)转为两两群体间比较的形式
先确保两个文件窗口一致，对PI值取比值
#取出两个文件PI值所在列
paste S1 D1 | cut -f 5,10 > s_w_pi_cal
#计算比值和log值
data <- read.table(file="s_w_pi_cal",header=F)
result <- -log10(data$V1/data$V2)
result1 <- as.data.frame(result)
write.table(result1,file="s_w_log10_pi_ratio.txt",sep="\t",row.names=F,col.names=F,quote=F)
#与含窗口区域的文件合并
paste S1 s_w_log10_pi_ratio.txt | cut -f 1-3,6 > s_w_log10_pi_ratio.txt

```
