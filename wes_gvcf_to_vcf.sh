#!/usr/bin/bash

## 这是多样本变异检测流程的后半部分，这个部分假设你已经用fastq_to_gvcf.sh为每个样本独立生成了其对应的比对结果和gvcf

#nohup sh wes_gvcf_to_vcf.sh PID21-056A,PID21-057,PID21-058 /data/home/test2/LiDingYang/WES_pedigree PID21-056-057-058 &
#nohup sh wes_gvcf_to_vcf.sh PID21-232,PID21-233,PID21-234 /data/home/test2/LiDingYang/WES_pedigree PID21-232-233-234 &
#nohup sh wes_gvcf_to_vcf.sh A140135A,A140136A,A140417A,A140418A /data/home/test2/LiDingYang/WES_pedigree A140 &

#nohup sh wes_gvcf_to_vcf.sh A140135A,A140136A,A140417A,A140418A,PID21-056A,PID21-057,PID21-058,PID21-232,PID21-233,PID21-234 /data/home/test2/presentation/WES_pedigree result &


# GATK的路径, 根据实际路径
gatk=gatk

#reference
reference=/data/home/test2/Data_practice/NGS/resource
GATK_bundle=/data/home/test2/Data_practice/NGS/resource

## shell执行参数
samples=$1  ## 所有的样本ID，用","分开
indir=$2  ## 输入目录的路径，这个输入路径要与fastq_to_gvcf.sh的输出路径完全相同
outname=$3 ## 设置输出文件名的前缀

outdir=$indir ## 输入和输出路径相同

# 设置群体变异检测结果的输出目录
if [ ! -d $outdir/population ]
then mkdir -p $outdir/population
fi

## 按照","，把所有样本ID拆分出来存至数组中
samples=$(echo $samples | tr "," "\n")

## 基于群体数据进行Joint genotyping, 同样有两种方式
## 第一，先合并所有的全gvcf结果，然后再统一进行GenotypeGVCFs，由于是全基因组，速度较慢
sample_gvcfs=""
for sample in $samples ; do 
	sample_gvcfs=${sample_gvcfs}"--variant $outdir/${sample}/gatk/${sample}.HC.g.vcf.gz "
done

time $gatk CombineGVCFs \
	-R $reference/Homo_sapiens_assembly38.fasta \
	${sample_gvcfs}\
  -O $outdir/population/${outname}.HC.g.vcf.gz && echo "** ${outname}.HC.g.vcf.gz done ** " && \
time $gatk GenotypeGVCFs \
	-R $reference/Homo_sapiens_assembly38.fasta \
	-V $outdir/population/${outname}.HC.g.vcf.gz \
	-O $outdir/population/${outname}.HC.vcf.gz && echo "** ${outname}.HC.vcf.gz done ** "


# ## 第二，按照染色体分开，合并每个样本对应染色体的gVCF，然后再各自对染色体进行Joint calling(注意，这里和GATK3不同，GATK4只能接受一个gvcf的输入，因此需要先合并),
# ## 最后合并各个染色体的Genotype结果, 这要求每个样本必须按照染色体输出gvcf，速度快
# chrom=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM )
# for i in ${chrom[@]}; do

#     sample_gvcfs=""
#     for sample in $samples ; do
# 	    sample_gvcfs=${sample_gvcfs}"-V $outdir/${sample}/gatk/${sample}.HC.${i}.g.vcf.gz \\"\n
#     done
#     time $gatk CombineGVCFs \
#         -R $reference/Homo_sapiens_assembly38.fasta \
#         ${sample_gvcfs} \
#         -O $outdir/population/${outname}.HC.${i}.g.vcf.gz && echo "** ${outname}.HC.${i}.g.vcf.gz done ** " && \
#     time $gatk GenotypeGVCFs \
#         -R $reference/Homo_sapiens_assembly38.fasta \
#         -V $outdir/population/${outname}.HC.${i}.g.vcf.gz \
#         -O $outdir/population/${outname}.HC.${i}.vcf.gz && echo "** ${outname}.HC.${i}.vcf.gz done ** " &
# done && wait
# merge_vcfs=""
# for i in ${chrom[@]}; do
#     merge_vcfs=${merge_vcfs}"-I $outdir/population/${outname}.HC.${i}.vcf.gz \\"\n
# done && time $gatk MergeVcfs ${merge_vcfs} -O $outdir/population/${outname}.HC.vcf.gz && echo "** MergeVcfs done **"


## VQSR, 由于评价SNP和Indel质量高低的标准是不同的，因此，需要分SNP和Indel这两种不同的模式，分别进行质控

# ## 这一步不是必须的，仅仅是为了提高速度，如果你的样本数目不多比如（500以下），那么可以选择忽略.如果你的样本很多（>500）,
# ## 考虑到处理速度，我们可以用“MakeSitesOnlyVcf”只把sites的信息提取出来，忽略VCF后面的个体信息，从而节省IO时间。
# ## 如果提取了Only sites，那么注意下面VQSR中的输入文件要替换为这里的${outname}.HC.sites.vcf.gz
# time $gatk VariantRecalibrator \
#     -I $outdir/population/${outname}.HC.vcf.gz \
#     -O $outdir/population/${outname}.HC.sites.vcf.gz && echo "** Only sites done **"


#使用SelectVariants，选出SNP
#time $gatk SelectVariants --select-type-to-include SNP -V $outdir/population/${outname}.HC.vcf.gz -O $outdir/population/${outname}.HC.snp.vcf.gz
# 为SNP作过滤
#time $gatk VariantFiltration -V $outdir/population/${outname}.HC.snp.vcf.gz --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "PASS" -O $outdir/population/${outname}.snp.filter.vcf.gz
## 首先是SNP mode
time $gatk VariantRecalibrator \
              -R $reference/Homo_sapiens_assembly38.fasta \
              -V $outdir/population/${outname}.HC.vcf.gz \
              --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf.gz \
              --resource:omin,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf.gz \
              --resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
              --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
              -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
              -mode SNP \
              -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
              --rscript-file $outdir/population/${outname}.HC.snps.plots.R \
              --tranches-file $outdir/population/${outname}.HC.snps.tranches \
              -O $outdir/population/${outname}.HC.snps.recal && \
time $gatk ApplyVQSR \
              -R $reference/Homo_sapiens_assembly38.fasta \
              -V $outdir/population/${outname}.HC.vcf.gz \
              --truth-sensitivity-filter-level 99.0 \
              --tranches-file $outdir/population/${outname}.HC.snps.tranches \
              --recal-file $outdir/population/${outname}.HC.snps.recal \
              --mode SNP \
              -O $outdir/population/${outname}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"
## 然后是Indel mode
time $gatk VariantRecalibrator \
              -R $reference/Homo_sapiens_assembly38.fasta \
              -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
              --resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
              -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
              --mode INDEL \
              --max-gaussians 6 \
              --rscript-file $outdir/population/${outname}.HC.snps.indels.plots.R \
              --tranches-file $outdir/population/${outname}.HC.snps.indels.tranches \
              -O $outdir/population/${outname}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
              -R $reference/Homo_sapiens_assembly38.fasta \
              -V $outdir/population/${outname}.HC.snps.VQSR.vcf.gz \
              --truth-sensitivity-filter-level 95.0 \
              --tranches-file $outdir/population/${outname}.HC.snps.indels.tranches \
              --recal-file $outdir/population/${outname}.HC.snps.indels.recal \
              -mode INDEL \
              -O $outdir/population/${outname}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${outname}.HC.VQSR.vcf.gz finish) done **"

# 可以被删除清理的文件，这不是必须执行的,可以保留.g.vcf.gz，原始HC.vcf.gz和完成质控的HC.VQSR.vcf.gz
# rm -f $outdir/population/${outname}.HC.snps.VQSR.vcf.gz 
# $outdir/population/${outname}.HC.*.recal

#使用SelectVariants，选出Indel
#time $gatk SelectVariants --select-type-to-include INDEL -V $outdir/population/${outname}.HC.vcf.gz -O $outdir/population/${outname}.HC.indel.vcf.gz
# 为Indel作过滤
#time $gatk VariantFiltration -V $outdir/population/${outname}.HC.indel.vcf.gz --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filter-name "PASS" -O $outdir/population/${outname}.indel.filter.vcf.gz








