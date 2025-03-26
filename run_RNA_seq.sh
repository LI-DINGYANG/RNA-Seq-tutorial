#!/bin/bash
#SBATCH -J rna
#SBATCH -o rna_%j.out
#SBATCH -e rna_%j.err
#SBATCH --nodes=1
#SBATCH --partition=warshel-cpu
#SBATCH --ntasks-per-node=4
# RNA-seq 自动化分析流程脚本
# 用于批量处理多个样本，使用共享数据资源，结果输出到用户目录
# 请设置你的学号
USER_ID="222051031"  # 更改成你的学号
# 设置目录
DATA_DIR="/home/222051031/RNASEQ"  # 共享数据资源目录,222051031不要改动
USER_DIR="/home/$USER_ID/RNASEQ"   # 用户结果目录
# 设置参数
THREADS=4
MEMORY="8G"
# 定义样本数组
SAMPLES=(
    "SRR26314914"
    "SRR26314915"
    "SRR26314919"
    "SRR26314921"
    "SRR26314922"
    "SRR26314923"
)

echo "=== 设置工作环境 ==="
echo "数据资源目录: $DATA_DIR"
echo "结果输出目录: $USER_DIR"

# 创建用户结果目录结构
echo -e "\n=== 创建用户目录结构 ==="
mkdir -p ${USER_DIR}/{results,trimedfq,report,qualimap}
echo "✓ 用户目录结构创建完成"

module load anaconda3/2024.06
conda activate RNA-seq

# 遍历处理每个样本
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\n\n=========================================="
    echo "🧬 开始处理样本: ${SAMPLE}"
    echo "==========================================\n"

    # 1. 原始数据质控
    echo -e "\n=== 步骤1: 原始数据质控 - ${SAMPLE} ==="
    if [ -f ${DATA_DIR}/rawfq/${SAMPLE}_1.fastq.gz ] && [ -f ${DATA_DIR}/rawfq/${SAMPLE}_2.fastq.gz ]; then
        fastqc ${DATA_DIR}/rawfq/${SAMPLE}_*.fastq.gz -o ${USER_DIR}/report/beforeQC -t ${THREADS}
        echo "✓ FastQC 完成 - ${SAMPLE}"
    else
        echo "⚠️ 警告: 找不到原始数据文件 ${DATA_DIR}/rawfq/${SAMPLE}_1.fastq.gz 和/或 ${DATA_DIR}/rawfq/${SAMPLE}_2.fastq.gz"
        echo "尝试查找替代命名格式..."
        
        # 尝试替代命名格式
        if [ -f ${DATA_DIR}/rawfq/${SAMPLE}.fastq.gz ] || [ -f ${DATA_DIR}/rawfq/${SAMPLE}_R1.fastq.gz ]; then
            FOUND_FILES=$(ls ${DATA_DIR}/rawfq/${SAMPLE}*.fastq.gz 2>/dev/null)
            echo "找到以下文件: ${FOUND_FILES}"
            fastqc ${DATA_DIR}/rawfq/${SAMPLE}*.fastq.gz -o ${USER_DIR}/report/beforeQC -t ${THREADS}
            echo "✓ FastQC 完成 (使用替代文件) - ${SAMPLE}"
        else
            echo "❌ 错误: 无法找到任何与 ${SAMPLE} 相关的原始数据文件，跳过此样本"
            continue
        fi
    fi

    # 2. 数据质量过滤
    echo -e "\n=== 步骤2: Fastp 质量过滤 - ${SAMPLE} ==="
    if [ -f ${DATA_DIR}/rawfq/${SAMPLE}_1.fastq.gz ] && [ -f ${DATA_DIR}/rawfq/${SAMPLE}_2.fastq.gz ]; then
        fastp -i ${DATA_DIR}/rawfq/${SAMPLE}_1.fastq.gz -I ${DATA_DIR}/rawfq/${SAMPLE}_2.fastq.gz \
              -o ${USER_DIR}/trimedfq/${SAMPLE}_1.trimmed.fastq.gz -O ${USER_DIR}/trimedfq/${SAMPLE}_2.trimmed.fastq.gz \
              -h ${USER_DIR}/report/${SAMPLE}_fastp.html -j ${USER_DIR}/report/${SAMPLE}_fastp.json
        echo "✓ Fastp 完成 - ${SAMPLE}"
    else
        echo "❌ 错误: 无法执行Fastp，找不到必要的输入文件，跳过此样本的后续步骤"
        continue
    fi

    # 3. 过滤后数据质控
    echo -e "\n=== 步骤3: 过滤后数据质控 - ${SAMPLE} ==="
    fastqc ${USER_DIR}/trimedfq/${SAMPLE}_*.trimmed.fastq.gz -o ${USER_DIR}/report/afterQC -t ${THREADS}
    echo "✓ 过滤后 FastQC 完成 - ${SAMPLE}"

    # 所有样本处理完成后，运行整体MultiQC报告
    echo -e "\n=== 生成所有样本的MultiQC报告 ==="
    multiqc ${USER_DIR}/report/beforeQC -o ${USER_DIR}/report/beforeQC/multiqc_report
    multiqc ${USER_DIR}/report/afterQC -o ${USER_DIR}/report/afterQC/multiqc_report
    echo "✓ MultiQC 完成"

    # 4. STAR比对
    echo -e "\n=== 步骤4: STAR比对 - ${SAMPLE} ==="
    STAR --runThreadN ${THREADS} \
         --genomeDir ${DATA_DIR}/genomedir \
         --readFilesIn ${USER_DIR}/trimedfq/${SAMPLE}_1.trimmed.fastq.gz ${USER_DIR}/trimedfq/${SAMPLE}_2.trimmed.fastq.gz \
         --outFileNamePrefix ${USER_DIR}/results/${SAMPLE} \
         --readFilesCommand zcat \
         --outSAMtype BAM Unsorted \
         --quantTranscriptomeBan Singleend \
         --outFilterType BySJout \
         --alignSJoverhangMin 8 \
         --outFilterMultimapNmax 20 \
         --alignSJDBoverhangMin 1 \
         --outFilterMismatchNmax 999 \
         --outFilterMismatchNoverReadLmax 0.04 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --alignMatesGapMax 1000000 \
         --quantMode TranscriptomeSAM \
         --outSAMattributes NH HI AS NM MD
    echo "✓ STAR 比对完成 - ${SAMPLE}"

    # 5. SAM处理
    echo -e "\n=== 步骤5: SAM处理 - ${SAMPLE} ==="
    samtools sort -@${THREADS} ${USER_DIR}/results/${SAMPLE}Aligned.out.bam > ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.bam
    samtools index -@${THREADS} ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.bam
    echo "✓ SAM排序和索引完成 - ${SAMPLE}"

    # 6. 比对统计
    echo -e "\n=== 步骤6: 比对统计 - ${SAMPLE} ==="
    samtools flagstat -@${THREADS} ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.bam > ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.flagstat
    echo "✓ 比对统计完成 - ${SAMPLE}"

    # 7. Qualimap质控
    echo -e "\n=== 步骤7: Qualimap BAM质控 - ${SAMPLE} ==="
    qualimap bamqc -bam ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.bam \
           -gff ${DATA_DIR}/genomedir/gencode.v47.annotation.gtf \
           -outdir ${USER_DIR}/qualimap/${SAMPLE}-bamqc-qualimap-report \
           --java-mem-size=${MEMORY}
    echo "✓ Qualimap BAM质控完成 - ${SAMPLE}"

    # 8. Qualimap RNA-seq分析
    echo -e "\n=== 步骤8: Qualimap RNA-seq分析 - ${SAMPLE} ==="
    qualimap rnaseq -bam ${USER_DIR}/results/${SAMPLE}Aligned.out.sorted.bam \
           -gtf ${DATA_DIR}/genomedir/gencode.v47.annotation.gtf \
           -outdir ${USER_DIR}/qualimap/${SAMPLE}-rnaseq-qualimap-report \
           --java-mem-size=${MEMORY}
    echo "✓ Qualimap RNA-seq分析完成 - ${SAMPLE}"

    multiqc ${USER_DIR}/qualimap -o ${USER_DIR}/qualimap/multiqc_report

    # 9. Salmon定量
    echo -e "\n=== 步骤9: Salmon转录本定量 - ${SAMPLE} ==="
    salmon quant -t ${DATA_DIR}/genomedir/GRCh38_no_alt_analysis_set_gencode.v47.transcripts.fa \
           --libType A \
           -a ${USER_DIR}/results/${SAMPLE}Aligned.toTranscriptome.out.bam \
           -o ${USER_DIR}/results/${SAMPLE}.salmon_quant \
           --gcBias --seqBias -p ${THREADS}
    echo "✓ Salmon定量完成 - ${SAMPLE}"

    echo -e "\n✅ 样本 ${SAMPLE} 分析完成!"
done


# 分析完成
echo -e "\n=========================================="
echo "🎉 全部6个样本的RNA-seq分析流程已完成!"
echo "结果文件位于: ${USER_DIR}/results"
echo "质控报告位于: ${USER_DIR}/report"
echo "综合MultiQC报告: ${USER_DIR}/report/multiqc_report"
echo "=========================================="
