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
USER_DIR="/home/$USER_ID/Assignment2"   # 用户结果目录
# 设置参数
THREADS=4
# 定义样本数组
SAMPLES=(
    "case1"
    "case2"
    "case3"
    "case4"
    "control1"
    "control2"
    "control3"
    "control4"
)

echo "=== 设置工作环境 ==="
echo "结果输出目录: $USER_DIR"

# 创建用户结果目录结构
echo -e "\n=== 创建用户目录结构 ==="
mkdir -p ${USER_DIR}/results
echo "✓ 用户目录结构创建完成"

#module load anaconda3/2024.06
#conda activate RNA-seq

#cd ${USER_DIR}
#wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.transcripts.fa.gz && gunzip gencode.vM36.annotation.gtf.gz
#wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M36/gencode.vM36.annotation.gtf.gz && gunzip gencode.vM36.transcripts.fa.gz

# 遍历处理每个样本
for SAMPLE in "${SAMPLES[@]}"; do
    echo -e "\n\n=========================================="
    echo "🧬 开始处理样本: ${SAMPLE}"
    echo "==========================================\n"


    #Salmon index
    echo -e "\n=== 步骤9: Salmon转录本定量 - ${SAMPLE} ==="
    #salmon index \
    #-t gencode.vM36.transcripts.fa \
    #-i salmon_index \
    #-p 12 \
    #--gencode

    #Salmon定量
    salmon quant -i /home/222051031/Assignment2/salmon_index -l A \
    -1 /home/yfwang/BIM3001_AS2/${SAMPLE}_1.fq.gz \
    -2 /home/yfwang/BIM3001_AS2/${SAMPLE}_2.fq.gz  \
    -p $THREADS --output ${USER_DIR}/results/${SAMPLE}.salmon_quant \
    --gcBias \
    --validateMappings


    echo "✓ Salmon定量完成 - ${SAMPLE}"

    echo -e "\n✅ 样本 ${SAMPLE} 分析完成!"
done



MERGED_FILE="${USER_DIR}/merged_salmon_counts_matrix.txt"
echo "=== 开始合并counts数据 ==="

# 创建表头
echo -n "GeneID" > ${MERGED_FILE}
for SAMPLE in "${SAMPLES[@]}"; do
    echo -n -e "\t${SAMPLE}" >> ${MERGED_FILE}
done
echo "" >> ${MERGED_FILE}

# 使用awk一次性处理所有数据
gawk -v outfile="${MERGED_FILE}" -v user_dir="${USER_DIR}" '
BEGIN {
    # 从第一个参数获取样本列表
    samples_count = split(ARGV[1], samples, ",")
    delete ARGV[1]  # 删除样本参数，以免被处理为文件
    
    # 读取所有样本数据
    for (s = 1; s <= samples_count; s++) {
        sample = samples[s]
        file = user_dir "/results/" sample ".salmon_quant/quant.sf"
        
        # 跳过标题行，读取基因和counts
        getline < file  # 跳过标题行
        while ((getline line < file) > 0) {
            split(line, fields, "\t")
            gene_id = fields[1]
            counts = fields[5]
            data[gene_id, sample] = counts
            all_genes[gene_id] = 1  # 记录所有出现的基因
        }
        close(file)
        print "已处理: " sample
    }
    
    # 输出所有基因的counts值
    for (gene_id in all_genes) {
        printf("%s", gene_id) >> outfile
        
        for (s = 1; s <= samples_count; s++) {
            sample = samples[s]
            if ((gene_id, sample) in data) {
                printf("\t%s", data[gene_id, sample]) >> outfile
            } else {
                printf("\tNA") >> outfile
            }
        }
        printf("\n") >> outfile
    }
}' "$(IFS=,; echo "${SAMPLES[*]}")"

gawk '{
     if (match($0, /transcript_id "([^"]+)"/, id) && match($0, /gene_name "([^"]+)"/, name)) {
         print id[1] "\t" name[1];
     }
 }' /home/222051031/Assignment2/gencode.vM36.annotation.gtf | sort -u > ${USER_DIR}/gene_id_to_symbol.txt

gawk 'BEGIN {FS=OFS="\t"}
     NR==FNR {map[$1]=$2; next} 
     $1 in map {$1=map[$1]} 
     1' ${USER_DIR}/gene_id_to_symbol.txt ${USER_DIR}/merged_salmon_counts_matrix.txt | gawk -F'\t' -v OFS='\t' '
NR==1 {print; next}  # 保留并打印第一行，然后跳到下一行
{
    sum = 0;
    for (i=2; i<=9; i++) {
        sum += $i + 0  # 强制转换为数值
    }
    if (sum >= 10) print  # 只保留总和≥10的行
}' > ${USER_DIR}/output.txt && mv ${USER_DIR}/output.txt ${USER_DIR}/merged_salmon_counts_matrix.txt



# 验证结果
GENE_COUNT=$(($(wc -l < ${MERGED_FILE}) - 1))
echo "✓ 合并完成！矩阵包含 ${GENE_COUNT} 个基因和 ${#SAMPLES[@]} 个样本"
echo "✓ 结果文件: ${MERGED_FILE}"
echo "前几行示例:"
head -n 10 ${MERGED_FILE} | column -t -s $'\t'

# 分析完成
echo -e "\n=========================================="
echo "🎉 全部8个样本的RNA-seq分析流程已完成!"
echo "结果文件位于: ${USER_DIR}/merged_salmon_counts_matrix.txt"
echo "=========================================="

