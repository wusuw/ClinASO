#!/bin/bash
gene_name=$1
user_seq=${2^^}
UUID_DIR=$3

gtf_file="/asodesigner/reference/human/genomic.gtf"
genome_fasta="/asodesigner/reference/human/genomic.fna"
snp_db_157="/snp/GCF_000001405.40.gz"

# 创建输出目录
mkdir -p "$UUID_DIR"

# Step 1: 提取基因坐标（使用您原来的逻辑）
awk -v gene="$gene_name" 'BEGIN{FS="\t"; OFS="\t"} $3=="gene" && $9 ~ "gene_id \"" gene "\"" {print $1,$4,$5,$7}' "$gtf_file" > "$UUID_DIR/gene_coordinates.txt"

# 检查是否找到基因
if [ ! -s "$UUID_DIR/gene_coordinates.txt" ]; then
    echo "错误: 未找到基因 $gene_name" >&2
    exit 1
fi

# Step 2: 提取基因组序列
bedtools getfasta -fi "$genome_fasta" -bed "$UUID_DIR/gene_coordinates.txt" -s -nameOnly | awk '{if(NR==1) print; else print toupper($0)}' > "$UUID_DIR/gene_sequence.fa"

# Step 3: 根据链方向处理用户序列
chrom=$(awk '{print $1}' "$UUID_DIR/gene_coordinates.txt")
strand=$(awk '{print $4}' "$UUID_DIR/gene_coordinates.txt")
gene_start=$(awk '{print $2}' "$UUID_DIR/gene_coordinates.txt")
gene_end=$(awk '{print $3}' "$UUID_DIR/gene_coordinates.txt")

if [ "$strand" == "+" ]; then
    # 正链：用户序列需要反向互补
    user_seq_rev=$(echo "$user_seq" | rev)
    user_seq_rc=$(echo "$user_seq_rev" | tr 'ATGCatgc' 'TACGtacg')
    echo ">user_sequence_rc" > "$UUID_DIR/user_seq_rc.fa"
    echo "$user_seq_rc" >> "$UUID_DIR/user_seq_rc.fa"
    search_seq="$user_seq_rc"
else
    # 负链：用户序列不需要反向互补
    echo ">user_sequence" > "$UUID_DIR/user_seq_rc.fa"
    echo "$user_seq" >> "$UUID_DIR/user_seq_rc.fa"
    search_seq="$user_seq"
fi

# Step 4: 比对和坐标转换
# 读取基因序列
gene_seq=$(tail -n +2 "$UUID_DIR/gene_sequence.fa" | tr -d '\n')

# 在基因序列中查找用户序列
start=$(echo "$gene_seq" | grep -aob "$search_seq" | head -1 | awk -F: '{print $1}')

if [ -z "$start" ]; then
    echo "错误: 在基因序列中未找到用户序列" >&2
    exit 1
fi

# grep返回的是0-based索引，转换为1-based
start=$((start + 1))
end=$((start + ${#search_seq} - 1))

# Step 5: 转换为基因组坐标
if [ "$strand" == "+" ]; then
    # 正链：从基因起始位置向后加
    genome_start=$((gene_start + start))
    genome_end=$((gene_start + end))
else
    # 负链：从基因末端位置向前减
    genome_start=$((gene_end - end))
    genome_end=$((gene_end - start))
fi

# 确保起始位置不大于结束位置
if [ "$genome_start" -gt "$genome_end" ]; then
    temp=$genome_start
    genome_start=$genome_end
    genome_end=$temp
fi

# 保存位置信息
echo "基因序列位置: $start-$end" > "$UUID_DIR/position_info.txt"
echo "基因组位置: $chrom:$genome_start-$genome_end" >> "$UUID_DIR/position_info.txt"
echo "链方向: $strand" >> "$UUID_DIR/position_info.txt"

# Step 6: 提取SNP
/root/miniconda3/bin/bcftools view -r "${chrom}:${genome_start}-${genome_end}" "$snp_db_157" > "$UUID_DIR/snps_output.vcf" 2>/dev/null
/root/miniconda3/bin/bcftools view -r "${chrom}:${genome_start}-${genome_end}" /snp/freq.vcf.gz > "$UUID_DIR/snpsfreq_output.vcf" 2>/dev/null

# 检查是否成功提取到SNP
if [ ! -s "$UUID_DIR/snps_output.vcf" ]; then
    echo "警告: 在指定区域未找到SNP信息" > "$UUID_DIR/snps_output.vcf"
fi

if [ ! -s "$UUID_DIR/snpsfreq_output.vcf" ]; then
    echo "警告: 在指定区域未找到SNP频率信息" > "$UUID_DIR/snpsfreq_output.vcf"
fi
echo $gene_name $user_seq $chrom $genome_start $genome_end $strand $UUID_DIR
python3 /snp/plot.py $gene_name $user_seq $chrom $genome_start $genome_end $strand $UUID_DIR

echo "处理完成q1"
echo "基因组位置: $chrom:$genome_start-$genome_end"
echo "输出文件保存在: $UUID_DIR" 
