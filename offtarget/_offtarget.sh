#!/bin/bash

# 检查参数数量
if [ "$#" -ne 3 ]; then
    echo "使用方法: $0 <GENE_NAME> <ASO_seq> <UUID_DIR>"
    exit 1
fi

# 获取参数
GENE_NAME="$1"
ASO_seq="$2"
UUID_DIR="$3"
# 错误处理函数
handle_error() {
    echo "错误: $1" >&2
    exit 1
}

# 打印步骤信息函数
print_step() {
    echo "步骤: $1"
}

# 创建必要的目录
mkdir -p "$UUID_DIR" || handle_error "无法创建目录 $UUID_DIR"

# 第一步：生成FASTA文件
print_step "生成FASTA文件 (1/1)"

# 检查ASO序列是否有效（只包含ATCG字符）
if [[ ! "$ASO_seq" =~ ^[ATCGatcg]+$ ]]; then
    handle_error "ASO序列包含无效字符: $ASO_seq"
fi

# 生成FASTA文件
FASTA_FILE="$UUID_DIR/aso.fasta"
{
    echo ">${GENE_NAME}_${ASO_seq}"
    echo "$ASO_seq"
} > "$FASTA_FILE"

# 检查文件是否创建成功
if [ ! -f "$FASTA_FILE" ]; then
    handle_error "无法创建FASTA文件: $FASTA_FILE"
fi

# 显示文件信息
echo "FASTA文件已创建: $FASTA_FILE"
echo "文件大小: $(stat -c%s "$FASTA_FILE" 2>/dev/null || echo "unknown") 字节"
echo "文件内容:"
cat "$FASTA_FILE"

# 第二步：运行RNAhybrid分析
print_step "运行RNAhybrid分析 (2/3)"

# 确保输出文件存在
OUTPUT_FILE="$UUID_DIR/output.txt"
touch "$OUTPUT_FILE"

# 运行RNAhybrid
/root/miniconda3/bin/conda run -n RNAhybrid RNAhybrid \
    -d theta -m 1000000 -b 1 -c \
    -t /offtarget/GRCh38.gene2.fasta \
    -q "$FASTA_FILE" >> "$OUTPUT_FILE" || echo "警告：RNAhybrid分析失败，但将继续分析"

# 检查RNAhybrid输出
if [ ! -s "$OUTPUT_FILE" ]; then
    echo "警告：RNAhybrid未生成有效输出: $OUTPUT_FILE，但将继续分析"
fi

echo "RNAhybrid分析完成，输出文件: $OUTPUT_FILE"

# 第三步：分析结果到文件
print_step "分析结果 (3/3)"
# 修复引号问题
/root/miniconda3/bin/python /offtarget/_plot.py "$UUID_DIR/output.txt" "$UUID_DIR/top15genes.pdf" || echo "警告：结果分析失败，但将继续分析"

echo "分析完成！结果保存在: $UUID_DIR"
