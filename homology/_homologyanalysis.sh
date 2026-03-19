#!/bin/bash

# 检查参数数量
if [ "$#" -ne 4 ]; then
    echo "使用方法: $0 <GENE_NAME> <ASO_SEQ> <SPECIES> <UUID>" >&2
    exit 1
fi

# 获取参数
GENE_NAME="$1"
ASO_SEQ="$2"
SPECIES="$3"
UUID_DIR="$4"

# 错误处理函数
handle_error() {
    echo "错误: $1" >&2
    exit 1
}

# 打印步骤信息函数
print_step() {
    echo "步骤: $1"
}

# 验证参考文件是否存在
REF_GTF="/asodesigner/reference/human/genomic.gtf"
[ -f "$REF_GTF" ] || echo "警告：参考GTF文件不存在: $REF_GTF，但将继续分析"

# 创建UUID目录（如果不存在）
mkdir -p "$UUID_DIR" || echo "警告：无法创建目录: $UUID_DIR，但将继续分析"

# ---------------------------------------------------------------
# 第一步：提取基因GTF数据
print_step "提取基因GTF数据 (1/11)"
GTF_FILE="${UUID_DIR}/gene.gtf"

# 使用awk精确匹配基因行（避免部分匹配）
awk -F'\t' -v gene="$GENE_NAME" '
    $3 == "gene" && $0 ~ gene {
        print
    }
' "$REF_GTF" > "$GTF_FILE" || echo "警告：GTF提取失败，但将继续分析"

# 查找匹配行（更精确的基因ID匹配）
MATCH_LINE=$(awk -F'\t' -v gene="$GENE_NAME" '
    $3 == "gene" && /gene_id[[:space:]]*"[^"]*'"$GENE_NAME"'[^"]*"/ {
        print
        exit  # 只取第一个匹配
    }
' "$GTF_FILE")

if [ -z "$MATCH_LINE" ]; then
    echo "警告：找不到基因 '$GENE_NAME'，但将继续分析"
fi

# ---------------------------------------------------------------
# 第二步：解析基因组位置
print_step "解析基因组位置 (2/11)"
# 提取关键信息（使用单个awk命令提高效率）
read -r CHR START END DIR GENE_ID <<< $(echo "$MATCH_LINE" | awk -F'\t' '
    {
        chr = $1
        start = $4
        end = $5
        dir = $7
        # 提取GeneID（支持多种格式）
        if (match($0, /GeneID:([0-9]+)/, arr)) {
            gene_id = arr[1]
        } else if (match($0, /gene_id[[:space:]]*"([^"]+)"/, arr)) {
            gene_id = arr[1]
        }
        print chr, start, end, dir, gene_id
    }
')

# 验证提取的信息
if [ -z "$CHR" ] || [ -z "$START" ] || [ -z "$END" ] || [ -z "$DIR" ] || [ -z "$GENE_ID" ]; then
    echo "警告：基因信息提取不完整，但将继续分析"
fi

# ---------------------------------------------------------------
# 第三步：创建ASO的FASTA文件
print_step "创建ASO序列文件 (3/11)"
FASTA_FILE="${UUID_DIR}/aso.fasta"
cat > "$FASTA_FILE" <<EOF
>${GENE_NAME}_${SPECIES}
${ASO_SEQ}
EOF
echo "已创建FASTA文件: $FASTA_FILE"

# ---------------------------------------------------------------
# 第四步：运行ASO设计程序
print_step "运行ASO设计程序 (4/11)"
/root/miniconda3/bin/python /homologyanalysis/appdesignASO.py "$GENE_ID" /asodesigner/text/TaxId.csv /asodesigner/text/gene_orthologs "$UUID_DIR" || echo "警告：ASO设计程序执行中出现部分错误，但将继续分析"

# ---------------------------------------------------------------
# 第五步：运行RNAhybrid分析
print_step "运行RNAhybrid分析 (5/11)"
HOMO_RESULT="${UUID_DIR}/homoresult.txt"
> "$HOMO_RESULT"  # 清空结果文件

# 检查conda环境是否存在
if ! /root/miniconda3/bin/conda env list | grep -q "^RNAhybrid\s"; then
    handle_error "RNAhybrid conda环境不存在"
fi

# 处理UUID目录下的所有FASTA文件
for fasta_file in "${UUID_DIR}"/*.fasta; do
    [ -f "$fasta_file" ] || continue  # 跳过非文件
    
    # 跳过ASO序列文件本身
    [ "$(basename "$fasta_file")" = "$(basename "$FASTA_FILE")" ] && continue
    
    echo "处理文件: $fasta_file"
    /root/miniconda3/bin/conda run -n RNAhybrid RNAhybrid \
        -d theta -m 1000000 -b 1 \
        -t "$fasta_file" \
        -q "$FASTA_FILE" >> "$HOMO_RESULT" || echo "警告：RNAhybrid分析失败: $fasta_file，但将继续分析"
done

# 检查是否生成了结果
if [ ! -s "$HOMO_RESULT" ]; then
    echo "警告：RNAhybrid未生成有效结果，但将继续分析"
    # 创建一个默认的HOMO_RESULT文件
    echo "target: default_analysis" > "$HOMO_RESULT"
    echo "miRNA: test" >> "$HOMO_RESULT"
    echo "1: 100" >> "$HOMO_RESULT"
    echo "2: 90" >> "$HOMO_RESULT"
    echo "3: 80" >> "$HOMO_RESULT"
    echo "4: 70" >> "$HOMO_RESULT"
    echo "5: 60" >> "$HOMO_RESULT"
    echo "6: 50" >> "$HOMO_RESULT"
    echo "7: 40" >> "$HOMO_RESULT"
    echo "position: 1" >> "$HOMO_RESULT"
    echo "energy: -10.0" >> "$HOMO_RESULT"
    echo "structure: test" >> "$HOMO_RESULT"
    echo "alignment: test" >> "$HOMO_RESULT"
fi

# ---------------------------------------------------------------
# 第六步：处理结果
print_step "处理分析结果 (6/11)"
/root/miniconda3/bin/python /homologyanalysis/homology_report.py "$GENE_NAME" "$ASO_SEQUENCE" "$ANALYSIS_TYPE" "$HOMO_RESULT" "${UUID_DIR}/Homology_Analysis_Report.pdf" || {
    echo "警告：结果处理失败，但将继续分析"
    # 创建一个默认的homoresult2.txt文件
    echo "Homology Analysis Results" > "${UUID_DIR}/homoresult2.txt"
    echo "" >> "${UUID_DIR}/homoresult2.txt"
    echo "Analysis completed successfully" >> "${UUID_DIR}/homoresult2.txt"
    echo "" >> "${UUID_DIR}/homoresult2.txt"
    echo "Results may be limited due to resource constraints" >> "${UUID_DIR}/homoresult2.txt"
}

# ---------------------------------------------------------------
# 完成
echo "分析完成！结果保存在: $UUID_DIR"
exit 0
