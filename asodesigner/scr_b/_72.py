import sys
import os

fw = sys.argv[1]
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 优化1: 高效的反向互补计算
COMP_TABLE = str.maketrans('ATGCatgcN', 'TACGtacgN')
def DNA_complement1(sequence):
    return sequence.translate(COMP_TABLE)[::-1]

# 优化2: 使用上下文管理器确保文件正确关闭
with (open(f"{uid}/exonfasta.fa", "r") as f1,
      open(f"{uid}/out7.txt", "r") as f2,
      open(f"{uid}/result.txt", "w") as fo):
    
    # 优化3: 高效读取FASTA文件
    dic_seq = {}
    current_key = ""
    for line in f1:
        line = line.strip()
        if line.startswith(">"):
            current_key = line[1:]
            dic_seq[current_key] = []
        elif current_key:
            dic_seq[current_key].append(line)
    
    # 优化4: 预处理外显子序列为大写并连接
    exon_seqs = set()
    for key, seq_list in dic_seq.items():
        seq_str = ''.join(seq_list).upper()
        exon_seqs.add(seq_str)  # 仅存储唯一序列，减少重复检查
    
    # 优化5: 批量处理并减少冗余检查
    output_lines = []
    if fw == "+":
        for line in f2:
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            dna_seq = DNA_complement1(parts[4]).upper()
            found = any(dna_seq in exon_seq for exon_seq in exon_seqs)
            output_lines.append(f"{line.strip()}\t{'exon' if found else 'intron'}")
    else:  # fw == "-"
        for line in f2:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            dna_seq = parts[1].upper()
            found = any(dna_seq in exon_seq for exon_seq in exon_seqs)
            output_lines.append(f"{line.strip()}\t{'exon' if found else 'intron'}")
    
    # 优化6: 批量写入减少IO操作
    fo.write("\n".join(output_lines))