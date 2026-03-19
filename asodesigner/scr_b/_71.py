import sys
import os

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 优化1: 使用高效的反向互补计算（虽然未使用但保留）
COMP_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')
def DNA_complement1(sequence):
    """优化的反向互补序列计算"""
    return sequence.translate(COMP_TABLE)[::-1]

# 优化2: 直接使用集合推导式减少内存使用
with open(f"{uid}/gene2.gtf", "r") as f1, \
     open(f"{uid}/exon.txt", "w") as fo:
    
    # 使用集合推导式直接创建唯一项
    unique_entries = set()
    
    for line in f1:
        parts = line.strip().split("\t")
        if len(parts) < 5 or parts[2] != "CDS":
            continue
            
        # 优化3: 避免多次索引访问和类型转换
        chrom = parts[0]
        try:
            start = int(parts[3]) - 16
            end = int(parts[4]) + 16
        except ValueError:
            continue  # 跳过无效行
            
        # 优化4: 直接添加格式化的字符串
        unique_entries.add(f"{chrom}\t{start}\t{end}")
    
    # 优化5: 批量写入减少IO操作
    fo.write("\n".join(unique_entries))
    if unique_entries:  # 确保最后有换行
        fo.write("\n")