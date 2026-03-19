import sys
import os
from collections import defaultdict, Counter

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 优化1: 使用defaultdict自动初始化列表
with open("/asodesigner/text/human_blstn_result.txt", "r") as f1, \
     open(f"{uid}/out5.txt", "r") as f2, \
     open(f"{uid}/out6.txt", "w") as fo:

    # 优化2: 使用更高效的数据结构
    dic_offgene = defaultdict(list)
    dic_genename = defaultdict(list)
    
    # 处理blast结果文件
    for line in f1:
        parts = line.strip().split("\t")
        if len(parts) < 4:
            continue
        
        seq_id = parts[0]
        gene = parts[3]
        offgene = parts[1].split("_")[0]  # 优化3: 避免重复split
        
        # 优化4: defaultdict自动处理键存在情况
        dic_offgene[seq_id].append(offgene)
        dic_genename[seq_id].append(gene)
    
    # 处理out5文件
    for line in f2:
        parts = line.strip().split("\t")
        if not parts:
            continue
            
        seq_id = parts[0]
        
        # 获取offgene信息
        offgenes = dic_offgene.get(seq_id, [])
        offagene = "/".join(set(offgenes)) if offgenes else "NA"  # 优化5: 去重
        
        # 获取genename计数
        genes = dic_genename.get(seq_id, [])
        gene_counter = Counter(genes)  # 优化6: 使用Counter高效计数
        
        # 优化7: 直接取值，使用默认值0
        counts = [
            gene_counter.get("16", 0),
            gene_counter.get("17", 0),
            gene_counter.get("18", 0),
            gene_counter.get("19", 0),
            gene_counter.get("16", 0)  # 注意这里是重复值
        ]
        
        # 优化8: 使用格式化字符串高效输出
        count_str = "\t".join(str(c) for c in counts)
        fo.write(f"{line.strip()}\t{count_str}\t{offagene}\n")