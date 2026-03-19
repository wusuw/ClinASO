import sys, os

nu_seq = int(sys.argv[1])
geneDirection = str(sys.argv[2])
chrn = sys.argv[3]
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 使用更高效的反向互补序列计算方法
COMP_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')

def DNA_complement1(sequence):
    """优化的反向互补序列计算"""
    return sequence.translate(COMP_TABLE)[::-1]

# 使用更高效的文件读取方式
with open(f"{uid}/out1.txt", "r") as f1, open(f"{uid}/out3.txt", "w") as fo:
    dic_seq = {}
    current_key = None
    
    # 逐行处理文件，避免加载大文件到内存
    for line in f1:
        line = line.strip()
        if line.startswith(">"):
            current_key = line[1:]
            dic_seq[current_key] = []
        elif current_key is not None:
            dic_seq[current_key].append(line)
    
    # 连接多行序列片段
    for key in dic_seq:
        dic_seq[key] = ''.join(dic_seq[key])
    
    a = 0
    # 优化循环和坐标计算
    for k, v in dic_seq.items():
        # 预处理序列（避免在循环中重复转换）
        if geneDirection == "-":
            seq = v.upper()
        else:
            seq = DNA_complement1(v.upper())
        
        # 预计算坐标范围
        start = int(k.split("-")[0])
        seq_len = len(seq)
        # 预生成坐标字符串模板
        coord_template = f"\t{chrn}\t{{}}\t{{}}\t"
        
        # 优化窗口序列生成
        for i in range(seq_len - nu_seq + 1):
            end_pos = start + i + nu_seq
            # 直接切片获取序列片段
            segment = seq[i:i+nu_seq]
            fo.write(f"ASO{a}{coord_template}{start+i}\t{end_pos}\t{segment}\n")
            a += 1