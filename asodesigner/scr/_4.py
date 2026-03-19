import sys, os

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 优化反向互补计算
COMP_TABLE = str.maketrans('ATGCatgc', 'TACGtacg')
def DNA_complement1(sequence):
    return sequence.translate(COMP_TABLE)[::-1]

# 优化资源管理，使用上下文管理器
with (open(f"{uid}/out3.txt", "r") as f1,
      open("/asodesigner/text/human_rnaseH_pwm.txt", "r") as f2,
      open(f"{uid}/out4.txt", "w") as fo,
      open(f"{uid}/forG.txt", "w") as fo2,
      open(f"{uid}/foroff.fa", "w") as fo3):
    
    # 预定义禁止序列
    imis = frozenset(["TCCA", "GCTC", "CCTG", "GGG"])
    len1 = 16
    
    # 预加载PWM文件
    dic_pwm = {}
    for line in f2:
        parts = line.strip().split("\t")
        if len(parts) > 1:
            # 预处理为浮点数
            dic_pwm[parts[0]] = tuple(float(x) for x in parts[1:14])
    
    # 预定义PWM位置权重
    PWM_LENGTH = 13
    DEFAULT_SCORES = [0.0] * PWM_LENGTH
    
    # 优化字典加载
    dic_B = {}
    for line in f1:
        parts = line.strip().split("\t")
        if len(parts) > 4:  # 确保有足够的字段
            dic_B[parts[0]] = parts
    
    n = 0
    # 预编译序列处理方法
    for k, v in dic_B.items():
        sequence = v[-1]
        asoSEQ = sequence.upper()  # 提前转换为大写避免重复操作
        
        # 计算偏好性得分
        total_pf = 0.0
        # 4个偏移位置的滑动窗口
        for offset in range(4):
            start_idx = 2 + offset
            end_idx = start_idx + PWM_LENGTH
            sub_seq = sequence[start_idx:end_idx]
            
            # 获取互补序列
            sw = DNA_complement1(sub_seq).upper()
            
            # 计算当前窗口的得分
            pf = 0.0
            for j in range(min(len(sw), PWM_LENGTH)):
                scores = dic_pwm.get(sw[j], DEFAULT_SCORES)
                if j < len(scores):
                    pf += scores[j]
            total_pf += pf
        
        # 检查禁止序列
        if any(banned in asoSEQ for banned in imis):
            continue
        
        # 计算GC含量
        g_count = asoSEQ.count('G')
        c_count = asoSEQ.count('C')
        asoGC = (g_count + c_count) / len(asoSEQ)
        
        # 更新计数并准备输出
        n += 1
        asoname = f"ASO{n}"
        # 避免直接修改原始列表，使用新列表
        output_line = v + [str(total_pf), str(asoGC)]
        fo.write("\t".join(output_line) + "\n")
        
        # 写入其他文件
        fo2.write(f">{k}\n{asoSEQ}\n")
        fo3.write(f">{k}\n{DNA_complement1(sequence)}\n")