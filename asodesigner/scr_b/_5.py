import sys
import os

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# 使用上下文管理器确保文件正确关闭
with open(f"{uid}/Gresult.txt", "r") as f1, \
     open(f"{uid}/out4.txt", "r") as f2, \
     open(f"{uid}/out5.txt", "w") as fo:

    # 优化1: 更高效的G值解析
    dic_g = {}
    lines = f1.readlines()
    for i in range(0, len(lines), 3):
        if i + 2 >= len(lines):
            break  # 确保有完整的3行
            
        # 提取ASO名称和G值
        aso_name = lines[i].strip()[1:]  # 直接去掉>号
        g_line = lines[i+2].strip()
        
        # 优化2: 更健壮的G值提取
        g_value = g_line.rsplit(')', 1)[0]  # 从右边开始分割更安全
        g_value = g_value.rsplit('(', 1)[-1].strip()
        
        # 优化3: 直接尝试转换为浮点数
        try:
            dic_g[aso_name] = float(g_value)
        except (ValueError, TypeError):
            continue

    # 优化4: 批量处理输出文件
    for line in f2:
        parts = line.strip().split("\t")
        if not parts:
            continue
            
        aso_name = parts[0]
        g_val = dic_g.get(aso_name)
        
        # 优化5: 使用条件表达式简化
        if g_val is not None:
            fo.write(f"{line.strip()}\t{g_val}\n")
        else:
            # 保留原始行即使缺少G值
            fo.write(f"{line.strip()}\tNA\n")