import numpy as np
import os, sys
import glob
from multiprocessing import Pool
import re
import pandas as pd

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
aso_length = int(sys.argv[-3])
# 碱基编码字典（优化内存占用）
base_to_code = {'A': 0, 'T': 1, 'C': 2, 'G': 3, 
                'a': 0, 't': 1, 'c': 2, 'g': 3}

def DNA_complement(sequence):
    """反向互补序列获取"""
    comp_dict = {"A": "T", "T": "A", "G": "C", "C": "G",
                 "a": "t", "t": "a", "g": "c", "c": "g"}
    return ''.join([comp_dict[base] for base in sequence[::-1]])

def process_sequence(current_seq, kmers, encoded):
    """处理单个序列的公共逻辑"""
    seq = ''.join(current_seq).upper()
    for i in range(len(seq)-aso_length+1):
        kmer = seq[i:i+aso_length]
        kmers.add(kmer)
        encoded.append([base_to_code.get(c, -1) for c in kmer])  # 使用get处理未知字符

def read_kmers_optimized(file):
    """读取FASTA文件并生成数值编码的k-mer矩阵"""
    kmers = set()
    encoded = []
    current_seq = []
    for line in file:
        if line.startswith(">"):
            if current_seq:
                process_sequence(current_seq, kmers, encoded)
                current_seq = []
        else:
            current_seq.append(line.strip())
    if current_seq:
        process_sequence(current_seq, kmers, encoded)
    return kmers, np.array(encoded, dtype=np.int8)

def get_max_match_vectorized(target, full_set, full_array):
    """向量化比对函数"""
    target_upper = target.upper()
    if target_upper in full_set:
        return 1.0
    try:
        target_enc = np.array([base_to_code.get(c, -1) for c in target_upper], dtype=np.int8)
        if len(target_enc) != aso_length:
            return 0.0
        matches = np.sum(full_array == target_enc, axis=1)
        return np.max(matches) / aso_length
    except (ValueError, TypeError):
        return 0.0

def auto_detect_species_files():
    """自动检测物种文件 - 在 /asodesigner/outfile/ 目录下查找"""
    # 指定搜索目录
    search_dir = f"{uid}"
    
    # 物种映射
    species_map = {
        'crab_eating_macaque': ['crab_eating_macaque'],
        'mouse': ['mouse'],
        'rat': ['rat'],
        'pig': ['pig'],
        'rabbit': ['rabbit'],
        "Guinea_pig": ["guinea_pig", "guineapig", "cavia"],
        "goose": ["goose", "anser"]
    }
    
    detected = {}
    
    # 确保目录存在
    if not os.path.exists(search_dir):
        print(f"警告: 目录 '{search_dir}' 不存在!")
        return detected
    
    # 获取目录下所有fasta文件
    fasta_files = glob.glob(os.path.join(search_dir, "*.*fa*"))
    
    # 匹配物种文件
    for species, keywords in species_map.items():
        matched_files = []
        
        # 遍历所有可能的文件
        for file_path in fasta_files:
            filename = os.path.basename(file_path).lower()
            
            # 检查关键字
            for keyword in keywords:
                if keyword.lower() in filename:
                    matched_files.append(file_path)
                    break
        
        # 如果找到匹配项，选择最新的文件
        if matched_files:
            # 按修改时间排序
            matched_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            detected[species] = matched_files[0]
    
    # 打印检测结果
    print("\n检测到的物种文件:")
    for species, path in detected.items():
        print(f"  - {species}: {os.path.basename(path)}")
    
    return detected

def load_species_data(file_path):
    """加载单个物种数据"""
    with open(file_path) as f:
        kmers, array = read_kmers_optimized(f)
    return (kmers, array)

def init_worker(shared_species_data):
    """初始化子进程的全局变量"""
    global species_data
    species_data = shared_species_data

def process_line_parallel(line):
    """并行处理单行数据"""
    parts = line.strip().split("\t")
    if len(parts) < 5:
        return None  # 如果格式不对，返回None
    
    # 提取原始数据
    row_id = parts[0]
    chr_info = parts[1]
    start_pos = parts[2]
    end_pos = parts[3]
    original_seq = parts[4]
    
    # 获取其他列（如果有）
    additional_cols = parts[5:] if len(parts) > 5 else []
    
    # 计算反向互补序列
    asoseq = DNA_complement(original_seq)
    
    # 计算同源性
    homology_results = {}
    for species_name, (species_set, species_arr) in species_data.items():
        similarity = get_max_match_vectorized(asoseq, species_set, species_arr)
        homology_results[species_name] = f"{similarity:.2f}"
    
    # 构建结果行
    result_row = [
        row_id, chr_info, start_pos, end_pos, original_seq
    ] + additional_cols + [
        homology_results.get('crab_eating_macaque', '0.00'),
        homology_results.get('mouse', '0.00'),
        homology_results.get('rat', '0.00'),
        homology_results.get('pig', '0.00'),
        homology_results.get('rabbit', '0.00'),
        homology_results.get('Guinea_pig', '0.00')
    ]
    
    return result_row

if __name__ == "__main__":
    # 自动检测并加载物种数据
    species_files = auto_detect_species_files()
    species_data = {}
    
    if not species_files:
        raise ValueError("未检测到物种数据文件！请检查 /asodesigner/outfile/ 目录")
    
    print("\n加载物种数据:")
    for species_name, file_path in species_files.items():
        print(f"  正在加载 {species_name}...")
        species_data[species_name] = load_species_data(file_path)
        print(f"  - {species_name} 已加载: {len(species_data[species_name][0])} 个 k-mer")
    
    # 定义表头
    headers = [
        "ID", "chr", "start", "end", "ASO sequence", "RNase H score", "GC content", 
        "ASO MFE", "0 mismatch genes","1 mismatch genes", "2 mismatch genes", "3 mismatch genes", 
        "4 mismatch genes", "mismatch genes name", "ASO position",
        "Macaca fascicularis (crab-eating macaque) homology",
        "Mus musculus (mouse) homology",
        "Rattus norvegicus (rat) homology",
        "Sus scrofa (pig) homology",
        "Oryctolagus cuniculus (rabbit) homology",
        "Cavia porcellus (Guinea pig) homology"
    ]
    
    # 启动多进程池
    cpu_count = max(2, os.cpu_count() - 1)  # 保留1个核心给系统
    with Pool(processes=cpu_count, 
             initializer=init_worker,
             initargs=(species_data,)) as pool:
        
        # 读取输入文件
        with open(f"{uid}/result.txt", "r") as f_in:
            lines = f_in.readlines()
        
        # 动态调整chunksize
        total_lines = len(lines)
        chunk_size = max(1000, total_lines // (cpu_count * 10))  # 每worker处理10个chunk
        
        print(f"\n开始处理: {total_lines} 行, 使用 {cpu_count} 个进程, chunk大小 {chunk_size}")
        
        # 处理所有行
        processed = 0
        results = []
        for result in pool.imap_unordered(process_line_parallel, lines, chunksize=chunk_size):
            if result is not None:
                results.append(result)
            processed += 1
            if processed % 1000 == 0:
                print(f"已处理: {processed}/{total_lines} 行 ({processed/total_lines*100:.1f}%)")
    
    # 创建DataFrame并保存为Excel
    df = pd.DataFrame(results)
    
    # 确保列数匹配
    if len(df.columns) > len(headers):
        # 如果结果列数多于表头，截断多余的列
        df = df.iloc[:, :len(headers)]
    elif len(df.columns) < len(headers):
        # 如果结果列数少于表头，添加空列
        for i in range(len(df.columns), len(headers)):
            df[i] = ""
    
    # 设置列名
    df.columns = headers
    
    # 保存为Excel文件
    output_filename = f"{uid}/ASO_AllCandidates_{sys.argv[-2]}.xlsx"
    df.to_excel(output_filename, index=False)
    output_filename_txt = f"{uid}/ASO_AllCandidates_{sys.argv[-2]}.txt"
    df.to_csv(output_filename_txt, sep='\t', index=False)
    
    print(f"\n处理完成! 结果已保存到 {output_filename}")
