import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor
import re
import traceback
import multiprocessing
import logging
import pickle

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 配置物种信息（俗名 -> 目录路径）
SPECIES_INFO = {
    "mouse": {
        "gtf": "/asodesigner/reference/mouse/GRCm39_genomic.gtf",
        "genome": "/asodesigner/reference/mouse/GRCm39_genomic.fna"
    },
    "rat": {
        "gtf": "/asodesigner/reference/rat/GRCr8_genomic.gtf",
        "genome": "/asodesigner/reference/rat/GRCr8_genomic.fna"
    },
    "guinea_pig": {
        "gtf": "/asodesigner/reference/guinea_pig/GCF_034190915.1_mCavPor4.1_genomic.gtf",
        "genome": "/asodesigner/reference/guinea_pig/GCF_034190915.1_mCavPor4.1_genomic.fna"
    },
    "pig": {
        "gtf": "/asodesigner/reference/pig/GCF_000003025.6_Sscrofa11.1_genomic.gtf",
        "genome": "/asodesigner/reference/pig/GCF_000003025.6_Sscrofa11.1_genomic.fna"
    },
    "crab_eating_macaque": {
        "gtf": "/asodesigner/reference/crab_eating_macaque/genomic.gtf",
        "genome": "/asodesigner/reference/crab_eating_macaque/genomic.fna"
    },
    "rabbit": {
        "gtf": "/asodesigner/reference/rabbit/GCF_964237555.1_mOryCun1.1_genomic.gtf",
        "genome": "/asodesigner/reference/rabbit/GCF_964237555.1_mOryCun1.1_genomic.fna"
    },
    "human": {
        "gtf": "/asodesigner/reference/human/genomic.gtf",
        "genome": "/asodesigner/reference/human/genomic.fna"
    }
}

# 使用进程安全的缓存管理器
class CacheManager:
    def __init__(self):
        # 使用普通字典而不是Manager.dict()
        self.genome_cache = {}
        self.gtf_cache = {}
    
    def get_genome_cache(self):
        return self.genome_cache
    
    def get_gtf_cache(self):
        return self.gtf_cache

# 解析物种映射文件
def load_species_mapping(species_file):
    species_map = {}
    with open(species_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                scientific_name, taxid, common_name = parts[:3]
                species_map[common_name] = taxid
    return species_map

# 解析同源文件，获取目标物种的基因ID
def get_ortholog_gene_ids(ortholog_file, human_gene_id, target_taxids, species_map):
    orthologs = {}
    with open(ortholog_file, 'r') as f:
        for line in f:
            cols = line.strip().split()
            if len(cols) < 5:
                continue
                
            if cols[1] == human_gene_id:
                target_taxid = cols[3]
                gene_id = cols[4]
                if target_taxid in target_taxids:
                    # 查找匹配的俗名
                    for common_name, taxid in species_map.items():
                        if taxid == target_taxid:
                            orthologs[common_name] = gene_id
                            break
    return orthologs

# 使用正则表达式快速提取 GeneID
gene_id_pattern = re.compile(r'GeneID:(\d+)')

# 从GTF文件中获取基因位置信息（缓存）
def get_gene_location_cached(gtf_file, target_gene_id, gtf_cache):
    if gtf_file not in gtf_cache:
        cache = {}
        try:
            with open(gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    if len(cols) < 9:
                        continue
                    if cols[2] == 'gene':
                        match = gene_id_pattern.search(cols[8])   
                        if match:
                            gene_id = match.group(1)
                            chrom = cols[0]
                            start = int(cols[3])
                            end = int(cols[4])
                            strand = cols[6]
                            cache[gene_id] = (chrom, start, end, strand)
            gtf_cache[gtf_file] = cache
        except Exception as e:
            logger.error(f"Error reading GTF file {gtf_file}: {str(e)}")
            return None
    return gtf_cache[gtf_file].get(target_gene_id)

# 从基因组中提取序列（缓存）
def extract_sequence_cached(genome_file, chrom, start, end, strand, genome_cache):
    if genome_file not in genome_cache:
        try:
            genome_cache[genome_file] = SeqIO.index(genome_file, "fasta")
        except Exception as e:
            logger.error(f"Error indexing genome file {genome_file}: {str(e)}")
            return None
            
    genome_index = genome_cache[genome_file]
    try:
        # 处理可能的染色体前缀（如'chr1' vs '1'）
        chrom_variants = [chrom, f"chr{chrom}", chrom.replace("chr", "")]
        for chrom_var in chrom_variants:
            if chrom_var in genome_index:
                seq = genome_index[chrom_var][start - 1:end]
                if strand == '-':
                    seq = seq.reverse_complement()
                return str(seq.seq)
        return None
    except KeyError:
        return None

# 更可靠的FASTA写入函数
def write_fasta_record(record, output_file):
    """手动写入FASTA记录，避免SeqIO.write的问题"""
    try:
        # 创建目录路径（如果不存在）
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        
        # 构建FASTA格式内容
        header = f">{record.id}"
        if record.description and record.description != record.id:
            header += f" {record.description}"
        
        # 格式化序列（每行80个字符）
        seq_str = str(record.seq)
        formatted_seq = "\n".join(seq_str[i:i+80] for i in range(0, len(seq_str), 80))
        
        # 写入文件
        with open(output_file, 'w') as f:
            f.write(f"{header}\n{formatted_seq}\n")
            
        return True
    except Exception as e:
        logger.error(f"Error writing to {output_file}: {str(e)}")
        traceback.print_exc()
        return False

# 单个物种处理函数
def process_species(args):
    """处理单个物种的序列提取任务"""
    species, human_gene_id, species_map, orthologs, output_dir = args
    
    try:
        # 创建本地缓存
        genome_cache = {}
        gtf_cache = {}
        
        # 获取当前物种的基因ID
        gene_id = orthologs.get(species)
        if not gene_id:
            logger.warning(f"No ortholog found for {species}")
            return False

        # 获取文件路径
        species_info = SPECIES_INFO.get(species)
        if not species_info:
            logger.error(f"No GTF/genome info for {species}")
            return False
            
        gtf_file = species_info["gtf"]
        genome_file = species_info["genome"]

        # 获取基因位置
        location = get_gene_location_cached(gtf_file, gene_id, gtf_cache)
        if not location:
            logger.error(f"Gene {gene_id} not found in {species} GTF file")
            return False

        # 提取序列
        chrom, start, end, strand = location
        sequence = extract_sequence_cached(genome_file, chrom, start, end, strand, genome_cache)
        if not sequence:
            logger.error(f"Chromosome {chrom} not found in {species} genome or empty sequence extracted")
            return False

        # 创建SeqRecord对象
        record_id = f"{species}_{gene_id}"[:200]  # 限制长度
        description = f"{species} ortholog of {human_gene_id}"[:200].replace('\n', '')
        
        record = SeqRecord(
            Seq(sequence),
            id=record_id,
            description=description
        )

        # 设置输出路径
        output_file = os.path.join(output_dir, f"{species}.fasta")
        logger.info(f"Processing {species}, writing to: {output_file}")

        # 使用自定义写入函数
        if write_fasta_record(record, output_file):
            logger.info(f"Successfully wrote FASTA for {species}")
            return True
        else:
            logger.error(f"Failed to write FASTA for {species}")
            return False

    except Exception as e:
        logger.error(f"Error processing {species}: {str(e)}")
        traceback.print_exc()
        return False

# 主函数
def main(human_gene_id, species_file, ortholog_file, output_dir):
    # 加载物种映射
    species_map = load_species_mapping(species_file)
    target_species = list(species_map.keys())
    target_taxids = [species_map[sp] for sp in target_species]

    # 获取同源基因
    orthologs = get_ortholog_gene_ids(ortholog_file, human_gene_id, target_taxids, species_map)
    orthologs["human"] = human_gene_id
    
    # 确保输出目录存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备任务参数
    tasks = []
    for species in target_species + ["human"]:
        # 跳过没有ortholog的物种（除了human）
        if species != "human" and species not in orthologs:
            logger.warning(f"Skipping {species}: no ortholog found")
            continue
            
        tasks.append((
            species, 
            human_gene_id, 
            species_map, 
            orthologs, 
            output_dir
        ))
    
    # 设置最大工作进程数（避免过多进程导致资源耗尽）
    max_workers = min(4, os.cpu_count() or 1)
    
    # 并行处理
    success_count = 0
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_species, task) for task in tasks]
        
        # 等待所有任务完成并收集结果
        for future in futures:
            try:
                result = future.result()
                if result:
                    success_count += 1
            except Exception as e:
                logger.error(f"Task failed: {str(e)}")
                traceback.print_exc()

    logger.info(f"Processing completed! Success: {success_count}/{len(tasks)} species")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract ortholog gene sequences")
    parser.add_argument("human_gene_id", help="Human gene ID (e.g., 255738)")
    parser.add_argument("species_file", help="Species mapping file (3 columns: scientific_name, taxid, common_name)")
    parser.add_argument("ortholog_file", help="Ortholog relationships file")
    parser.add_argument("output_dir", help="Output directory for FASTA files")
    args = parser.parse_args()

    main(args.human_gene_id, args.species_file, args.ortholog_file, args.output_dir)
