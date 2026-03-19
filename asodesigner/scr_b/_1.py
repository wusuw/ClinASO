import sys
import os
import bisect

output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

# Read files using context managers for automatic closing
with open(f"{uid}/gene.fa", "r") as f1, \
     open(f"{uid}/snp.gtf", "r") as f2, \
     open(f"{uid}/out1.txt", "w") as fo:

    # 1. Preprocess SNP positions - sort once
    snp_positions = sorted(int(line.split("\t")[2]) for line in f2)

    # 2. Read gene sequences efficiently
    gene_dict = {}
    current_seq = []
    current_key = None
    
    for line in f1:
        line = line.strip()
        if line.startswith(">"):
            if current_key is not None:
                gene_dict[current_key] = ''.join(current_seq)
            current_key = line.split(":")[1]
            current_seq = []
        else:
            current_seq.append(line)
    if current_key is not None:
        gene_dict[current_key] = ''.join(current_seq)

    # 3. Process each sequence with optimized splitting
    for region, seq in gene_dict.items():
        try:
            start, end = map(int, region.split("-"))
            seq_length = len(seq)
            
            # Get relevant SNPs using binary search
            left_index = bisect.bisect_left(snp_positions, start + 1)
            right_index = bisect.bisect_right(snp_positions, end - 1)
            region_snps = snp_positions[left_index:right_index]
            
            # Build split points (relative to sequence start)
            split_points = [0]
            split_points.extend(pos - start for pos in region_snps)
            split_points.append(seq_length)  # Use actual sequence length
            
            # Generate valid segments
            for i in range(1, len(split_points)):
                seg_length = split_points[i] - split_points[i-1]
                if seg_length >= 20:
                    seg_start = start + split_points[i-1]
                    seg_end = start + split_points[i] - 1  # Corrected end position
                    segment = seq[split_points[i-1]:split_points[i]]
                    fo.write(f">{seg_start}-{seg_end}\n{segment}\n")
                    
        except ValueError:
            # Skip malformed entries
            continue

# Retained unused function for compatibility
def DNA_complement1(sequence):
    comp_dict = {"A":"T", "T":"A", "G":"C", "C":"G",
                 "a":"t", "t":"a", "g":"c", "c":"g"}
    return ''.join(comp_dict[base] for base in sequence)[::-1]