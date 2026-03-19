# use SNP split sequenc
import sys,os
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
f1 = open(f"{uid}/gene2.gtf","r")
fo = open(f"{uid}/exon.txt","w")


def DNA_complement1(sequence):              #反向互补序列获取
    # 构建互补字典
    comp_dict = {
        "A":"T",
        "T":"A",
        "G":"C",
        "C":"G",
        "a":"t",
        "t":"a",
        "g":"c",
        "c":"g",
    }
    #求互补序列
    sequence_list = list(sequence)
    sequence_list = [comp_dict[base] for base in sequence_list]
    string = ''.join(sequence_list)[::-1]
    return string
lis1 = []
for line in f1:
    seq =line.strip().split("\t")
    if seq[2] == "CDS":
        a = seq[0]+"\t"+str(int(seq[3])-16)+"\t"+str(int(seq[4])+16)
        lis1.append(a)
lis1 = set(lis1)
for i in  lis1:
    fo.write(i+"\n")
