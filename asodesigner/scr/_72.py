import sys,os


fw = sys.argv[1]
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
f1 = open(f"{uid}/exonfasta.fa","r")
f2 = open(f"{uid}/out7.txt","r")
fo = open(f"{uid}/result.txt","w")

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
        "N":"N",}
    #求互补序列
    sequence_list = list(sequence)
    sequence_list = [comp_dict[base] for base in sequence_list]
    string = ''.join(sequence_list)[::-1]
    return string
dic_seq = {}
for line in f1:
    if line.startswith(">"):
        name = line.strip().replace(">","")
        dic_seq[name] = ""
    else :
        dic_seq[name] += line.strip()
for line in f2 :
    if fw == "+":
        a = 0
        seq = line.strip().split("\t")
        dna_seq = DNA_complement1(seq[4])
        for k,v in dic_seq.items():
            v = v.upper()
            if dna_seq in v:
                fo.write(line.strip()+"\t"+"exon"+"\n")
                a = 1
                break
        if a == 0:
            fo.write(line.strip()+"\t"+"intron"+"\n")
    if fw == "-":
        a = 0
        seq = line.strip().split("\t")
        dna_seq = seq[1]
        for k,v in dic_seq.items():
            v = v.upper()
            if dna_seq in v:
                fo.write(line.strip()+"\t"+"exon"+"\n")
                a = 1
                break
        if a == 0:
            fo.write(line.strip()+"\t"+"intron"+"\n")        
