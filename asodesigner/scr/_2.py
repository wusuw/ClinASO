import sys,os


nu_seq = int(sys.argv[1])           # bin_length
geneDirection = str(sys.argv[2])         # gene direaction
chrn = sys.argv[3]
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
f1 = open(f"{uid}/out1.txt","r")
fo = open(f"{uid}/out3.txt","w")     #change file


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

dic_seq = {}
for line in f1:
    if line.startswith(">"):
        name = line.strip().replace(">","")
        dic_seq[name] = ""
    else :
        dic_seq[name] += line.strip()
a = 0
for k,v in dic_seq.items():
    for i in range(len(v)-nu_seq+1):
        if geneDirection == "-":
            sequence = v[i:i+nu_seq].upper()
        if geneDirection == "+":
            sequence = DNA_complement1(v[i:i+nu_seq].upper())
        fo.write("ASO"+str(a)+"\t"+chrn+"\t"+str(int(k.split("-")[0])+i)+"\t"+str(int(k.split("-")[0])+i+nu_seq)+"\t"+sequence+"\n")
        a+=1
