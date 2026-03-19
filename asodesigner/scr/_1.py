# use SNP split sequence
import sys,os
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

f1 = open(f"{uid}/gene.fa","r")
f2 = open(f"{uid}/snp.gtf","r")
fo = open(f"{uid}/out1.txt","w")


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
for line in f2:
    seq = line.strip().split("\t")
    lis1.append(seq[2])


dic = {}
for line in f1:
    line = line.strip()
    if line.startswith(">"):
        name =line.split(":")[1]
        dic[name]= ""
    else:
        dic[name] = line

for k,v in dic.items():
    star = k.split("-")[0]
    end = k.split("-")[1]
    lis2 = []
    lis2.append(0)
    for i in lis1:
        if int(star) < int(i) < int(end):
           lis2.append(int(i)-int(star))
    lis2.append(int(end)-int(star))

    a = 0
    for j in range(len(lis2)):
        if a >0:
            name = ">"+str(int(star)+lis2[j-1])+"-"+str(int(star)+lis2[j]) 
            sequence = v[lis2[j-1]:lis2[j]]
            if len(sequence)>=20:
                fo.write(name+"\n"+sequence+"\n")
        a += 1
