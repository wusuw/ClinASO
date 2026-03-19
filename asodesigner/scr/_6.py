import sys,os
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
aso_len = int(sys.argv[-2])
targetgene = str(sys.argv[-3])
f1 = open(f"{uid}/human_blstn_result.txt","r")       #blast result
f2 = open(f"{uid}/out5.txt","r")       
fo = open(f"{uid}/out6.txt","w")       

def countnum(l):
    d = {}
    for i in l:
        if i not in d:
            d[i] = 1
        else:
            d[i] += 1
    return d

dic_genename = {}
dic_offgene = {}
for line in f1:
    seq = line.strip().split("\t")
    gene = str(seq[3])
    offgene = seq[1].split("_")[0]
    if offgene == targetgene:
        continue
    if seq[0] in dic_offgene.keys():
        dic_offgene[seq[0]].append(offgene)
    if seq[0] not in dic_offgene.keys():
        dic_offgene[seq[0]] = []
        dic_offgene[seq[0]].append(offgene)
    if seq[0] in dic_genename.keys():
        dic_genename[seq[0]].append(gene)
    if seq[0] not in dic_genename.keys():
        dic_genename[seq[0]]=[]
        dic_genename[seq[0]].append(gene)
for line in f2:
    seq = line .strip().split("\t")
    offagene = str("/".join(dic_offgene.get(seq[0],"NA")))
    offTgene = dic_genename.get(seq[0],"NA")
    offgenenum =countnum(offTgene)
    fo.write(line.strip()+"\t"+str(offgenenum.get(str(aso_len),"0"))+"\t"+str(offgenenum.get(str(aso_len-1),"0"))+"\t"+str(offgenenum.get(str(aso_len-2),"0"))+"\t"+str(offgenenum.get(str(aso_len-3),"0"))+"\t"+str(offgenenum.get(str(aso_len-4),"0"))+"\t"+offagene+"\n")
fo.close()
