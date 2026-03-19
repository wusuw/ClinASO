import sys,os
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]

f1 =open(f"{uid}/Gresult.txt","r")          #G.fa
f2 = open(f"{uid}/out4.txt","r")         #ASO_prdeca.txt
fo = open(f"{uid}/out5.txt","w")


nu = 1
dic_g = {}
for line in f1:
    if nu%3 == 1:
        aso_name =  line.strip().replace(">","")
    if nu%3 == 0:
        g = line.strip().split("(")[-1].split(")")[0].strip()
        dic_g[aso_name] = g
    nu += 1
for line in f2:
    try:
        seq = line.split("\t")
        v = float(dic_g.get(seq[0]))
        fo.write(line.strip()+"\t"+str(v)+"\n")
    except:
        continue
