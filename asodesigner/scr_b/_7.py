import sys,os
output_dir = "/asodesigner/outfile/"
uid = sys.argv[-1]
f1 = open(f"{uid}/out6.txt","r")
fo = open(f"{uid}/out7.txt","w")

for line in f1:
    seq = line.split("\t")
    if  "N/A" in seq[-1]:
        continue
    fo.write(line)
