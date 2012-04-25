#coding: utf-8
import os
from Bio.Seq import Seq
'''
1.V.polyspora Position 1
2.T.phaffii Position 1
3.X.blattae Position 1 (draft annotation)
4.N.dairenensis Position 1
5.N.castellii Position 1
6.S.castellii Position 1 (old annotation)
7.X.naganishii Position 2 (draft annotation)
8.X.africana Position 2 (draft annotation)
9.C.glabrata Position 1
10.S.bayanus Position 1
11.S.cerevisiae Position 1
12.Anestral Gene Order
13.Z.rouxii
14.T.delbrueckii
15.K.lactis
16.E.gossypii 
17.E.cymbalariae
18.L.kluyveri
19.L.thermotolerans
20.L.waltii
21.S.cerevisiae Position 2
22.S.bayanus Position 2
23.C.glabrata Position 2
24.X.africana Position 2 (draft annotation)
25.X.naganishii Position 2 (draft annotation)
26.S.castellii Position 2 (old annotation)
27.N.castellii Position 2
28.N.dairenensis Position 2
29.X.blattae Position 2 (draft annotation)
30.T.phaffii Position 2
31.V.polyspora Position 2
'''
pillar_tab_file = "/Users/bingwang/VimWork/db/Pillars.tab"
f = open(pillar_tab_file)
gene_tree = {}
aln_groups = {}
for i,line in enumerate(f):
    line = line.strip()
    genename = line.split("\t")
    aln_groups[i] = [name for name in genename if name != "---"]
    evolve_pattern = ""
    for a in genename:
        has = "F" if a == "---" else "T"
        evolve_pattern += has
    try:
        gene_tree[evolve_pattern].append(i)
    except:
        gene_tree[evolve_pattern] = [i]

genome_tab_path = "/Users/bingwang/VimWork/db/genome_tab/"
gene2sp = {}
gene2pos = {}
tab_files = [a for a in os.listdir(genome_tab_path) if a.endswith(".tab")]
for tab_file in tab_files:
    sp_name = tab_file[:tab_file.find("_")]
    if sp_name != "Ancestor":
        f = open(genome_tab_path + tab_file)
        for line in f:
            elements = line.split("\t")
            gene2sp[elements[0]] = sp_name
            gene2pos[elements[0]] = [elements[5],elements[7]]
    else:
        f = open(genome_tab_path + tab_file)
        for line in f:
            name = line.split("\t")[0]
            gene2sp[name] = sp_name

sequence_path = "/Users/bingwang/VimWork/db/sequence/"
scaffold2sq = {}
sp2scaffold = {}
sequence_files = [a for a in os.listdir(sequence_path) if a.endswith(".fsa")]
for sequence_file in sequence_files:
    sp_name = sequence_file[:sequence_file.find("_")]
    f = open(sequence_path + sequence_file)
    for line in f:
        if line[0] == ">":
            name = line[1:line.find(" ")]
            try:
                sp2scaffold[sp_name]
                if sp2scaffold[sp_name] != name[:name.rfind("_")+1]:
                    print "Error", sp2scaffold[sp_name], "!= ", name[:name.rfind("_")+1]
            except:
                sp2scaffold[sp_name] = name[:name.rfind("_")+1]
            scaffold2sq[name] = ""
        else:
            scaffold2sq[name] += line.strip()

print "Finish Reading"

def get_sequence(sequence,pos,comp=None):
    def reverse(s):
        a = ""
        for c in s:
            if c.upper() == "A":
                a = "T" + a
            elif c.upper() == "T":
                a = "A" + a
            elif c.upper() == "C":
                a = "G" + a
            elif c.upper() == "G":
                a = "C" + a
        return a
    if comp == None:
        comp = False
    seq = ""
    for sites in pos:
        [start,end] = sites.split("..")
        start = int(start)
        end = int(end)
        if comp == False:
            seq += sequence[start-1:end]
        else:
            seq = reverse(sequence[start-1:end]) + seq
    return seq

gene2sq = {}
for i,name in enumerate(gene2sp):
    if name[:3] != "Anc":
        scaf_name = sp2scaffold[gene2sp[name]] + gene2pos[name][0]
        comp = True if gene2pos[name][1].count("complement") > 0 else False
        pos = gene2pos[name][1][gene2pos[name][1].find("(")+1:\
                gene2pos[name][1].rfind(")")].split(",")
        gene2sq[name] = get_sequence(scaffold2sq[scaf_name],pos,comp)
        print i * 1.0 / len(gene2sp)

f = open("/Users/bingwang/VimWork/db/YGOB_data.fsa","w")
for name in gene2sq:
    f.write(">"+name+"\n")
    f.write(gene2sq[name]+"\n")


