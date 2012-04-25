#coding: utf-8
import math
import sys
import numpy as np
sys.path.append("/Users/bingwang/VimWork/")
import lib.func.read as read
import lib.func.paml as paml
WORKPATH = "/Users/bingwang/VimWork/"

class Species():
    def __init__(self,name):
        self.name = name
        self.gene = {}
        f = open(WORKPATH+'db/'+self.name+"AA.fasta")
        for line in f:
            line = line.strip()
            if line[0] == ">":
                name = line[1:]
            else:
                self.gene[name] = line
        f.close()

def branch_calculte(filename):
    f=open(WORKPATH+filename)
    f.readline()
    paml_rate = {}
    for line in f:
        line = line.strip()
        paml_out = paml.read_line(line)
        sp0 = round(((paml_out.sp1_sp0 + paml_out.sp2_sp0) - paml_out.sp2_sp1)/2,4)
        sp1 = round(((paml_out.sp2_sp1 + paml_out.sp1_sp0) - paml_out.sp2_sp0)/2,4)
        sp2 = round(((paml_out.sp2_sp1 + paml_out.sp2_sp0) - paml_out.sp1_sp0)/2,4)
        if math.isnan(sp0) or math.isnan(sp1) or math.isnan(sp2):
            pass
        else:
            sp0 = sp0 if sp0 > 0 else 0.0001
            sp1 = sp1 if sp1 > 0 else 0.0001
            sp0_devide_sp1 = sp0/sp1
            sp0_devide_sp1 = 16 if sp0_devide_sp1 > 16 else sp0_devide_sp1
            sp0_devide_sp1 = 0.0625 if sp0_devide_sp1 < 0.0625 else sp0_devide_sp1
            paml_rate[paml_out.sp0name]=np.log2(sp0_devide_sp1)    #paml-rate structure
    return paml_rate

def draw_scatter():
    import lib.draw.scatter_hist_test as scatter
    kwal_sklu = []
    cgla_scer = []
    cgla_sbay = []
    agos_klac = []
    for name in Bg_Genes:
        kwal_sklu.append(Kwal_Sklu_Klac[Scer.orth_Kwal[name]])
        cgla_scer.append(Cgla_Scer_Scas[Scer.orth_Cgla[name]])
        cgla_sbay.append(Cgla_Sbay_Scas[Scer.orth_Cgla[name]])
        agos_klac.append(Agos_Klac_Sklu[Scer.orth_Agos[name]])
    scatter.scatter(kwal_sklu,cgla_scer,WORKPATH+"AftGalLoss/co_Kwal_Cgla(Scer).png",\
            "lg2(Kwal/Sklu) vs lg2(Cgla/Scer)")
    scatter.scatter(kwal_sklu,cgla_sbay,WORKPATH+"AftGalLoss/co_Kwal_Cgla(Sbay).png",\
            "lg2(Kwal/Sklu) vs lg2(Cgla/Sbay)")
    scatter.scatter(kwal_sklu,agos_klac,WORKPATH+"AftGalLoss/co_Kwal_Agos.png",\
            "lg2(Kwal/Sklu) vs lg2(Agos/Klac)")
    scatter.scatter(cgla_scer,cgla_sbay,WORKPATH+"AftGalLoss/co_Cgla(Scer)_Cgla(Sbay).png",\
            "lg2(Cgla/Scer) vs lg2(Cgla/Sbay)")
    scatter.scatter(cgla_scer,agos_klac,WORKPATH+"AftGalLoss/co_Cgla(Scer)_Agos.png",\
            "lg2(Cgla/Scer) vs lg2(Agos/Klac)")
    scatter.scatter(cgla_sbay,agos_klac,WORKPATH+"AftGalLoss/co_Cgla(Sbay)_Agos.png",\
            "lg2(Cgla/Sbay) vs lg2(Agos/Klac)")

def go_term_query():
    f = open(WORKPATH+"/Agos/term_query_Genes.txt","w")
    f.write("Bg_Genes\tlog2(Agos/Sklu)\tlog2(Agos/Kwal)\tlog2(Agos/Klac)\tlog2(Agos/Klac)\n")
    for name in Bg_Genes:
        f.write(name+"\t"+str(Agos_Sklu_Kwal[Scer.orth_Agos[name]])\
               +"\t"+str(Agos_Kwal_Scer[Scer.orth_Agos[name]])  \
               +"\t"+str(Agos_Klac_Kwal[Scer.orth_Agos[name]]) \
               +"\t"+str(Agos_Klac_Sklu[Scer.orth_Agos[name]])+"\n")

def p_value(n_hit,n,p):
    pv = 0
    for i in range(n_hit,n+1,1):
        pv += ((math.factorial(n)/(math.factorial(i)* \
                math.factorial(n-i)))*(p**i)*((1-p)**(n-i)))
    return pv


########main##########
Agos_Sklu_Kwal = branch_calculte("Agos/Agos_Sklu_Kwal.txt")
Agos_Kwal_Scer = branch_calculte("Agos/Agos_Kwal_Scer.txt")
Agos_Klac_Kwal = branch_calculte("Agos/Agos_Klac_Kwal.txt")
Agos_Klac_Sklu = branch_calculte("Agos/Agos_Klac_Sklu.txt")

Scer = Species("Scer")
Scer.orth_Agos = read.one_orth("Scer","Agos")

Bg_Genes = []
for name in Scer.gene:
    try:
        Agos_Sklu_Kwal[Scer.orth_Agos[name]]
        Agos_Kwal_Scer[Scer.orth_Agos[name]]
        Agos_Klac_Kwal[Scer.orth_Agos[name]]
        Agos_Klac_Sklu[Scer.orth_Agos[name]]
        Bg_Genes.append(name)
    except:
        continue

draw_scatter()
go_term_query()

