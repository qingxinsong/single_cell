#This script will generate normalized expression level of each gene in each cell using slope and intercept from ERCC level data
import sys
import re
import os
import math
IN_GLM=open('GLM_param.txt')
glm={}
tmp=IN_GLM.readline()
for line in IN_GLM:
	tmp=line.split('\t')
	glm[tmp[0]]=[float(tmp[1]),float(tmp[2])]
dir=os.listdir('raw_data/gene_raw_reads')
gene={}
tissue=[]
for i in sorted(dir):
	IN=open('raw_data/gene_raw_reads/'+i)
	cell=i.split('_')[0]
	tissue.append(cell)
	for line in IN:
		tmp=line.strip().split('\t')
		if float(tmp[1])==0:
			gene_level=float(tmp[1])
		else:
			gene_level=2**((math.log(float(tmp[1]))-glm[cell][0])/glm[cell][1])
		if tmp[0] in gene:
			gene[tmp[0]].append(gene_level)
		else:
			gene[tmp[0]]=[gene_level]
	IN.close()
OUT=open('Normalized_gene_abundance.txt',mode='w')
OUT.write('gene_id'+'\t'+'\t'.join(tissue)+'\n')
for i in sorted(gene.keys()):
	OUT.write(i+'\t'+'\t'.join([str(m) for m in gene[i]])+'\n')