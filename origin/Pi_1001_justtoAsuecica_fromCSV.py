#Get Pi of each 1001 accession just to Asuecica using a CSV file created from the VCF
from optparse import OptionParser
import itertools
parser = OptionParser()
parser.add_option("-c", "--csv", dest="csv", help="combined csv file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")

(options, args) = parser.parse_args()

out=open('%s/%s.distance'%(options.dir, options.csv), 'w')

from itertools import combinations

for line in open('%s/%s'%(options.dir,options.csv),'r'):
	line=line.strip().split(',')
	if 'chr' in line:
		header=line
		accAt=line[2:972]
		accAs=line[972:]
		#print(accAs)
		#break
		comb=list(itertools.product(accAt, accAs))	
		totalL,dist=[0]*len(comb),[0]*len(comb)
	else:
		gen=line[2:]
		genANC=line[2:972] 
		genAS=line[972:]
		count=0
		num_05 = sum(1 for g in genAS if g != 'NA' and float(g) == 0.5)
		num_05AT = sum(1 for g in genANC if g != 'NA' and float(g) == 0.5)
		if num_05 > 5 or num_05AT >= 0.05 * len(genANC):
			continue
		for i,j in list(itertools.product(accAt, accAs)):
			ind_i, ind_j=accAt.index(i),accAs.index(j)
			gen_i,gen_j=genANC[ind_i],genAS[ind_j]
			if gen_i=='NA' or gen_j=='NA':
				pass
			elif gen_i==gen_j:
				totalL[count]+=1
			else:
				dist[count]+=abs(float(gen_i)-float(gen_j))
			count+=1
out.write('acc1\t%acc2\taln\ttotal_dist\n')
for i in range(0, len(comb)):
	out.write('%s\t%s\t%s\t%s\n'%(list(comb[i])[0], list(comb[i])[1], totalL[i], dist[i]))
out.close()

