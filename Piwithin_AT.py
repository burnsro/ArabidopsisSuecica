from itertools import combinations
from numpy import *
import vcf
import pysam
import tabix
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help='VCF.gz, tabixed', default="")
parser.add_option("-c", "--chr", dest="chr", help='', default="")                  
parser.add_option("-s", "--st", dest="st", help='', default="")                  
parser.add_option("-e", "--end", dest="end", help='', default="")                  
parser.add_option("-n", "--num", dest="num", help='line number in the interval file', default="")
parser.add_option("-o", "--out", dest="out", help='output directory', default="")

(options,args) = parser.parse_args()

vcf_reader=vcf.Reader(filename='%s'%options.vcf)


L,Lal,P=0,0,0
L=int(options.end)-int(options.st)
for record in vcf_reader.fetch(options.chr, int(options.st), int(options.end)):
	gen=[]
	for sample in record.samples:
		gen.append(sample['GT'])
	gen=gen[9:39]
	gen=list(filter(lambda a: a != './.', gen))
	if (gen.count('0/1')+gen.count('1/0'))<float(len(gen))/3:
	
		gen=list(''.join([x[0::2] for x in gen]))
		if len(set(gen))==1:
			Lal+=1
			#print (gen)
		elif record.var_type=='snp' and len(list(set(gen)))==2:
			Lal+1
			if len(gen)<15:##allow for no missing data only
				pass
			else:
				#print (gen)
				Ncomp,N=0,0
				for s1,s2 in combinations(gen,2):
					Ncomp+=1
					if s1!=s2:
						N+=1
				P+=N/float(Ncomp)
		else:
			pass
out=open('%s/gene%s_%s_%s_%s.piwithinAT'%(options.out,options.num,options.chr, options.st, options.end), 'w')
out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(options.num,options.chr, options.st, options.end, L,Lal,P,P*100/float(Lal)))
out.close()
