from optparse import OptionParser
from numpy import unique
# pipeline params
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="combined vcf file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")
(options, args) = parser.parse_args()

out=open('%s/%s.ALfixed2'%(options.dir, options.vcf), 'w')
for line in open('%s/%s'%(options.dir, options.vcf)):
	if 'chr' in line:
		line=line.strip().split(',')
		#gen=line[48:74]
		#print(gen)
		#break
		#pass
	else:
		line=line.strip().split(',')
		gen=line[46:87] #outgroup 26 Alyratas
		#print (gen)
		#break
		gen=filter(lambda a: a != 'NA', gen)
		gen2=line[2:46] + line[87:88] #biallelic SNP in As and At
		gen2=filter(lambda a: a != 'NA', gen2)
		if len(set(gen))==1 and len(set(gen2))>1:		
			gen=line[46:87]
			gen2=line[2:46] + line[87:88]
			gen_new=unique(list(gen))
			if gen.count('NA')<20 and gen2.count('NA')<20:
				out.write('%s\t%s\n'%('\t'.join(line[0:2]), '\t'.join(gen_new)))
			else:
				pass
out.close()
