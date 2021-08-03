import sys
from optparse import OptionParser
from numpy import unique
# pipeline params
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="combined vcf file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")
(options, args) = parser.parse_args()

out=open('%s/%s.ATfixed'%(options.dir, options.vcf), 'w')
for line in open('%s/%s'%(options.dir, options.vcf)):
	if 'chr' in line:
		#line=line.strip().split(',')
		#gen=line[2:46]+line[87:88]
		#print(gen)
		#break
		pass
	else:
		gen=line.strip().split(',')[2:32]#outgroup
		#print (gen.count("NA"))
		#break
		gen=filter(lambda a: a != None, gen)
		gen=filter(lambda a: a != './././.', gen)
		gen=filter(lambda a: a != 'NA', gen)
		gen2=line.strip().split(',')[32:48] + line.strip().split(',')[50:59] + line.strip().split(',')[98:99] #biallelic SNP in As and At
		gen2=filter(lambda a: a != None, gen2)
		gen2=filter(lambda a: a != './././.', gen2)	
		if len(set(gen))==1 and len(set(gen2))>1:	
			gen=line.strip().split(',')[2:32]
			gen2=line.strip().split(',')[32:48] + line.strip().split(',')[50:59] + line.strip().split(',')[98:99] 
			gen_new=unique(list(gen))
			line2=line.strip().split(',')
			if gen.count("NA")==0 and gen2.count("NA")==0 and len(set(gen))==1 and len(set(gen2))>1:
				out.write('%s\t%s\n'%('\t'.join(line2[0:2]), '\t'.join(gen_new)))
			else:
				pass
out.close()
