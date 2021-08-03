from optparse import OptionParser

# pipeline params
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="combined vcf file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")

(options, args) = parser.parse_args()                                                    

out=open('%s/%s.csv'%(options.dir, options.vcf), 'w')
for line in open('%s/%s'%(options.dir, options.vcf)):
	if '##' in line:
		pass
	elif '#CHR' in line:
		names=line.strip().split()[9:97]
		#print names
		#break
		out.write('chr,pos,%s,ann\n'%(','.join(names)))
	elif line.strip().split()[0] in ['Asue_scaffold1','Asue_scaffold2','Asue_scaffold3','Asue_scaffold4','Asue_scaffold5']:
		gen=line.strip().split()[9:97]
		gen=[x.split(':')[0] for x in gen]
		gen_new=[]
		for g in gen:
			if g=='./.' or g =='./././.':
				gen_new.append('NA')
			else:
				g=[int(x) for x in g.split('/')]
				gen_new.append(str(sum(g)/float(len(g))))
		gen_new=','.join(gen_new)
		out.write('%s,%s,%s,%s\n'%(line.strip().split()[0], line.strip().split()[1], gen_new,line.strip().split()[7]))
	else:
		pass
