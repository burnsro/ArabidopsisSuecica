out=open('ATside_forpi.csv', 'w')
for line in open('As_SNPs2020.Asue_At.Al2At.nonVars.BiSNPS.vcf'):
	if '##' in line:
		pass
	elif '#CHR' in line:
		names=line.strip().split()[9:]
		out.write('chr,pos,%s\n'%(','.join(names)))
	elif line.strip().split()[0] in ['Asue_scaffold1', 'Asue_scaffold2', 'Asue_scaffold3', 'Asue_scaffold4', 'Asue_scaffold5']:
		gen=line.strip().split()[9:]
		gen=[x.split(':')[0] for x in gen]
		gen_new=[]
		for g in gen:
			if g=='./.' or g =='./././.':
				gen_new.append('NA')
			else:
				g=[int(x) for x in g.split('/')]
				gen_new.append(str(sum(g)/float(len(g))))
		test=filter(lambda a: a != 'NA', gen_new)
		if len(set(test))>=1:
			gen_new=','.join(gen_new)
			out.write('%s,%s,%s\n'%(line.strip().split()[0], line.strip().split()[1], gen_new))
	else:
		break
