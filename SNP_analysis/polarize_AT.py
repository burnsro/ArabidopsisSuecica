import sys
from optparse import OptionParser
parser=OptionParser()
parser.add_option('-c', '--chr', dest='chr', help='chromosome number', default='')
(options, args)=parser.parse_args()

chr=options.chr

all_coord,rev_coord, rev_gen=[],[], []
for line in open('As_SNPs2020.Asue_At.Al2At.BiSNPS.vcf.ANN.csv.ALfixed2'):
	line=line.strip().split('\t')
	all_coord.append('_'.join(line[:2]))
	if float(line[2])==1.0 and line[0]==chr:
		rev_coord.append('_'.join(line[:2]))
		rev_gen.append(line[2])
		#print (int(line[2][-1]), int(chr[-1]))
	elif float(line[0].split('scaffold')[-1])>float(chr.split('scaffold')[-1]):
		break
	else:
		pass
#print (len(rev_coord))

out=open('As_SNPs2020.Asue_At.Al2At.BiSNPS.vcf.ANN.csv.%s.POL'%chr, 'w')
for line in open('As_SNPs2020.Asue_At.Al2At.BiSNPS.vcf.ANN.csv'):
	mychr=line.strip().split(',')[0]
	coord='_'.join(line.strip().split(',')[:2])
	if '##' in line:
		pass
	elif 'chr' in line:
		names=line.strip().split(',')[2:46]+line.strip().split(',')[87:88]
		#print (names)
		#break
		out.write('coord\t%s\trev\tann\n'%('\t'.join(names)))
		#print (names[:30], names[30:])
		#break
	elif mychr==chr and coord in all_coord:
		gen=line.strip().split(',')[2:46]+line.strip().split(',')[87:88]
		#gen=[x.split(":")[0] for x in gen]
		genAT,genAS=gen[:30], gen[30:]
		if genAT.count('NA')<12 and genAS.count('NA')<float(len(genAS))/10:
			coord='_'.join(line.strip().split(',')[:2])
			if coord in rev_coord:
				ind=rev_coord.index(coord)
				rev=float(rev_gen[ind])
				#print (coord, rev)
				#break
				#print (genAT)
				genAT_pol=[]
				for g in genAT:
					if g=='NA':
						genAT_pol.append('NA')
					else:
						#g=[int(x) for x in g.split('/')]
						#print(g)
						#break
						g=abs(float(g) - rev)
						g=sum([float(g)])/(len([float(g)]))
						genAT_pol.append(g)
				genAT_pol='\t'.join([str(x) for x in genAT_pol])
				#print (genAT_pol)
				genAS_pol=[]
				for g in genAS:
					if g=='NA':
						genAS_pol.append('NA')
					else:
						#g=[int(x) for x in g.split('/')]
						g=abs(float(g) - rev)
						g=sum([float(g)])/(len([float(g)]))
						genAS_pol.append(g)
				genAS_pol='\t'.join([str(x) for x in genAS_pol])
			else:
				rev='-'
				genAT_pol=[]
				for g in genAT:
					if g=='NA':
						genAT_pol.append('NA')
					else:
						g=sum([float(g)])/float(len([float(g)]))
						genAT_pol.append(g)
				genAT_pol='\t'.join([str(x) for x in genAT_pol])
				genAS_pol=[]
				for g in genAS:
					if g=='NA':
						genAS_pol.append('NA')
					else:
						g=sum([float(g)])/float(len([float(g)]))
						genAS_pol.append(g)
				genAS_pol='\t'.join([str(x) for x in genAS_pol])
			eff=line.strip().split('|')[1:3]
			#if 'ANN' in eff:
			out.write('%s\t%s\t%s\t%s\t%s\n'%(coord, genAT_pol, genAS_pol,rev, eff))
			#print (line[0])
	elif int(mychr.split('scaffold')[-1])>int(chr.split('scaffold')[-1]):
		break
	else:
		pass
out.close()
