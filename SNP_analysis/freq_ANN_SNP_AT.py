out=open('ATside_freqV2.txt', 'w')
out.write('freqAT,freqAS,ann1,ann2\n')
for line in open('AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPsONLY.vcf_clean.ann.ALpolarisedATsub.Asue_scaffoldAT.POL.ANN'):
	if 'coord' in line:
		names=line.strip().split()
		#print (names[4:13]+names[17:23])
	else:
		gen1, gen2=line.strip().split()[4:13]+line.strip().split()[17:23], line.strip().split()[31:46]
		#if gen2.count('0.5')<:
		#	print (gen2.count('0.5'), gen2)
		gen1, gen2=list(filter(lambda a: a != 'NA', gen1)),list(filter(lambda a: a != 'NA', gen2))
		if gen2.count('0.5')<float(len(gen2))/3:
			freq1, freq2=sum([float(x) for x in gen1])/len(gen1), sum([float(x) for x in gen2])/len(gen2)
			eff=line.strip().split()[-1].split('|')
			#print (eff[1], eff[2])
			if eff[2]!='MODIFIER' and freq1!=freq2:
				out.write('%s,%s,%s,%s\n'%(freq1, freq2, eff[1], eff[2]))
out.close()
