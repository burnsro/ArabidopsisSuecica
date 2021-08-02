out=open('/groups/nordborg/projects/suecica/015SNPs/AApolarize/AAside_freq.txt', 'w')
out.write('freqAA,freqAS,ann1,ann2\n')
for line in open('/groups/nordborg/projects/suecica/015SNPs/AApolarize/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.biallelicSNPsONLY.vcf_clean.ann.ATpolarisedAAsub.ann.Asue_scaffoldAA.POL'):
	if 'coord' in line:
		names=line.strip().split()
	else:
		gen1, gen2, gen3 =line.strip().split()[1:16], line.strip().split()[16:25], line.strip().split()[58:64]+line.strip().split()[69:75]+line.strip().split()[82:83]+line.strip().split()[93:94]+line.strip().split()[95:97]+line.strip().split()[98:99]
		#print (gen3)
		#break
		gen1, gen2, gen3 =list(filter(lambda a: a != 'NA', gen1)),list(filter(lambda a: a != 'NA', gen2)), list(filter(lambda a: a != 'NA', gen3))
		#print (len(gen3))
		#break
		if gen2.count('0.5')<float(len(gen2))/8 and len(gen1) != 0 and len(gen2)!=0 and len(gen3)!=0:
			freq1, freq2, freq3 =sum([float(x) for x in gen1])/len(gen1), sum([float(x) for x in gen2])/len(gen2), sum([float(x) for x in gen3])/len(gen3)
			eff=line.strip().split()[-1].split('|')
		#print (eff[1], eff[2])
		#break
			if eff[2]!='MODIFIER' and freq1!=freq2 and freq1!=freq3 and freq2!=freq3:
				out.write('%s,%s,%s,%s,%s\n'%(freq1, freq2, freq3, eff[1], eff[2]))
		else:
			pass
out.close()
