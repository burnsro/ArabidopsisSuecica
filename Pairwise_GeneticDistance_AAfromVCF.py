from itertools import combinations

out=open('AAsideALL.pairwise.txt', 'w')
for line in open('AarenosaAsuecica.combined.nonvariantsbiallelicSNPs.filt.ds.vcf'):
	if '##' in line:
		pass
	else:
		if '#CHROM' in line:
			pass
			namesAS=line.strip().split()[9:23]
			namesAA=line.strip().split()[23:]
			acc=namesAS+namesAA
			comb=list(combinations(acc,2))
			totalL,dist,distSh,distPrAS,distPrANC=[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb)
		else:
			genAS=line.strip().split()[9:23]
			genAA=line.strip().split()[23:]
			genAS,genAA=[x.split(':')[0] for x in genAS],[x.split(':')[0] for x in genAA]
			genAA_l,genAS_l=len(genAA), len(genAS)
			genAAold,genASold=genAA,genAS
			genAA,genAS=list(filter(lambda a: a != './.', genAA)),list(filter(lambda a: a != './.', genAS))
			genAA,genAS=list(filter(lambda a: a != './././.', genAA)),list(filter(lambda a: a != './././.', genAS))
			if genAS.count('0/1')<1 and len(genAA)>0.5*genAA_l and len(genAS)>0.5*genAS_l:
				#genAA,genAS=[x.split('/') for x in genAA], [x.split('/') for x in genAS]
				genAA=[item for sublist in [list(x[0::2] for x in genAA)] for item in sublist]
				genAS=[item for sublist in [list(x[0::2] for x in genAS)] for item in sublist]
				genAA, genAS=list(''.join(genAA)), list(''.join(genAS))
				nonV,PrAS,PrANC,Sh=0,0,0,0

				if set(genAA)==set(genAS) and len(set(genAA))==1:
					nonV=1
					##non variant
				if len(list(set(genAA)))>1 and len(list(set(genAS)))>1:
					Sh=1
				##shared
				if len(list(set(genAA)))==1 and len(list(set(genAS)))>1:
				##private
					PrAS=1
				if len(list(set(genAA)))>1 and len(list(set(genAS)))==1:
					PrANC=1
				#print (nonV,Sh, PrAS, PrANC,list(set(genAA)),list(set(genAS)))
#################################################################################
				count=0
				gen=genASold+genAAold
				for i,j in list(combinations(acc,2)):
					ind_i, ind_j=acc.index(i),acc.index(j)
					gen_i,gen_j=gen[ind_i],gen[ind_j]
					if gen_i=='./.' or gen_j=='./.' or gen_i=='./././.' or gen_j=='./././.':
						pass
					elif gen_i==gen_j:
						totalL[count]+=1
					else:
						totalL[count]+=1
						g1,g2=gen_i.split('/'), gen_j.split('/')
						between=[(a,b) for a in g1 for b in g2]
						P=0
						for a in between:
							if len(set(list(a)))>1:
								P+=1
						Pi=float(P)/len(between)
						dist[count]+=Pi
						if Sh==1:
							distSh[count]+=Pi
						if PrAS==1:
							distPrAS[count]+=Pi
						if PrANC==1:
							distPrANC[count]+=Pi
					count+=1
out.write('acc1\t%acc2\taln\ttotal_dist\tshared_dist\tprivAS_dist\tprivANC_dist\n')
for i in range(1, len(comb)):
	out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(list(comb[i])[0], list(comb[i])[1], totalL[i], dist[i], distSh[i], distPrAS[i], distPrANC[i]))
out.close()
