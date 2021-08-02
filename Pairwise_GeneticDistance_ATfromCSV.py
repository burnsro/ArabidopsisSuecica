from itertools import combinations

out=open('ATside_forpi.geneticdistance', 'w')
for line in open('ATside_forpi.csv'):  #ATsideV3.csv
	line=line.strip().split(',')
	if 'chr' in line:
		header=line
		acc=line[2:]
		comb=list(combinations(acc,2))
		totalL,dist,distSh,distPrAS,distPrANC=[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb)
	else:
		gen=line[2:]
		genANC,genAS=gen[0:30], gen[30:44]+gen[96:97]
		genANC_l,genAS_l=30, 15		
		genANC,genAS=list(filter(lambda a: a != 'NA', genANC)),list(filter(lambda a: a != 'NA', genAS))	
		if gen.count('0.5')<20:
			nonV,PrAS,PrANC,Sh=0,0,0,0	
			if set(genANC)==set(genAS) and len(set(genANC))==1:
				nonV=1	
			if len(list(set(genANC)))>1 and len(list(set(genAS)))>1:	
				Sh=1	
			if len(list(set(genANC)))==1 and len(list(set(genAS)))>1:
				PrAS=1
			if len(list(set(genANC)))>1 and len(list(set(genAS)))==1:
				PrANC=1
#################################################################################
			count=0
			for i,j in list(combinations(acc,2)):
				ind_i, ind_j=acc.index(i),acc.index(j)
				gen_i,gen_j=gen[ind_i],gen[ind_j]
				if gen_i=='NA' or gen_j=='NA':
					pass
				elif gen_i==gen_j:
					totalL[count]+=1
				else:
					dist[count]+=abs(float(gen_i)-float(gen_j))
					if Sh==1:
						distSh[count]+=abs(float(gen_i)-float(gen_j))
					if PrAS==1:
						distPrAS[count]+=abs(float(gen_i)-float(gen_j))
					if PrANC==1:
						distPrANC[count]+=abs(float(gen_i)-float(gen_j))
				count+=1	
out.write('acc1\t%acc2\taln\ttotal_dist\tshared_dist\tprivAS_dist\tprivANC_dist\n')
for i in range(1, len(comb)):
	out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(list(comb[i])[0], list(comb[i])[1], totalL[i], dist[i], distSh[i], distPrAS[i], distPrANC[i]))
out.close()
