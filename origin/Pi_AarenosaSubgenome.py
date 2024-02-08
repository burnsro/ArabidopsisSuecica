from optparse import OptionParser
from itertools import combinations
# pipeline params
parser = OptionParser()
parser.add_option("-v", "--vcf", dest="vcf", help="combined vcf file", default="")
parser.add_option("-d", "--dir", dest="dir", help="directory of combined vcf file", default="")

(options, args) = parser.parse_args()

acc=[]
comb=[]

out=open('%s/%s.Pi'%(options.dir, options.vcf), 'w')
for line in open('%s/%s'%(options.dir, options.vcf), 'r'):
	if '##' in line:
		pass
	else:
		if '#CHROM' in line:
			pass
			namesAS=line.strip().split()[30:44] + line.strip().split()[53:55]
			namesAA=line.strip().split()[9:30] + line.strip().split()[44:53] + line.strip().split()[55:]
			#print(namesAS)
			#print(namesAA)
			#break
			acc=namesAS+namesAA
			comb=list(combinations(acc,2))
			#print(comb)
			#break
			totalL,dist,distSh,distPrAS,distPrANC, distdiff=[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb),[0]*len(comb), [0]*len(comb)
		else:
			genAS=line.strip().split()[30:44] + line.strip().split()[53:55]
			genAA=line.strip().split()[9:30] + line.strip().split()[44:53] + line.strip().split()[55:]
			genAS,genAA=[x.split(':')[0] for x in genAS],[x.split(':')[0] for x in genAA]
			genAA_l,genAS_l=len(genAA), len(genAS)
			#print(genAA_l) 
			#print(genAS_l)
			#break
			genAAold,genASold=genAA,genAS		
			genAA,genAS=list(filter(lambda a: a != './.', genAA)),list(filter(lambda a: a != './.', genAS))
			genAA,genAS=list(filter(lambda a: a != './././.', genAA)),list(filter(lambda a: a != './././.', genAS))
			genAA=[item for sublist in [list(x[0::2] for x in genAA)] for item in sublist]
			#print(genAA)
			#break
			genAS=[item for sublist in [list(x[0::2] for x in genAS)] for item in sublist]
			genAA, genAS=list(''.join(genAA)), list(''.join(genAS))
			nonV,Sh,PrAS,PrANC, diff=0,0,0,0,0
			if len(set(genAA))==1 and len(set(genAS))==1 and set(genAA)==set(genAS):#set(genAA)==set(genAS) and len(set(genAA))==1:
				nonV=1
					##non variant
			if len(set(genAA))>1 and len(set(genAS))>1:
				Sh=1
				##shared
			if len(set(genAA))==1 and len(set(genAS))>1:
				##private
				PrAS=1
				PrANC=0
			if len(set(genAA))>1 and len(set(genAS))==1:
				PrANC=1
				PrAS=0
			if len(set(genAA))==1 and len(set(genAS))==1 and set(genAA)!=set(genAS):
				diff=1
			#print (nonV,Sh, PrAS, PrANC, diff,len(set(genAA)),len(set(genAS)))
#################################################################################
			count=0
			gen=genASold+genAAold
			for i,j in list(combinations(acc,2)):
				#print(acc)
				ind_i, ind_j=acc.index(i),acc.index(j)
				#print(ind_i)
				#print(ind_j)
				#break
				gen_i,gen_j=gen[ind_i],gen[ind_j]
				#print(gen_i)
				#print(gen_j)
				g1,g2=gen_i.split('/'), gen_j.split('/')
				#print(g1,g2)
				#print(sum(map(int,g1)), sum(map(int,g2)))
				#break
				#if len(set(list(g1)))==len(set(list(g2))):
				if './.' in gen_i or './././.' in gen_i or './.' in gen_j or './././.' in gen_j:
					pass
				elif len(set(list(g1)))==1 and len(set(list(g2)))==1 and sum(map(int,g1))==sum(map(int,g2)):
					totalL[count]+=1
				else:
					between=[(a,b) for a in g1 for b in g2]
					P=0
					for a in between:
						if len(set(a))>1:
							P+=1
					Pi=float(P)/len(between)
					dist[count]+=Pi
					if Sh==1:
						PrAS=0
						PrANC=0
						diff=0
						distSh[count]+=Pi
					elif PrAS==1:
						PrANC=0
						Sh=0
						diff=0
						distPrAS[count]+=Pi
					elif PrANC==1:
						PrAS=0
						Sh=0
						diff=0
						distPrANC[count]+=Pi
					elif diff==1 and Sh==0 and PrAS==0 and PrANC==0:
						Sh=0
						PrAS=0
						PrANC=0
						distdiff[count]+=Pi
					else:
						pass
				count+=1
out.write('acc1\t%acc2\taln\ttotal_dist\tshared_dist\tprivAS_dist\tprivANC_dist\tdistdiff\n')
for i in range(len(comb)):
	out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(list(comb[i])[0], list(comb[i])[1], totalL[i], dist[i], distSh[i], distPrAS[i], distPrANC[i], distdiff[i]))
out.close()
