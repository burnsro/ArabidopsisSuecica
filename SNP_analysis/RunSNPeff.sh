#!/usr/bin/env bash                                                                                 
#SBATCH --nodes=1
#SBATCH --ntasks=1  # 14 physical cores per task
#SBATCH --mem=44G   # 64GB of RAM
#SBATCH --qos=medium
#SBATCH --time=0-16:00:00
#SBATCH --output=%A_%a.SNPeffAT.stdout

snpEff='/groups/nordborg/projects/suecica/005scripts/001Software/snpEff'
out='/groups/nordborg/projects/suecica/015SNPs'

vcf='/path/to/vcf'
TMP='/scratch-cbe/users/robin.burns/tmp'
ml java/1.8.0_212

java -Djava.io.tmpdir=$TMP -Xmx40g -jar ${snpEff}/snpEff.jar ann -c ${snpEff}/snpEff.config -v asuecica $vcf -ud 500 > $vcf.ANN

######################################
#module load tabix
#bgzip -f $out/$vcf.ann
#tabix -p vcf $out/$vcf.ann.gz

