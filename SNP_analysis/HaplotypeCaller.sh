#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=4  # 14 physical cores per task
#SBATCH --mem=34G   # 64GB of RAM
#SBATCH --qos=short
#SBATCH --time=0-02:40:00
#SBATCH --output=%A_%a.suesnps.stdout
# SBATCH --array=1-13

ml gatk/3.8-1-java-1.8
ml vcftools/0.1.16-foss-2018b-perl-5.28.0

out='/scratch-cbe/users/robin.burns/021MapAs_2020/parents'
ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/release/Asue_genome_210620.full.fasta'



#ref='/scratch-cbe/users/robin.burns/014Fullgenome/Asue_genome_210620.full.fasta'
#out='/scratch-cbe/users/robin.burns/014Fullgenome/map_bam'
#samples=$out'/samples'
#export l=$SLURM_ARRAY_TASK_ID\p
#line=`sed -n $l $samples`
#acc=`echo $line | cut -d ' ' -f1`
#echo $acc
#ploidy=`echo $line | cut -d ' ' -f2`
#echo $ploidy

samples=$out'/chromosomes.txt'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
chrom=`echo $line | cut -d ' ' -f1`
echo $chrom


TMP='/scratch-cbe/users/robin.burns/tmp'
cd $out
#==================================== STEP1 =========================================
# ==== call SNPs with GATK - Haplotype Caller per each Chromosome per individual ====
#echo "for the given input file ${inputfile}, the output vcf file is $output_vcf_file"
#chrom="Asue_scaffold12"
#ploidy=2
#mv ${out}/${acc}.intervals.${chrom} ${out}/${acc}.${chrom}.intervals
#acc='ASS3a'
#java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R ${ref} \
#  -T HaplotypeCaller \
#   -I ${out}/${acc}.filt.bam\
#   --emitRefConfidence GVCF \
#   -o ${out}/${acc}.${chrom}.g.vcf \
#  --sample_ploidy ${ploidy} \
#   -variant_index_type LINEAR \
#   -variant_index_parameter 128000 \
#   -nct 8 \
#  -L ${chrom}
# 
#==================================== STEP2 ===============================================
# ==== combine g.vcf files of all individuals for a given chromosome with CombineGVCFs ====
#input files are *.g.vcf files obtained from the script above
#cd $out
#ls ASS3*$chrom.g.vcf > ASS3.$chrom.list
#ls AS150*$chrom.g.vcf > AS150.$chrom.list
#
#cat ASS3.$chrom.list AS150.$chrom.list > parents.$chrom.list
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R ${ref} \
#   -T CombineGVCFs \
#   -V ${out}/parents.${chrom}.list \
#   -o ${out}/ASS3AS150.${chrom}.combined.g.vcf 

#=========================================== STEP3 =======================================================
# ==== convert g.vcf file of a given chr to a vcf with GenotypeGVCFs for joint genotyping of all inds ====
#input_file2=${out}/Phylo2AT.combined.g.vcf
#output_file2=${out}/Phylo2AT.combined.vcf
##echo "for the input file $input_file2, the output file is $output_file2"
#java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R $ref \
#   -T GenotypeGVCFs \
#   --variant ${out}/ASS3AS150.${chrom}.combined.g.vcf \
#   -o ${out}/ASS3AS150.${chrom}.combined.vcf  \
#   --includeNonVariantSites \
#   -nt 6

#================================================ STEP4 ==================================================
# ==== merge vcf files from all chrs to produce a single vcf prior to SNP filtration with CatVariants ====
#CatVariants concatenates vcf files of non-overlapping genome intervals, with the same set and ORDER of samples
#java -Xmx30g -Djava.io.tmpdir=$TMP -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#   -R $ref \
#   --variant ${out}/ASS3AS150.Asue_scaffold1.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold2.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold3.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold4.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold5.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold6.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold7.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold8.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold9.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold10.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold11.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold12.combined.vcf \
#   --variant ${out}/ASS3AS150.Asue_scaffold13.combined.vcf \
#   -out ${out}/ASS3AS150.Asue_scaffoldAll.combined.vcf \
#   -assumeSorted \
#   -log cat.all.chrs.log \
#   --logging_level ERROR \
#
#java -Xmx80g -Djava.io.tmpdir=$TMP -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#   -R $ref \
#   --variant ${out}/Asuecica.Asue_scaffold1.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold2.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold3.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold4.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold5.combined.vcf \
#   -out ${out}/Asuecica.Asue_scaffoldAT.combined.vcf \
#   -assumeSorted \
#   -log cat.all.chrs.log \
#   --logging_level ERROR \


#java -Xmx80g -Djava.io.tmpdir=$TMP -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#	-R $ref \
#	--variant ${out}/ASS3AS150.Asue_scaffold1.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold2.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold3.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold4.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold5.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold6.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold7.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold8.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold9.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold10.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold11.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold12.combined.vcf \
#	--variant ${out}/ASS3AS150.Asue_scaffold13.combined.vcf \
#	-out ${out}/Parents.ASS3aAS150.Asue_scaffolds.combined.vcf \
#	-assumeSorted \
#	   -log cat.all.chrs.log \
#	   --logging_level ERROR \
### -assumeSorted flag requires input variant files to be ordered according to the appearance of intervals in the ref genome
### Chr_01 < Chr_02 < Chr_03 < .... < Chr_07
### like simple cat command, it appends the nth file to the end of n-1 th file


##################### NOW THE FILTERING STEP ###################################

# =============== CALL ONLY NON-VARIANT SITES ===============
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#  -R $ref \
#   -T SelectVariants \
#   --variant ${out}/Parents.ASS3aAS150.Asue_scaffolds.combined.vcf \
#   -o $out/Parents.ASS3aAS150.Asue_scaffolds.combined.nonvariants.vcf \
#   -selectType NO_VARIATION \
#   -nt 16

# ============== CALL BIALLELIC VARIANT SITES ===============
# in addition to excluding non-variant sites, restrict variant sites to biallelic SNPs
#java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R $ref \
#   -T SelectVariants \
#   --variant ${out}/ASS3AS150.Asue_scaffoldAll.combined.vcf \
#   -o $out/ASS3AS150.Asue_scaffoldAll.combined.biallelicSNPS.vcf \
#   -selectType SNP \
#   --restrictAllelesTo BIALLELIC \
#   -nt 4

# ============ COMBINE NON PLUS BIVARIANT SITES ============
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R ${ref} \
#   -T CombineVariants \
#   --variant $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.vcf \
#   --variant $out/AarenosaChip.Asue_scaffoldAll.nonvariants.biallelic.combined.vcf \
#   -o $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.AAgeneINT.vcf \
#   -L /groups/nordborg/projects/suecica/mygenePOS_AAgenesAsuecica.intervals \
#   -genotypeMergeOptions UNIQUIFY \
#   --disable_auto_index_creation_and_locking_when_reading_rods \
#   -nt 16

#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#	  -R $ref \
#	  -T SelectVariants \
#	  --variant $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.AAgeneINT.vcf \
#	  -o $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.AAgeneINT.newsamples.vcf \
#          -selectType NO_VARIATION \
#           -L /groups/nordborg/projects/suecica/mygenePOS_AAgenesAsuecica.intervals \
#	   -nt 16

#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#	      -R $ref \
#	      -T SelectVariants \
#	      --variant $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.AAgeneINT.vcf \
#             -o $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.biallelic.AAgeneINT.newsamples.vcf \
#	      -selectType SNP \
#	      --restrictAllelesTo BIALLELIC \
#	      -L /groups/nordborg/projects/suecica/mygenePOS_AAgenesAsuecica.intervals \
#	      -nt 16

#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#	   -R ${ref} \
#	   -T CombineVariants \
#	   --variant $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.biallelic.AAgeneINT.newsamples.vcf \
#	   --variant $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.AAgeneINT.newsamples.vcf \
#	   -o $out/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariantsbiallelicSNPs.2020.AAgeneINT.newsamples.vcf \
#	   --assumeIdenticalSamples \
#	   --disable_auto_index_creation_and_locking_when_reading_rods \
#	   -L /groups/nordborg/projects/suecica/mygenePOS_AAgenesAsuecica.intervals \
#	   -nt 16

cd $out
myvcf='ASS3AS150.Asue_scaffoldAll.combined.biallelicSNPS.vcf'
CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
ml anaconda3/2019.03
ml anaconda3/2019.03
source activate RobinCondaSCRATCH/  
cd $out
ml java/1.8.0_212
#export PATH='/groups/nordborg/projects/suecica/005scripts/001Software/snpEff':$PATH
#snpEff='/groups/nordborg/projects/suecica/005scripts/001Software/snpEff'

#java -Djava.io.tmpdir=$TMP -Xmx60g -jar ${snpEff}/snpEff.jar ann -c ${snpEff}/snpEff.config -v asuecica $out/$myvcf -ud 1000 > $out/$myvcf.ann

vcftools --vcf ${myvcf} --max-missing 0 --mac 2 --minQ 60 --recode --recode-INFO-all --out ${myvcf}_filtermismq60mac1

vcftools --vcf ${myvcf}_filtermismq60mac1.recode.vcf --minDP 20 --recode --recode-INFO-all --out ${myvcf}_filtermismq60mac1dp5

vcftools --vcf ${myvcf}_filtermismq60mac1dp5.recode.vcf --max-meanDP 80 --recode --recode-INFO-all --out ${myvcf}_filtermismq60mac1dp580

#vcftools --vcf Asuecica.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.filtermismq60mac3dp580.recode.vcf --keep cross.txt --recode --recode-INFO-all --out Asuecica.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.filtermismq60mac3dp580_Cross

#awk -F '\t' '$5!="\*" {print $0}' ${myvcf}_filtermismq60mac1dp580.recode.vcf > ${myvcf}_filtermismq60mac1dp580.recode.noasterisk.vcf 

#vcftools --vcf ${myvcf}_filtermismq60mac1dp580.recode.noasterisk.vcf --maf 0.5 --recode --recode-INFO-all --out ${myvcf}_filtermismq60mac1dp580noastmafhalf

#cat Parents.ASS3aAS150.Asue_scaffolds.combined.nonvariants.biallelicSNPs.vcf_filtermismq60mac1dp580noastmafhalf.recode.vcf | grep -v "##" | awk -F '\t' '{print $1FS$2FS$4FS$5}' > Parents.ASS3aAS150.Asue_scaffolds.1119.Good_positions.txt

#java -Xmx60G -Djava.io.tmpdir=$LOCAL_DISK -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R $ref \
#   -T VariantsToTable \
#   --variant $out/Asuecica.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.filtermismq60mac3dp580_Cross.recode.vcf \
#   -o $out/Asuecica.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.filtermismq60mac3dp580_Cross.table.txt \
#   -F CHROM -F POS -F REF -F ALT -F FILTER -F NO-CALL \
#   --genotypeFields GT \
#   --allowMissingData

#scripts='/groups/nordborg/projects/suecica/005scripts/old_cluster/Popgen'
#myvcf='/scratch-cbe/users/robin.burns/003SNPs/AsuecicaPhyloAthaliana.Asue_scaffoldAll.combined.nonvariants.biallelicSNPs.vcf'
#cd /scratch-cbe/users/robin.burns/003SNPs/
#CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
#CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
#cd /groups/nordborg/projects/suecica/005scripts/001Software
#ml anaconda3/2019.03
#source activate RobinCondaSCRATCH/
#cd /scratch-cbe/users/robin.burns/003SNPs/
#python $scripts/CleanVCF.py --vcf ${myvcf}


#perl -i -p -e 's/:\S*//g' ${myvcf}_clean
#python /home/GMI/robin.burns/scripts/TE/replacetet.py -v Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_clean
#python /home/GMI/robin.burns/scripts/TE/replacedip.py -v Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_clean.AlleleCount

#python VCF_Athpart.py --vcf Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_clean.AlleleCount.AlleleCount --dir ${out}
#python VCF_Aarpart.py --vcf Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_clean.AlleleCount.AlleleCount --dir ${out}

#awk -F '\t' '$1~"#" || $1=="Asue_scaffold1" || $1=="Asue_scaffold2" || $1=="Asue_scaffold3" || $1=="Asue_scaffold4" || $1=="Asue_scaffold5" {print $0}' ${myvcf}_clean > ${myvcf}_clean_Athpart

#awk -F '\t' '$1~"#" || $1=="Asue_scaffold6" || $1=="Asue_scaffold7" || $1=="Asue_scaffold8" || $1=="Asue_scaffold9" || $1=="Asue_scaffold10" || $1=="Asue_scaffold11" || $1=="Asue_scaffold12" || $1=="Asue_scaffold13" {print $0}' ${myvcf}_clean > ${myvcf}_clean_Aarpart

#python selectcol.py
#python filtvcf2.py
#python CombinedVCFtodiploid.py --vcf Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_Athpart_dsfilt --dir ${out}
#cd $out
#python VCFstats.py --vcf Phylo2AS.Asue_scaffoldall.combined.nonvariants.biallelic.SNPs.vcf_Aarpart_dsfilt.split.vcf --dir $out > test.txt
