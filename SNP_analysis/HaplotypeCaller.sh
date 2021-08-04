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

#First call snps using individuals and ploidy then switch to per chromosome
#Go through step by step

samples=$out'/samples'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
acc=`echo $line | cut -d ' ' -f1`
echo $acc
ploidy=`echo $line | cut -d ' ' -f2`
echo $ploidy

#Switch to chromosomes after step1
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
chrom="Asue_scaffold1" #change
ploidy=2
mv ${out}/${acc}.intervals.${chrom} ${out}/${acc}.${chrom}.intervals
acc='ASS3a'
java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
   -R ${ref} \
  -T HaplotypeCaller \
   -I ${out}/${acc}.filt.bam\
   --emitRefConfidence GVCF \
   -o ${out}/${acc}.${chrom}.g.vcf \
  --sample_ploidy ${ploidy} \
   -variant_index_type LINEAR \
   -variant_index_parameter 128000 \
   -nct 8 \
  -L ${chrom}
# 
#==================================== STEP2 ===============================================
# ==== combine g.vcf files of all individuals for a given chromosome with CombineGVCFs ====
#input files are *.g.vcf files obtained from the script above
#cd $out
#ls *$chrom.g.vcf > $chrom.list


#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R ${ref} \
#   -T CombineGVCFs \
#   -V ${out}/${chrom}.list \
#   -o ${out}/Asuecica.${chrom}.combined.g.vcf 

#=========================================== STEP3 =======================================================
# ==== convert g.vcf file of a given chr to a vcf with GenotypeGVCFs for joint genotyping of all inds ====
#input_file2=${out}/Phylo2AT.combined.g.vcf
#output_file2=${out}/Phylo2AT.combined.vcf
##echo "for the input file $input_file2, the output file is $output_file2"
#java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R $ref \
#   -T GenotypeGVCFs \
#   --variant ${out}/Asuecica.${chrom}.combined.g.vcf \
#   -o ${out}/Asuecica.${chrom}.combined.vcf  \
#   --includeNonVariantSites \
#   -nt 6

#================================================ STEP4 ==================================================
# ==== merge vcf files from all chrs to produce a single vcf prior to SNP filtration with CatVariants ====
#CatVariants concatenates vcf files of non-overlapping genome intervals, with the same set and ORDER of samples
#java -Xmx30g -Djava.io.tmpdir=$TMP -cp $EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
#   -R $ref \
#   --variant ${out}/Asuecica.Asue_scaffold1.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold2.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold3.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold4.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold5.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold6.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold7.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold8.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold9.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold10.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold11.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold12.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold13.combined.vcf \
#   -out ${out}/Asuecica.Asue_scaffoldAll.combined.vcf \
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
#   -R $ref \
#   --variant ${out}/Asuecica.Asue_scaffold6.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold7.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold8.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold9.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold10.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold11.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold12.combined.vcf \
#   --variant ${out}/Asuecica.Asue_scaffold13.combined.vcf \
#   -out ${out}/Asuecica.Asue_scaffoldAA.combined.vcf \
#   -assumeSorted \
#   -log cat.all.chrs.log \
#   --logging_level ERROR \


##################### NOW THE FILTERING STEP ###################################

# =============== CALL ONLY NON-VARIANT SITES ===============
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#  -R $ref \
#   -T SelectVariants \
#   --variant ${out}/Asuecica.Asue_scaffolds.combined.vcf \
#   -o $out/Asuecica.Asue_scaffolds.combined.nonvariants.vcf \
#   -selectType NO_VARIATION \
#   -nt 16

# ============== CALL BIALLELIC VARIANT SITES ===============
# in addition to excluding non-variant sites, restrict variant sites to biallelic SNPs
#java -Xmx30g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R $ref \
#   -T SelectVariants \
#   --variant ${out}/Asuecica.Asue_scaffoldAll.combined.vcf \
#   -o $out/Asuecica.Asue_scaffoldAll.combined.biallelicSNPS.vcf \
#   -selectType SNP \
#   --restrictAllelesTo BIALLELIC \
#   -nt 4

# ============ COMBINE NON PLUS BIVARIANT SITES ============
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#   -R ${ref} \
#   -T CombineVariants \
#   --variant $out/Asuecica.Asue_scaffoldAll.combined.biallelicSNPS.vcf \
#   --variant $out/Asuecica.Asue_scaffolds.combined.nonvariants.vcf \
#   -o $out/Asuecica.Asue_scaffolds.combined.nonvariants.biallelicSNPs.vcf \
#  --assumeIdenticalSamples \
#   --disable_auto_index_creation_and_locking_when_reading_rods \
#   -nt 16

#Get subgenomes from VCF file
#vcf='As_SNPs2020.Asue_scaffoldAll.nonVars.BiSNPS.vcf'
#awk -F '\t' '$1~"#" || $1=="Asue_scaffold1" || $1=="Asue_scaffold2" || $1=="Asue_scaffold3" || $1=="Asue_scaffold4" || $1=="Asue_scaffold5" {print $0}' $raw/${vcf} > $out/${vcf}_Ath

#awk -F '\t' '$1~"#" || $1=="Asue_scaffold6" || $1=="Asue_scaffold7" || $1=="Asue_scaffold8" || $1=="Asue_scaffold9" || $1=="Asue_scaffold10" || $1=="Asue_scaffold11" || $1=="Asue_scaffold12" || $1=="Asue_scaffold13" {print $0}' ${raw}/${vcf} > ${out}/${vcf}_Aar

#downsample columns to remove A. thaliana acc in the A. arenosa subgenome VCF and A. arenosa acc in the A. thaliana VCF
#because everything was mapped altogether to reduce bias in mapping
#downsample Alyrata to just european A. lyrata petreae
#cd $out
#python selectcol_AA.py
#python selectcol_AT.py
#python selectcol_AL.py


#to polarize map Athaliana to A. arenosa subgenome only and call SNPs
#map European Alyrata.lyratapetreae (more individuals than A.arenosa) to A. thaliana and call SNPs

#Combine vcfs with Al2At and At2AA


#refdir='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/release'
#java -Xmx80g -Djava.io.tmpdir=/scratch-cbe/users/robin.burns/tmp -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#                      -R ${refdir}/Asue_genome_210620.full.Popte2.AA.fasta \
#                      -T CombineVariants \
#                      --variant $out/As_SNPs2020.Asue_scaffoldAll.nonVars.BiSNPS.DS.vcf_Aar \
#                      --variant ${out}/At_SNPS2020.AsAA.Asue_scaffoldAA.combined.NonVarBiSNPs.vcf \
#                      -o $out/As_SNPs2020.AsAA.forpolarizing.nonVars.BiSNPs.vcf \
#                      -nt 16 \



#java -Xmx80g -Djava.io.tmpdir=/scratch-cbe/users/robin.burns/tmp -jar $EBROOTGATK/GenomeAnalysisTK.jar \
#                       -R ${refdir}/Asue_genome_210620.full.Popte2.AT.fasta \
#                      -T CombineVariants \
#                      --variant $out/As_SNPs2020.Asue_scaffoldAll.nonVars.BiSNPS.DS.vcf_Ath\
#                      --variant ${out}/Al_SNPS2020.AsAT.Asue_scaffoldAT.combined.NonVarBiSNPs.DS.vcf\
#                      -o $out/As_SNPs2020.AsAT.forpolarizing.nonVars.BiSNPs.vcf\
#                      -nt 16 \

#vcfAA=As_SNPs2020.AsAA.forpolarizing.nonVars.BiSNPs.vcf
#python filtCombVCF.py -v ${vcfAA} -d ${out} #remove missing sites 


#Run SNPeff

#Now make csv
#vcfAA=As_SNPs2020.AsAA.forpolarizing.nonVars.BiSNPs.vcf.filt.ANN
#python vcf2csv_AA_ann.py -v ${vcfAA} -d ${out}

#vcfAT=As_SNPs2020.AsAT.forpolarizing.nonVars.BiSNPs.vcf.filt.ANN
#python vcf2csv_AT_ann.py -v ${vcfAT} -d ${out}

#vcfAT=As_SNPs2020.AsAT.forpolarizing.nonVars.BiSNPs.vcf.filt.ANN.csv
#python get_fixed_SNPs_positionsAthsubgenome.py -v ${vcfAT} -d ${out}

#vcfAA=As_SNPs2020.AsAA.forpolarizing.nonVars.BiSNPs.vcf.filt.ANN.csv
#python get_fixed_SNPs_positionsAarsubgenome.py -v ${vcfAA} -d ${out}

#Now polarize 
#per chromosome
samples=$out'/chromosomes.txt'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
chrom=`echo $line | cut -d ' ' -f1`
echo $chrom


python polarize.py -c ${chrom}
python polarizeAA.py -c ${chrom}
