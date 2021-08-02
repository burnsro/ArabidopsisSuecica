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



