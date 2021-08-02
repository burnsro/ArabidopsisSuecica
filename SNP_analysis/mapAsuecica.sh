#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=8  # 14 physical cores per task
#SBATCH --mem=74GB   # 64GB of RAM
#SBATCH --qos=short
#SBATCH --time=0-06:00:00
#SBATCH --output=%A_%a.Asuecica.stdout
#SBATCH --array=1-15

ml bwa/0.7.17-foss-2018b
ml samtools/1.9-foss-2018b

out='/scratch-cbe/users/robin.burns/003Asuecica'
ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/release/Asue_genome_210620.full.fasta'

raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/Asuecica'

samples=$raw'/samples_ASdownsampled.txt'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
acc=`echo $line | cut -d ' ' -f1`
echo $acc
cd $out


CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH/


samtools faidx $ref
bwa index $ref
bwa mem -t 12 -M -U 15 -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $ref $raw/${acc}.1.fastq.gz $raw/${acc}.2.fastq.gz > $out/$acc.sam
#increased penalty for unpaired reads
samtools view -@ 12  -bh -t $ref.fai -o $out/$acc.bam $out/$acc.sam
samtools sort -@ 12 -o $out/$acc.sort.bam $out/$acc.bam
samtools index $out/$acc.sort.bam
samtools rmdup $out/$acc.sort.bam $out/$acc.rmdup.bam
samtools index $out/$acc.rmdup.bam
samtools view -F 256 -q 5 -f 3 -h -b $out/$acc.rmdup.bam > $out/$acc.filt.bam
#primary aligned reads, unique and mapped in proper pair
samtools index $out/$acc.filt.bam
#final file for SNP calling

################################################################################
#################Splitting subgenomes################################
cd $out
ml bamtools/2.5.1-foss-2018b
bamtools split -in $out/$acc.filt.bam -reference

bamtools merge -in ${acc}.filt.REF_Asue_scaffold1.bam -in ${acc}.filt.REF_Asue_scaffold2.bam -in ${acc}.filt.REF_Asue_scaffold3.bam -in ${acc}.filt.REF_Asue_scaffold4.bam -in ${acc}.filt.REF_Asue_scaffold5.bam -out ${acc}.AT.bam
samtools index ${acc}.AT.bam

bamtools merge -in ${acc}.filt.REF_Asue_scaffold6.bam -in ${acc}.filt.REF_Asue_scaffold7.bam -in ${acc}.filt.REF_Asue_scaffold8.bam -in ${acc}.filt.REF_Asue_scaffold9.bam -in ${acc}.filt.REF_Asue_scaffold10.bam -in ${acc}.filt.REF_Asue_scaffold11.bam -in ${acc}.filt.REF_Asue_scaffold12.bam -in ${acc}.filt.REF_Asue_scaffold13.bam -out ${acc}.AA.bam
samtools index ${acc}.AA.bam


samtools view -h ${acc}.AT.bam > ${acc}.AT.sam
samtools view -h ${acc}.AA.bam > ${acc}.AA.sam
samtools -@ 8 sort -n ${acc}.AT.bam -o ${acc}.AT.sort.bam
samtools -@ 8 sort -n ${acc}.AA.bam -o ${acc}.AA.sort.bam
samtools view -h ${acc}.AT.sort.bam > ${acc}.AT.sort.sam
samtools view -h ${acc}.AA.sort.bam > ${acc}.AA.sort.sam
script='/groups/nordborg/projects/suecica/005scripts/002scripts201718/my_scripts'
perl $script/samToFastq.pl $out/$acc.AT.sam
perl $script/samToFastq.pl $out/$acc.AA.sam
ml  bioperl/1.7.2-foss-2018b-perl-5.28.0
perl $script/split1to2.pl $out/$acc.AT.sort.sam.fastq $out/$acc.AT.1.fastq $out/$acc.AT.2.fastq 101
perl $script/split1to2.pl $out/$acc.AA.sort.sam.fastq $out/$acc.AA.1.fastq $out/$acc.AA.2.fastq 101

###################################################################################################
