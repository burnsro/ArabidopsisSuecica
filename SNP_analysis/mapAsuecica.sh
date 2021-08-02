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
#ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/Asue_genome_210620.full.Popte2.fasta'

out='/scratch-cbe/users/robin.burns/003Asuecica'
ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/release/Asue_genome_210620.full.fasta'

#raw='/groups/nordborg/projects/nordborg_rawdata/Phylogenomics/phylogenomics/RAWdata/Phylogenomics/renamedRAW'
raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/Asuecica'

#raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/AthAnc'
#raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/Artificial_suecica'
#raw='/groups/nordborg/projects/suecica/012ChIP/001fastq'
#raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/Artificial_suecica/Asa1_generations'

samples=$raw'/samples_ASdownsampled.txt'
#samples=$raw'/samples_parents.txt'
#samples=$raw'/samples_arenosa.txt'
#samples=$raw'/Aa_input.txt'
#samples=$raw'/samples_AT30.txt'
#samples=$raw'/samples_AaAt'
#samples=$raw'/diploids_lyrata.txt'
#samples=$out'/mysamples'
export l=$SLURM_ARRAY_TASK_ID\p
line=`sed -n $l $samples`
acc=`echo $line | cut -d ' ' -f1`
#acc='ASS3a'
echo $acc
#ploidy=`echo $line | cut -d ' ' -f2`
#echo $ploidy
cd $out
#acc='Aarenosaarenosa4'
CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH/

#cd /scratch-cbe/users/robin.burns/021MapAs_2020/coverage_plots
#samtools faidx $ref
#bwa index $ref
#bwa mem -t 12 -M -U 15 -R '@RG\tID:'$acc'\tSM:'$acc'\tPL:Illumina\tLB:'$acc $ref $raw/${acc}.1.fastq.gz $raw/${acc}.2.fastq.gz > $out/$acc.sam
#samtools view -@ 12  -bh -t $ref.fai -o $out/$acc.bam $out/$acc.sam
#samtools sort -@ 12 -o $out/$acc.sort.bam $out/$acc.bam
#samtools index $out/$acc.sort.bam
#samtools rmdup $out/$acc.sort.bam $out/$acc.rmdup.bam
#samtools index $out/$acc.rmdup.bam
#samtools view -F 256 -q 5 -f 3 -h -b $out/$acc.rmdup.bam > $out/$acc.filt.bam  
#samtools index $out/$acc.filt.bam

#rm $out/$acc.bam $out/$acc.sam $out/$acc.sort.bam $out/$acc.rmdup.bam $out/$acc.sort.bam.bai $out/$acc.rmdup.bam.bai
cd $out
#ml bamtools/2.5.1-foss-2018b
#bamtools split -in $out/$acc.filt.bam -reference

#bamtools merge -in ${acc}.filt.REF_Asue_scaffold1.bam -in ${acc}.filt.REF_Asue_scaffold2.bam -in ${acc}.filt.REF_Asue_scaffold3.bam -in ${acc}.filt.REF_Asue_scaffold4.bam -in ${acc}.filt.REF_Asue_scaffold5.bam -out ${acc}.AT.bam
#samtools index ${acc}.AT.bam

#bamtools merge -in ${acc}.filt.REF_Asue_scaffold6.bam -in ${acc}.filt.REF_Asue_scaffold7.bam -in ${acc}.filt.REF_Asue_scaffold8.bam -in ${acc}.filt.REF_Asue_scaffold9.bam -in ${acc}.filt.REF_Asue_scaffold10.bam -in ${acc}.filt.REF_Asue_scaffold11.bam -in ${acc}.filt.REF_Asue_scaffold12.bam -in ${acc}.filt.REF_Asue_scaffold13.bam -out ${acc}.AA.bam
#samtools index ${acc}.AA.bam


#samtools view -h ${acc}.AT.bam > ${acc}.AT.sam
#samtools view -h ${acc}.AA.bam > ${acc}.AA.sam
#samtools -@ 8 sort -n ${acc}.AT.bam -o ${acc}.AT.sort.bam
#samtools -@ 8 sort -n ${acc}.AA.bam -o ${acc}.AA.sort.bam
#samtools view -h ${acc}.AT.sort.bam > ${acc}.AT.sort.sam
#samtools view -h ${acc}.AA.sort.bam > ${acc}.AA.sort.sam
script='/groups/nordborg/projects/suecica/005scripts/002scripts201718/my_scripts'
perl $script/samToFastq.pl $out/$acc.AT.sam
perl $script/samToFastq.pl $out/$acc.AA.sam
ml  bioperl/1.7.2-foss-2018b-perl-5.28.0
perl $script/split1to2.pl $out/$acc.AT.sort.sam.fastq $out/$acc.AT.1.fastq $out/$acc.AT.2.fastq 101
perl $script/split1to2.pl $out/$acc.AA.sort.sam.fastq $out/$acc.AA.1.fastq $out/$acc.AA.2.fastq 101


#rm ${acc}*filt.REF_Asue_scaffold*.bam



#mv $out/$acc.filt.bam $out2/$acc.filt.bam
#mv $out/$acc.rmdup.bam $out2/$acc.rmdup.bam
#mv $out/$acc.sort.bam $out2/$acc.sort.bam
#mv $out/$acc.sort.bam.bai $out2/$acc.filt.bam.bai
#ml r/3.5.1-foss-2018b
#R --vanilla --slave --no-save --no-restore --args $out/${acc}.filt.bam < /groups/nordborg/projects/suecica/005scripts/get_rDNA_coverage.R 

#Popte2
#TMP='/scratch-cbe/users/robin.burns/tmp'
#hier='/groups/nordborg/projects/suecica/002Annotation/002Asuecica/002TEs/AS_RM_TEs_hier_131119.hier'
#ml java/1.8.0_212
#cd $out
#popte2='/groups/nordborg/projects/suecica/005scripts/001Software/popte2-v1.10.04.jar'
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} ppileup --bam ${acc}.rmdup.bam --map-qual 15 --hier ${hier} --output ${acc}.ppileup.gz
# identify TE signatures
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} identifySignatures --ppileup ${acc}.ppileup.gz --mode separate --output ${acc}.signatures --min-count 2
# calculate frequencies
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} frequency --ppileup ${acc}.ppileup.gz --signature ${acc}.signatures --output ${acc}.freqsig
# write TE output 
#java -Xmx60g -Djava.io.tmpdir=$TMP -jar ${popte2} pairupSignatures --signature ${acc}.freqsig --ref-genome $ref --hier ${hier} --min-distance -200 --max-distance 500 --output ${acc}.teinsertions  

#Lost genes
#ml bedtools/2.25.0-foss-2018b
#ml bedtools/2.27.1-foss-2018b
#atgenes='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/at_lostgenesAS.gff'
#algenes='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/at_lostgenesAS.gff'
#together='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/atal_lostgenesAS.gff'
#chr1='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/atal_lostgenesAS.chr1.gff'
#justalgenes='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/lsA_ALgene.bed'
#atall='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/AT.clean2.bed'
#alall='/groups/nordborg/projects/suecica/002Annotation/001_2016_ATAL/AL.clean2.bed'
#cd /scratch-cbe/users/robin.burns/004lostgenes
#raw='/groups/nordborg/projects/suecica/016_SuecicaAnalysis_20162017/map/bam'
#samples=$raw'/mysamples'
#raw='/groups/nordborg/projects/nordborg_rawdata/Phylogenomics/phylogenomics/PhylogenomicsMap/map_out_AL/bam'
#samples='/groups/nordborg/projects/nordborg_rawdata/Phylogenomics/phylogenomics/RAWdata/Phylogenomics/renamedRAW/samples_arenosa.txt'
#raw='/groups/nordborg/projects/nordborg_rawdata/Asuecica/short_reads/A_suecica/RAWdata/AthAnc/'
#samples=$raw/'New_indAthAnc.txt'
#export l=$SLURM_ARRAY_TASK_ID\p
#line=`sed -n $l $samples`
#acc=`echo $line | cut -d ' ' -f1`
#echo $acc
#acc='AS150'
#ml bedtools/2.25.0-foss-2018b
#cd /scratch-cbe/users/robin.burns/004lostgenes
#bedtools coverage -b ${raw}/${acc}.filt.REF_Chr1.bam -a $chr1 > ${acc}.lostgenes.bedtoolscov_1 
#cd ${raw}
#samtools index ${raw}/${acc}.filt.bam
#cd /scratch-cbe/users/robin.burns/004lostgenes
#acc='Alyratalyrata2'
#samtools index ${raw}/${acc}.realigned.sort.bam
#samtools bedcov ${alall} ${raw}/${acc}.filt.bam  -Q 30 -j > ${acc}.alall.samtoolscov 
#cd /groups/nordborg/projects/transposons/Robin/centmapped
#CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
#CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
#ml anaconda3/2019.03
#conda create -p /groups/nordborg/projects/suecica/005scripts/001Software/RobinCondaSCRATCH python=2.7 #or 3.6 just run this command once then source it
#cd /groups/nordborg/projects/suecica/005scripts/001Software
#ml anaconda3/2019.03
#source activate RobinCondaSCRATCH/
#cd /groups/nordborg/projects/transposons/Robin/centmapped
#ml samtools/0.1.20-foss-2018b
#ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/Asue_genome.HiCGeneticMap.270919.fasta'
#samtools mpileup -f ${ref} $acc.filt.bam > $acc.mpileup.txt


#ml bedtools/2.25.0-foss-2018b
#bedtools genomecov -ibam $acc.filt2.bam -bga -split | awk -F '\t' '$4==0 {print $0}' > ${acc}.zerocov.bed

#samtools index ${acc}.filt.bam
#cd ${raw}
#cd /groups/nordborg/projects/transposons/Robin/centmapped
#cd /groups/nordborg/projects/nordborg_rawdata/Phylogenomics/phylogenomics/PhylogenomicsMap/map_out_AT/filtBAM
#acc='ASS3a'
#ml bamtools/2.5.1-foss-2018b
#bamtools split -in ${acc}.filt.bam -reference
#bamtools split -in 14318.filt.bam -reference




#Crossmap
#ml bbmap/38.26-foss-2018b
#pileup.sh in=$out/${acc}.filt.bam out=$out/${acc}.reads.coverage.txt overwrite=true


#cp /groups/nordborg/projects/suecica/005scripts/get_percentageRNAreads_crossmapped*.R $out
#cd $out
#R --vanilla --slave --no-save --no-restore --args $out/${acc}.reads.coverage.txt < get_percentageRNAreads_crossmapped.R 

#R --vanilla --slave --no-save --no-restore --args $out/${acc}_eagle/${acc}.reads.coverage.txt < get_percentageRNAreads_crossmapped_AT.R 

CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'
#ml anaconda3/2019.03
#conda create -p /groups/nordborg/projects/suecica/005scripts/001Software/RobinCondaSCRATCH python=2.7 #or 3.6 just run this command once then source it
cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH/

#samtools index $out/$acc.filt.bam
#ml samtools/0.1.20-foss-2018b
#samtools mpileup -f ${ref} -l /groups/nordborg/projects/suecica/007GeneticMap/Parents_goodpositions.filtvcftoolshet.txt $out/$acc.filt.bam > $out/$acc.mpileup
#cd $out
#grep -w "Asue_scaffold1" $acc.mpileup > $acc.mpileup.1
#grep -w "Asue_scaffold2" $acc.mpileup > $acc.mpileup.2 
#grep -w "Asue_scaffold3" $acc.mpileup > $acc.mpileup.3
#grep -w "Asue_scaffold4" $acc.mpileup > $acc.mpileup.4
#grep -w "Asue_scaffold5" $acc.mpileup > $acc.mpileup.5
#grep -w "Asue_scaffold6" $acc.mpileup > $acc.mpileup.6
#grep -w "Asue_scaffold7" $acc.mpileup > $acc.mpileup.7
#grep -w "Asue_scaffold8" $acc.mpileup > $acc.mpileup.8
#grep -w "Asue_scaffold9" $acc.mpileup > $acc.mpileup.9
#grep -w "Asue_scaffold10" $acc.mpileup > $acc.mpileup.10
#grep -w "Asue_scaffold11" $acc.mpileup > $acc.mpileup.11
#grep -w "Asue_scaffold12" $acc.mpileup > $acc.mpileup.12
#grep -w "Asue_scaffold13" $acc.mpileup > $acc.mpileup.13 
#rm $out/$acc.mpileup
