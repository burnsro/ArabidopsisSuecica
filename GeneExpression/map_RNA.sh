#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=16  # 14 physical cores per task
#SBATCH --mem=64G   # 64GB of RAM
#SBATCH --qos=short
#SBATCH --time=0-06:00:00
#SBATCH --output=%A_%a.STAR.stdout
#SBATCH --array=1-20

CONDA_ENVS_PATH='/groups/nordborg/projects/suecica/005scripts/001Software'
CONDA_PKGS_DIRS=/opt/anaconda/pkgs='/groups/nordborg/projects/suecica/005scripts/001Software'

cd /groups/nordborg/projects/suecica/005scripts/001Software
ml anaconda3/2019.03
source activate RobinCondaSCRATCH/
ml bamtools/2.5.1-foss-2018b

raw='/groups/nordborg/projects/suecica/006RNAseq/001fastq'
out='/scratch-cbe/users/robin.burns/012RNA/mapAS'
temp='/scratch-cbe/users/robin.burns/tmp'
samples=$raw'/mixparents'
acc=$(awk "NR==$SLURM_ARRAY_TASK_ID" $samples)
echo $acc
ann='/groups/nordborg/projects/nordborg_rawdata/Asuecica/RobinPacbio/Asue_annot/June2020/proteinortho/augustus_hints3sourcesAs_June2020.aug.renamed.gff'
ref='/groups/nordborg/projects/suecica/001Assembly/004Asuecica/Asue_genome_210620.full.fasta'


mkdir -p $out/Asue_$acc
cd $out/Asue_$acc


#mkdir $out/mygenome
#STAR --runThreadN 16 \
#  --runMode genomeGenerate \
#  --genomeDir $out/mygenome \
#  --genomeFastaFiles ${ref}  \
#  --sjdbGTFfile ${ann} \
#  --sjdbOverhang 124 \
#  --genomeSAindexNbases 12

STAR --runThreadN 16 \
	--genomeDir ${out}/mygenome \
	--readFilesIn $raw/$acc.1.fastq.gz $raw/$acc.2.fastq.gz \
	--readFilesCommand zcat \
       --outFilterMultimapNmax 1 \
	--outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix ${acc}_STARPaired1hitSorted \
	--outReadsUnmapped Fastx \
	--outSAMorder Paired \
	--outSAMprimaryFlag OneBestScore \
	--quantMode GeneCounts \
  --limitBAMsortRAM 6804733228	


