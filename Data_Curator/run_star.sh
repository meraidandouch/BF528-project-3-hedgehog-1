#!/bin/bash

# Modules loaded and organized in conda environment

for var in SRR1177966 SRR1177969 SRR1177970 SRR1177993 SRR1177994 SRR1177995 SRR1177998 SRR1178001 SRR1178003

do
	STAR --readFilesIn samples/${var}_1.fastq.gz samples/${var}_2.fastq.gz --genomeDir /project/bf528/project_3/reference/rn4_STAR --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix star_output/${var}
done
