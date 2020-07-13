#!/bin/bash

#!/bin/bash
#SBATCH -J lionfish_bam2fastx
#SBATCH -o lionfish.out
#SBATCH -e lionfish.err
#SBATCH -N 4
#SBATCH -n 8
#SBATCH -t 00:30:00
#SBATCH -p normal

module load smrtanalysis

bam2fastq -o m54336U_200609_173116.subreads m54336U_200609_173116.subreads.bam
