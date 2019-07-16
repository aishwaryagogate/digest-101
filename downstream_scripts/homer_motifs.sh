#!/bin/bash

#SBATCH --job-name=homer-motif                           # job name
#SBATCH --partition=super                                # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                        # number of nodes requested by user
#SBATCH --time=0-48:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=serialJob.%j.out                        # standard output file name
#SBATCH --error=serialJob.%j.time                        # standard error output file name
#SBATCH --mail-user=aishwarya.gogate@utsouthwestern.edu      # specify an email address
#SBATCH --mail-type=ALL                                  # send email when job status change (start, end, abortion and etc.)


PATH=$PATH:/project/GCRB/Banaszynski_lab/shared/software/homer/bin
cd /project/GCRB/Banaszynski_lab/shared/software/homer/bin

./changeNewLine.pl /project/GCRB/Banaszynski_lab/shared/path-to-your-peaks-file/filename.bed
./findMotifsGenome.pl /project/GCRB/Banaszynski_lab/shared/path-to-your-peaks-file/filename.bed /project/GCRB/Banaszynski_lab/shared/software/homer/data/genomes/mm10 /project/GCRB/Banaszynski_lab/shared/path-to-output-folder/outputfolder -size -500,500 -mask
