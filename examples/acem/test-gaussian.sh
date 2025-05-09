#!/bin/bash

#SBATCH --job-name=test-gaussian
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --mail-type=NONE
#SBATCH --mem-per-cpu=2200
# 
hostname

 
#set -xv
export g16root=/cluster/software/gaussian
export GAUSS_SCRDIR=/scratch/boittier/test-gaussian
mkdir -p /scratch/boittier/test-gaussian
source $g16root/g16/bsd/g16.profile

$g16root/g16/g16 /pchem-data/meuwly/boittier/home/acem/test-gaussian.com /scratch/boittier/test-gaussian/test-gaussian.out

formchk test.chk test.fchk 
cubegen 0 Density test.fchk Density.cube -3 h
cubegen 0 Potential test.fchk ESP.cube -3 h


# don't delete the result file if not able to copy to fileserver 
cp /scratch/boittier/test-gaussian/test-gaussian.out /pchem-data/meuwly/boittier/home/acem/test-gaussian.out 
status=$?
if [ $status -eq 0 ] 
then 
   rm -rf /scratch/boittier/test-gaussian
else
   host=`/bin/hostname`
   /usr/bin/Mail -v -s "Error at end of batch job" $USER <<EOF

At the end of the batch job the system could not copy the output file
	$host:/scratch/boittier/test-gaussian/test-gaussian.out
to
	/pchem-data/meuwly/boittier/home/acem/test-gaussian.out
Please copy this file by hand or inform the system manager.

EOF
 
fi
