#!/bin/bash
#set -x
########################################################################
#   Grid Engine Job Wrapper
########################################################################
#$ -N L55Q4-
#$ -pe orte 50
#$ -q jed
#$ -j y
#$ -m e
#$ -M yf160@physics.rutgers.edu
#$ -v LD_LIBRARY_PATH,OMP_NUM_THREADS,MODULEPATH
########################################################################
# DON'T remove the following line!
source $TMPDIR/sge_init.sh
########################################################################
export OMP_NUM_THREADS=1
export SMPD_OPTION_NO_DYNAMIC_HOSTS=1
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64/:/opt/intel/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/
module load intel/18.0 intel/ompi
source makeinput
for file in *.in
	do
		echo $file
                mpirun -n $NSLOTS ./get_moment $file 6 >& $file.scr
		rm $file
	done
tar -cf data.tar *.dat
tar -rf data.tar *.log
#######tar -cf data_out.tar ?.??.dat
rm *.dat
rm *.log
