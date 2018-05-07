#!/bin/bash
#set -x
########################################################################
#   Grid Engine Job Wrapper
########################################################################
#$ -N POT55M2
#$ -pe orte 50
#$ -q jed
#$ -j y
#####$ -m e
#####$ -M yf160@physics.rutgers.edu
#$ -v LD_LIBRARY_PATH,OMP_NUM_THREADS,MODULEPATH
########################################################################
# DON'T remove the following line!
source $TMPDIR/sge_init.sh
########################################################################
export OMP_NUM_THREADS=1
export SMPD_OPTION_NO_DYNAMIC_HOSTS=1
export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64/:/opt/intel/lib/intel64/:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64/
module load intel/18.0 intel/ompi
gunzip data.tar.gz
tar -xf data.tar
for file in *_0000.dat
	do
		export filebasename=`basename $file _0000.dat`
		echo $filebasename
                mpirun -n $NSLOTS ./get_rho_parallel $filebasename .false. 1.0 0 299 0 0 500 >& Rho$filebasename.scr
		# input, FFT?, Emax, RLZmin, RLZmax, Der=0or1, ForceNc, SetDataPoints
	done
tar -cf data.tar *.dat
tar -rf data.tar *.log
rm *.log
rm *_????.dat
tar -cf data_out_Qoff.tar *.dat
rm *.dat
