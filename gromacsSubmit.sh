#!/bin/bash

# Issue directives to scheduler:

#$ -N GMX
#$ -cwd
#$ -o output.txt
#$ -j y
#$ -pe shared 16
#$ -l mf=1G
#$ -l m_arch=INTEL
#$ -q all.q,all-HiPri.q,all-LoPri.q

# Run commands needed here:

. /etc/profile
module load intel/compiler/64/14.0/2013_sp1.3.174
module load intel/mkl/64/11.1/2013_sp1.3.174
module load gromacs/intel/double/4.6.5

#source /home/dsponsel/bin/gromacs/bin/GMXRC
#mdrun_d -deffnm nvt

#source /home/ncrump/programs/gromacs-4.6.5/build/bin/GMXRC
#mdrun_mpi_d -deffnm nvt

OMP_NUM_THREADS=16
mdrun_d -v -deffnm npt0_T0_P101
