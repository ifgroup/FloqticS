#!/bin/bash

date

procs=2

hdim=6                                            #number of wannier functions
direction=1                #direction of polarization of drive and probe; for x use value = 1,for y=2, for z=3
n_floquet=150              #number of Floquet basis, total Floquet Fourier basis = 2*n_floquet + 1
amplitude=`python -c "print(0.1/51.4220674763)"`  #in atomic units
OMGA=`python -c "print(0.4/27.211386245988)"`     #in atomic units
nblk=64                                           #change only when you know what it means
fermi=`python -c "print(-3.0/27.211386245988)"`   #highest energy of the filled band in atomic units
nkpt=500                                          #number of k points sampling the BZ
total_commutator_terms=22                         #total commutator terms 
unit_cell_volume=`python -c "print(4.7168)"`      #in atomic units, comes from QE, or = a in one-dimension solid

echo $amplitude $total_commutator_terms  $OMGA

echo $hdim > input.txt
echo $direction >> input.txt
echo $n_floquet >> input.txt
echo $amplitude >> input.txt
echo $OMGA >> input.txt
echo $nblk >> input.txt
echo $fermi >> input.txt
echo $nkpt >> input.txt
echo $total_commutator_terms >> input.txt
echo $unit_cell_volume >> input.txt


module load intel
module load elpa
module load impi
module load lapack 
module load blas

date

echo 'electronic structure start'
gfortran  -O3 electronic_structure_generator.f90 -llapack -lblas

 ./a.out

echo 'electronic structure done'
date
rm ./a.out

mpiifort -O3  spectrum_calculator.f90 -o compiled -mkl=cluster  -I/software/elpa/2017.05.002/b1/include/elpa-2017.05.002/modules/ -L/software/elpa/2017.05.002/b1/lib -lelpa

echo 'absorption start'
date

mpirun -np $procs ./compiled > spectrum.txt
rm compiled

echo "absorption end"
date

python3 lorentzfit.py

