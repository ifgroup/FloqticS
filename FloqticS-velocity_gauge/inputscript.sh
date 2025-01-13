## Requirements: python,fortran compiler, impi, ELPA
##
##
##
## Input description:
## nbands:         number of bands at a k-vector (same for all k-vectors)
## flqfunctions:   number of time periodic function to be taken (total no.:  integers \in [-flqc,flqc])
## Thus, total dimension of the Floquet Hamiltonian for each k-vector will be = nbands*(2*flqc+1)
## edrive:         input drive field amplitude in Volt/Angstrom
## omega:          input drive field photon energy in eV 
## nkvectors:      number of k-vectors in the Brillouin zone
## totvol:         total volume of crystal (in atomic units)
## processors:     number of availabel processors to run the job
##
##
##
##########################################################




########################################################
#################  enter  input here  #######################################
#######################################################

nbands=11
flqfunctions=600
edrive=0.2                        #### in V/A
omega=0.5                         #### in eV
nkvectors=500
totvol=4000.0

processors=40

########################################################
#########################  input ends #######################################
#######################################################



edrive=`python -c "print($edrive/51.4220674763)"`
omega=`python -c "print($omega/27.211386245988)"`

echo $nbands > input.txt
echo $flqfunctions >> input.txt
echo $edrive >> input.txt
echo $omega >> input.txt
echo $nkvectors >> input.txt
echo $totvol >> input.txt


########################################################
###################  enter command to compile the code #######################################
###################  enter the location to access ELPA directory #############################
#######################################################

echo "Job started at"
date 


mpiifort FloqticS.f90 -o compiled -mkl=cluster  -I/software/elpa/2017.05.002/b1/include/elpa-2017.05.002/modules/ -L/software/elpa/2017.05.002/b1/lib -lelpa

mpirun -np $processors ./compiled > spec.txt

rm compiled
rm input.txt

echo "Stick spectra computation completed"
echo "stick spectra in spec.txt" 
date

python3 lorentzfit.py

echo "Absorption coefficient computation completed"
echo "Broadened spectrum in absorptionspectrum.txt"

echo "Thank You for using this code"
date

