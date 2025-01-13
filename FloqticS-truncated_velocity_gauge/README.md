# FloqticS code
 First-principle based Floquet engineering of solids in the velocity gauge
author: Vishal Tiwari, Ignacio Franco

Code: VIshal TIwari 
 

 
This is the README file containing the information to sucessfully run the code.

#The code requires the following input:
Number of processors
k-vectors $`(\mathbf{k})`$
Wannier function hoppings
Wannier function position operator matrix elements
Direction of polarization of drive laser $`(\hat{\mathbf{e}}_{\text{d}} \cdot \mathbf{p}_{u\mathbf{k},v\mathbf{k}})`$ and probe laser $`(\hat{\mathbf{e}}_{\text{p}} \cdot \mathbf{p}_{u\mathbf{k},v\mathbf{k}})`$ 
Lattice vectors
Reciprocal space vectors
Volume of unit cell
Field amplitude
Field photon energy
Fermi energy for filled bands

It outputs the field free spectra and the  laser-dressed absorption spectrum $`(A(\omega))`$, that is absorption coefficient as a function of probe laser photon energy. 
 
 
 
#Requirements: Python, FORTRAN compiler, Intel MPI, [Eigenvalue soLver for Petaflops Application (ELPA) package](https://elpa.mpcdf.mpg.de/)


# Input files description:
 1) `kpoints.txt`       kvectors (in atomic units) format -->  kx | ky | kz |
 2) `parameter.dat `    Wannier function hopping and dipole matrix elements from Wannier90 `wann_hr.dat` and `wann_r.dat` files
 3) `submit.sh`         Input script to specify the  parameters of system and laser
 4) `electronic_structure_generator` Code that generates intermediary inputs, band energies and truncated momentum matrix
 5) `spectrum_calculator`    Code that calculates the absorption spectrum of laser-dressed solid for given drive laser
 
The `parameter.dat` file contains the real part of the tight-binding parameters of the solid considered. It is created by combining the output from Wannier90 files `wann_hr.dat` and `wann_r.dat` in the format
 $` x_index y_index z_index Wannier_m Wannier_n Hopping_real r_x_real r_y_real r_z_real`$

 For eg. if you have 6 Wannier functions, 4 nearest neighbour hoppings and 500 k-vectors:
    `kvector.txt` will have 3 columns and 500 rows, each row one k-vector
    `parameter.txt` will have 6*6*(4+4+1) rows and 9 columns. 

Make sure to enter the correct command to compile the code using ELPA and mpiifort compiler or equivalent in submit.sh file according to your terminal and location of ELPA directory. 
 
 Output files description: 
 1) `bandeng.txt` - contains the band energies (in eV) you have given corresponding to the k-vectors (for sanity check) 
    format 
    kx | ky | kz | energy for band 1 | energy for band 2 | ...  <br>
 2) `spectrum.txt` -  contains the laser-dressed spectrum for the drive amplitude and photon energy specified in the following format:
    $`kx \hspace{0.5cm} |  \hspace{0.5cm} ky \hspace{0.5cm} | \hspace{0.5cm} kz \hspace{0.5cm} |\hspace{0.5cm}  \hbar\omega \hspace{0.5cm} |\hspace{0.5cm}  A(\omega) \hspace{0.5cm} |  \hspace{0.5cm} \alpha \hspace{0.5cm} | \hspace{0.5cm} \beta \hspace{0.5cm} | \hspace{0.5cm} n \hspace{0.5cm} |\hspace{0.5cm}  E_{\alpha}\hspace{0.5cm}  | \hspace{0.5cm} E_{\beta} \hspace{0.5cm} |  \hspace{0.5cm} \hspace{0.5cm} |\hspace{0.5cm} \mathcal{P}^{(n)}_{\alpha \beta}|^{2} \hspace{0.5cm} | \hspace{0.5cm} \Lambda_{\alpha \beta} \hspace{0.5cm} |\hspace{0.5cm}  \Lambda_{\beta \alpha} \hspace{0.5cm} `$ <br>
    The transitions will be ordered according to the k-vectors. <br>
    $`\mathcal{P}`$ is the Fourier component of MME, $` \Lambda `$ is the population factor. For exact meaning of the above variables refer to the original paper.
    If there are errors in running the code they will also be printed in spectrum.txt file <br>
    Units: 
    $`kx,ky,kz`$                   -->  Å$`^{-1} `$ <br>
    $`\hbar\omega   `$             --> eV  <br>
    $`A(\omega)  `$                --> $`cm^{-1} `$ <br>
    $`E_{\alpha,\beta} `$            --> eV  <br>
    $`|P^{(n)}_{\alpha \beta}|^{2}`$ --> (eV fs Å$`^{-1})^{2} `$ <br>
 3) `fieldfreespectra.txt` -  contains the spectrum of pristine solid in the format
    $`kx \hspace{0.5cm} |  \hspace{0.5cm} ky \hspace{0.5cm} | \hspace{0.5cm} kz \hspace{0.5cm} |\hspace{0.5cm}  \hbar\omega \hspace{0.5cm} |\hspace{0.5cm}  A(\omega) \hspace{0.5cm} `$ <br>
 
 
 Protocol to follow:
 1) Make sure you have kpoints.txt and parameter.dat file for inputs
 2) Enter details in the submit.sh file
 3) Enter the lattice vectors and reciprocal space lattice vectors in electronic\_structure\_generator.f90 file
 4) Run the code, the command to compile needs to be adjusted according to your machine, an example method is mentioned in the `submit.sh` file.
 5) Make sure the quasienergies obtained are converged w.r.t the number of Floquet Fourier basis, number of commutator terms and so does the  absorption spectrum is converged
 6) You can run `lorentzfit.py` code to obtain the absorption spectrum where peaks are broadened using a Lorentzian function or use your own broadening function as well. The input  describing the broadening and the grid for spectrum is to be input in the file `lorentzfit.py`.

The folder `test/` contains example input files used to compute the optical absorption in the paper for $`E_{d}=0.1`$ V/Å  and $`\hbar\Omega=0.4`$ eV with 6 Wannier functions bands and 500 k-points in BZ.
To run,> ./submit.sh

For further questions or any issues regarding any part of the code please feel free to contact: Vishal Tiwari (vtiwari3@ur.rochester.edu) or Ignacio Franco (ifranco@rochester.edu)
