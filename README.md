# FloqticS code
 Floquet theory and computational method for the optical absorption of laser-dressed solids
author: Vishal Tiwari, Bing Gu, Ignacio Franco

Code: VIshal TIwari
 

 
This is the README file containing the information to sucessfully run the code.

The code takes the k-vectors $`(\mathbf{k})`$, band energies at corresponding k-vector $`(\epsilon_{u\mathbf{k}})`$, momentum matrix elements (MME) at corresponding k-vectors among the Bloch states along the direction of polarization of drive laser $`(\hat{\mathbf{e}}_{\text{d}} \cdot \mathbf{p}_{u\mathbf{k},v\mathbf{k}})`$ and probe laser $`(\hat{\mathbf{e}}_{\text{p}} \cdot \mathbf{p}_{u\mathbf{k},v\mathbf{k}})`$ at corresponding  k-vectors, and initial occupation of bands ($`\bar{n}_{u\mathbf{k}}`$) and provides the laser-dressed absorption spectrum $`(A(\omega))`$, that is absorption coefficient as a function of probe laser photon energy. For meaning of notations see paper.
 
 
 
Requirements: Python (only for postprocessing), FORTRAN compiler, Intel MPI, [Eigenvalue soLver for Petaflops Application (ELPA) package](https://elpa.mpcdf.mpg.de/)


 Input files description:
 1) `kvector.txt`         kvectors (in atomic units) format -->  kx | ky | kz |
 2) `bandeng.txt `        energy (in atomic units) of all bands in a single row 
 3) `dipoledrive.txt`     MME along drive laser polarization (in atomic units) 
 4) `dipoleprobe.txt `    MME along probe laser polarization (in atomic units)
 5) `initialoccup.txt`    the initial occupation of the bands (either 0 (empty) or 1 (filled)), $`\bar{n}_{u\mathbf{k}}`$, ordered according to the k-vectors
 6) `inputscript.sh`      Input script to specify the  parameters of system and laser
 
 For eg. if you have 11 bands and 500 k-vectors:
    `kvector.txt` will have 3 columns and 500 rows, each row one k-vector
    `bandeng.txt` will have 500 rows and 11 columns, each row has energy of 11 bands for a k-vector
    `dipoledrive.txt` will have 11x500 rows and 11 columns, first 11 rows and 11 columns represent MME for a k-vector
    `dipoleprobe.txt` will have 11x500 rows and 11 columns, first 11 rows and 11 columns represent MME for a k-vector
    `initialoccup.txt` will have 500 rows and 11 columns, the occupation in each row must be ordered according to the band <br>
 MME are complex numbers so if a MME is for e.g. 2+3i  then it should be written in input file dippoledrive.txt as (2,3). This is the FORTRAN complex number format. <br> Make sure to enter the correct command to compile the code using ELPA and mpiifort compiler or equivalent in inputscript.sh file according to your computer and location of ELPA directory. 
 
 Output files description: 
 1) `band.txt` - contains the band energies (in eV) you have given corresponding to the k-vectors (for sanity check) 
    format 
    kx | ky | kz | energy for band 1 | energy for band 2 | ...  <br>
 2) `quasienergy.txt` - contains the quasienergies (in eV) withing the fundamental FBZ. The number of quasienergies is equal to the number of bands in the computation and arranged according to the k-vectors
    Format:
    kx | ky | kz | quasienergy for mode 1 | quasienergy for mode 2 | $` ... `$    <br>
 3) `spec.txt` -  contains the laser-dressed spectrum for the drive amplitude and photon energy specified in the following format:
    $`kx \hspace{0.5cm} |  \hspace{0.5cm} ky \hspace{0.5cm} | \hspace{0.5cm} kz \hspace{0.5cm} |\hspace{0.5cm}  \hbar\omega \hspace{0.5cm} |\hspace{0.5cm}  A(\omega) \hspace{0.5cm} |  \hspace{0.5cm} \alpha \hspace{0.5cm} | \hspace{0.5cm} \beta \hspace{0.5cm} | \hspace{0.5cm} n \hspace{0.5cm} |\hspace{0.5cm}  E_{\alpha}\hspace{0.5cm}  | \hspace{0.5cm} E_{\beta} \hspace{0.5cm} |  \hspace{0.5cm} \hspace{0.5cm} |\hspace{0.5cm} \mathcal{P}^{(n)}_{\alpha \beta}|^{2} \hspace{0.5cm} | \hspace{0.5cm} \Lambda_{\alpha \beta} \hspace{0.5cm} |\hspace{0.5cm}  \Lambda_{\beta \alpha} \hspace{0.5cm} `$ <br>
    The transitions will be ordered according to the k-vectors. <br>
    $`\mathcal{P}`$ is the Fourier component of MME, $` \Lambda `$ is the population factor. For exact meaning of the above variables refer to the original paper.
    If there are errors in running the code they will also be printed in spec.txt file <br>
    Units: 
    $`kx,ky,kz`$                   -->  Å$`^{-1} `$ <br>
    $`\hbar\omega   `$             --> eV  <br>
    $`A(\omega)  `$                --> $`m^{-1} `$ <br>
    $`E_{\alpha,\beta} `$            --> eV  <br>
    $`|P^{(n)}_{\alpha \beta}|^{2}`$ --> (eV fs Å$`^{-1})^{2} `$ <br>
 
 The k-vectors are not used in the calculation but instead used as a index for transitions
 
 Protocol to follow:
 1) Generate the band energies and MME in direction of probe and drive laser for the solid.
 2) Arrange the band energies to create the bandeng.txt file as mentioned above.
 3) Arrange the MME input files, dipoledrive.txt and dipoleprobe.txt as mentioned above. (make sure the name of input files are what are mentioned above)
 4) Enter the system and laser parameters in the inputscript.sh file.
 5) Run the code, the command to compile needs to be adjusted according to your machine, an example method is mentioned in the `inputscript.sh` file.
 6) Make sure the quasienergies obtained are converged w.r.t the number of time-periodic functions  and the absorption spectra is converged w.r.t the number of bands and k-vectors.
 7) You can run `lorentzfit.py` code to obtain the absorption spectrum where peaks are broadened using a Lorentzian function or use your own broadening function as well. The input  describing the broadening and the grid for spectrum is to be input in the file `lorentzfit.py`.

The folder `test/` contains example input files used to compute the optical absorption in the paper for $`E_{d}=0.2`$ V/Å  and $`\hbar\Omega=0.5`$ eV with 11 bands and 500 k-points in BZ.

For further questions or issues regarding the code please contact: Vishal Tiwari (vtiwari3@ur.rochester.edu)
