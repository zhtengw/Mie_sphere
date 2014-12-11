Mie_eff.for is the main program, edit it to set the parameters, and then compile and link it together with other files, for example:

$ ifort -O3 -o Mie_eff *.for *.f

The code uses getnk subroutine to calculate the refractive index of Au or Ag. 

Before running the program, copy nk(gold).txt or nk(silver).txt to nk.txt to set the metal you use. 

The program will read the discrete data in nk.txt and use the Piecewise Cubic Hermite Interpolating Polynomial (pchip) to get the refractive index of the metal at specific wavelength (eV). The discrete refractive index data in nk(gold).txt and nk(silver).txt comes from P. B. Johnson and R. W. Christy, Phys. Rev. B 6, 4370 (1972). 

Therefore, if you copy the binary program file to another place, make sure you take nk.txt with it.

