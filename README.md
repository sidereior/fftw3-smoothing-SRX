<h1> This is part of @mesuOSU 's M-THEMES Project. Read through my project description and walkthrough here: https://tinyurl.com/alexFFTW </h1>

(1)This model includes the lattice resistance term (tau_0) in the velocity of dislocations, currently it is set as zero, but if needed then please change the value of tau0xf term in the Cu_p.sx file, the code reads this term as the lattice resistance stress.
(2)The elastic constants are hard coded inside the code in the init.c file. The elastic constants for Cu and Al are given as a temperature dependent term. So, whenever you want to use other metals, please change the elastic constants inside the init.c file, changing in input file won't work. (Update) Or you cn switch off the elastic constants hard coded in the init.c, then it can read from the input file.
(3)To run the code with parallel processor, make sure dimension along X axis of your microstructure (CellDim[0] in code) must be divisible by the total number of processors used. For example, the X dimension is 64, then you can use 2,4,16,32, or 64 processors to run the code, but not with 3,9,21, or 30 processors. This is done in this way because the messege passing (MPI) is done with blocks of data and the blocks are of same of memory size and segmented along the X axis.
This is one of the way of doing, you can change the way want to send and receive the blocks, for that you need to change the code.
(4)For DFT calculations, select the Cell dimensions such a way that it can be factorized into the lowest prime number (i.e. 2). For example, dimensions like 64*64*64 or 128*128*128 works best with DFT. However, 27*27*27, 125*125*125 or 125*27*27 will also work (maybe less accuracy). 
Deformed positions can be visualized using Paraview and the vectors written in the new_position_...iout files. convert this particular iout file using position.c code, as following, 
mpicc position.c -o position.out
./position.out new_position...iout new_vector...txt
use new vectors in you vtk file/HDF5 files. And use "wrap by vectors" in paraview
