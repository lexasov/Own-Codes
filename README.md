# Own-Codes
Uploaded Personal Codes

The Rasterization Code is purposed to turbulence tests with roughly uniform density in the box.

Some prerequisites to the newly added heFFTe library to compute the 3D FFT. -> (https://github.com/af-ayala/heffte)

The following cmake command is what I used for building the heFFTe library:
cmake -D CMAKE_BUILD_TYPE=Debug -D BUILD_SHARED_LIBS=ON -D CMAKE_INSTALL_PREFIX=!!!!/PATH/TO/INSTALL!!!! -D Heffte_ENABLE_AVX=ON -D Heffte_ENABLE_FFTW=ON ..
make
make install

heFFTe depends on fftw library if one wants to run it on the CPU.
How to install in Ubuntu: sudo apt-get install libfftw3-dev

heFFTe also depends on MPI.
Ubuntu install MPI: sudo apt install mpich

After the installation is done, either add the heffte library path to the PATH environment variable or export it to LD_LIBRARY_PATH when running the Rasterization code.

Makefile includes the compilation command.

Note: The Rasterization code requires unlimited stack in the operating system for the time being. We can fix that later.
(It can be set to unlimited with the following command:  ulimit -s unlimited)