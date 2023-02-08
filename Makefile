CC = g++
CFLAGS = -Wall -g
#-DWRITE_TO_FILE
# Also export heffte to ld_lib_path
 
rasterize: 
	mpicxx Spectra.cpp gridding_clean.hpp -o spectra -fopenmp -g -DWRITE_TO_FILE -I/home/simseko/Repositories/lib_installed/. -I/usr/include/mpi/. -L/home/simseko/Repositories/lib_installed/lib/. -lheffte -L/usr/lib/x86_64-linux-gnu/. -lmpi -lfftw3 -lfftw3f

clean:
	rm -f rasterize spectra ./*.o ./*.gch