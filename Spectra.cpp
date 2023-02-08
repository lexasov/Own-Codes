#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "gridding_clean.hpp"
#include "file_operations.hpp"
//#include "curlv.hpp"

using namespace std;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  // GENERAL CONFIGURATIONS
  // std::string sphexa_filename = "/home/simseko/Repositories/Data/turb_out_gold.txt";
  // std::string v_filename = "Gv_sphexa_100.txt";
  // std::string spectra_filename = "turb_spectra_100.txt";
  // int n = 100;

  std::string sphexa_filename = "turb_50.txt";
  std::string v_filename = "Gv_sphexa_50.txt";
  std::string spectra_filename = "turb_spectra_50.txt";
  int n = 50;
  double xmin = 0.0, xmax = 1.0;
  double halfbox = 0.5 * (xmax - xmin); // Define parameters
  bool centered_at_0 = true;
  double partperpixel = 8.; // Doesn't work approximately lower than 6
  double Lbox = xmax - xmin;

  size_t npart = std::pow(n, 3);
  double mass = 1. / npart;
  size_t npixels = 2 * n; // updated to 2 since Axel said it is enough
  
  size_t npixels3 = npixels * npixels * npixels;
  double *xpos, *ypos, *zpos, *h, *rho, *vx, *vy, *vz, *p, *divv, *curlv, *c;
  double dummy;
  double *Gv; // [npixels3];

  int Kmax=std::ceil(std::sqrt(3)*npixels*0.5);
  double E[Kmax],k_center[Kmax];

  // Allocate memory
  xpos  = new double[npart];
  ypos  = new double[npart];
  zpos  = new double[npart];
  h     = new double[npart];
  rho   = new double[npart];
//   vx    = new double[npart];
//   vy    = new double[npart];
//   vz    = new double[npart];
  p     = new double[npart];
  divv  = new double[npart];
  curlv = new double[npart];
  c     = new double[npart];
  Gv    = (double *)malloc(npixels3 * sizeof(double));

  auto start = chrono::steady_clock::now();

  read_sphexa_file(sphexa_filename, npart, centered_at_0, xmin, xmax,
                   xpos, ypos, zpos, divv, rho);

  auto end = chrono::steady_clock::now();
  std::cout << "Reading file took: "
            << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            << " ms" << std::endl;

  start = chrono::steady_clock::now();

  gridding_spectra_clean(npart, Lbox, xpos, ypos, zpos, divv, rho, mass, partperpixel, npixels, Gv, E, k_center);
  
  end = chrono::steady_clock::now();
  std::cout << "Gridding took: "
            << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            << " ms" << std::endl;


  start = chrono::steady_clock::now();

  write_spectra_file(v_filename, spectra_filename, npixels, Kmax, Gv, E, k_center);

  end = chrono::steady_clock::now();
  std::cout << "Writing to file took: "
            << chrono::duration_cast<chrono::milliseconds>(end - start).count()
            << " ms" << std::endl;

  delete[] xpos; 
  delete[] ypos; 
  delete[] zpos; 
  delete[] h;     
  delete[] rho;  
//   delete[] vx;   
//   delete[] vy;   
//   delete[] vz;   
  delete[] p;    
  delete[] divv; 
  delete[] curlv;
  delete[] c;    
  delete[] Gv;   

  MPI_Finalize();
}