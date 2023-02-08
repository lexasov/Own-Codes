#include <iostream>
#include <fstream>
#include <cmath>
#include "hdf5.h"

#define "turb_50.h5"

void read_sphexa_file(std::string in_filename, uint64_t npart, bool is_centered_at_zero, double xmin, double xmax,
                     double *xpos, double *ypos, double *zpos, double *divv, double *rho)
{
  double halfbox = 0.5 * (xmax - xmin);
  double dummy;
  std::ifstream infile;

// For reading vx, vy, vz from the file
//   double *vx, *vy, *vz;

//   vx = new double[npart];
//   vy = new double[npart];
//   vz = new double[npart];

  infile.open(in_filename);

  if (!infile.is_open())
  {
    std::cout << "Input file not found" << std::endl;
    exit(0);
  }
  std::cout << "Reading Data from file..." << std::endl;

  for (uint64_t k = 0; k < npart; k++)
  {
    // The file has the following values: x,y,z,h,rho,vx,vy,vz,p,divv,curlv,c
    infile >> xpos[k] >> ypos[k] >> zpos[k] >> dummy >> rho[k] >> dummy >> dummy >> dummy >> dummy >> divv[k] >> dummy >> dummy;

    if (is_centered_at_zero)
    {
      xpos[k] = xpos[k] + halfbox;
      ypos[k] = ypos[k] + halfbox; // Griding is done in box starting at xmin=0;
      zpos[k] = zpos[k] + halfbox;
    }
  }

  infile.close();
// For reading vx, vy, vz from the file
//   delete[] vx;
//   delete[] vy;
//   delete[] vz;
}

void write_spectra_file(std::string v_filename, std::string spectra_filename, uint64_t npixels, int Kmax,
                   double* Gv, double* E, double* k_center)
{
  uint64_t kstart = 0;
  uint64_t kend = npixels * npixels * npixels;
  kstart = npixels / 2 * (npixels * npixels);
  kend = (npixels / 2 + 1) * (npixels * npixels);

  std::cout << "Writing in file..." << std::endl;

  // std::ofstream vfile;
  // vfile.open(v_filename, std::ios::trunc);
  // for (uint64_t k = kstart; k < kend; k++)
  // {
  //   if (k % 1000000 == 0)
  //   {
  //     std::cout << "cell " << k << " of " << npixels * npixels << std::endl;
  //   }
  //   vfile << std::setprecision(8) << std::scientific << Gv[k] << std::endl;
  //   // divfile << std::setprecision(8) << std::scientific << Gv2[k] << std::endl;
  //   // curlfile << std::setprecision(8) << std::scientific << Gv3[k] << std::endl;
  //   // rhofile << std::setprecision(8) << std::scientific << Gv4[k] << std::endl;
  // }
  // vfile.close();

  std::ofstream spectra;
  spectra.open(spectra_filename, std::ios::trunc);
  for (int i = 0; i < Kmax; i++)
  {
    spectra << std::setprecision(8) << std::scientific << k_center[i] << ' ' << E[i] << std::endl;
  }
  spectra.close();
}