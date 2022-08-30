#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "Gridding.hpp"
//#include "curlv.hpp"


using namespace std;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int n = 100;
  double xmin = 0.0, xmax = 1.0;
  double halfbox = 0.5 * (xmax - xmin); // Define parameters
  bool centered_at_0 = true;
  double partperpixel = 8.; // Doesn't work approximately lower than 6

  double Lbox = xmax - xmin;
  int npart = std::pow(n, 3);
  double mass = 1. / npart;
  int npixels = 2 * n; // updated to 2 since Axel said it is enough
  int npixels3 = npixels * npixels * npixels;
  // double xpos[npart], ypos[npart], zpos[npart], v[npart], ro[npart], vx[npart], vy[npart], vz[npart];
  double *xpos, *ypos, *zpos, *v, *ro, *vx, *vy, *vz;
  double dummy;
  // double h[npart],u[npart],radius[npart],curlv[npart],p[npart],cs[npart],divv[npart], ax[npart], ay[npart], az[npart];
  // double Gvx[npixels3],Gvy[npixels3],Gvz[npixels3];
  double Gv[npixels3];

  int Kmax=std::ceil(std::sqrt(3)*npixels*0.5);
  double E[Kmax],k_center[Kmax];

  std::ifstream infile;
  std::ofstream v_xyz;
  std::ofstream vfile;

  // Allocate memory
  xpos = new double[npart];
  ypos = new double[npart];
  zpos = new double[npart];
  v    = new double[npart];
  ro   = new double[npart];
  vx   = new double[npart];
  vy   = new double[npart];
  vz   = new double[npart];

  infile.open("sphexa_100_18s.txt"); // N=1024000
  // v_xyz.open("v_xyz.txt",std::ios::trunc);
  vfile.open("v_sphexa_100_18s_2.txt", std::ios::trunc);

  if (!infile.is_open())
  {
    std::cout << "Input file not found" << std::endl;
    return -1;
  }
  std::cout << "Reading Data from file..." << std::endl;
  auto start = chrono::steady_clock::now();
  for (int k = 0; k < npart; k++)
  {
    // infile >> n >> xpos[k] >> ypos[k] >> zpos[k] >> h[k] >> u[k] >> ro[k] >> vx[k] >> vy[k] >> vz[k];
    // infile.ignore(13*29); // SPHYNX
    // infile >>  xpos[k] >> ypos[k] >> zpos[k] >> h[k] >>  ro[k] >> vx[k] >> vy[k] >> vz[k] >> cs[k] >> p[k] >> u[k] >> neighbours[k]
    //        >> divv[k] >> curlv[k] >> ax[k] >> ay[k] >> az[k]; //sphexa
    infile >> xpos[k] >> ypos[k] >> zpos[k] >> dummy >> ro[k] >> vx[k] >> vy[k] >> vz[k] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy; // sphexa
    // infile >>  xpos[k] >> ypos[k] >> zpos[k] >> u[k] >> ro[k] >> vx[k] >> vy[k] >> vz[k] ; //sphexa
    v[k] = std::sqrt(vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k]);

    if (centered_at_0)
    {
      xpos[k] = xpos[k] + halfbox;
      ypos[k] = ypos[k] + halfbox; // Griding is done in box starting at xmin=0;
      zpos[k] = zpos[k] + halfbox;
    }
  }
  auto end = chrono::steady_clock::now();

  std::cout << "Reading file took: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

  start = chrono::steady_clock::now();
  gridding(npart, Lbox, xpos, ypos, zpos, v, ro, mass, partperpixel, npixels, Gv);
  // gridding3(npart,Lbox,xpos,ypos,zpos,vx,vy,vz,promro,mass,partperpixel,npixels,Gvx,Gvy,Gvz);
  end = chrono::steady_clock::now();
  std::cout << "Gridding took: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;

  

#ifdef WRITE_TO_FILE
  std::cout << "Writing in file..." << std::endl;

  start = chrono::steady_clock::now();
  for (int k = 0; k < npixels3; k++)
  {
    if (k % 1000000 == 0)
    {
      std::cout << "cell " << k << " of " << npixels3 << std::endl;
    }
    // v_xyz << std::setprecision(16) << std::scientific <<  Gvx[k] << "  " << Gvy[k] << "   " << Gvz[k] << "   " << std::endl;
    vfile << std::setprecision(8) << std::scientific << Gv[k] << std::endl;
  }
  end = chrono::steady_clock::now();
  std::cout << "Writing to file took: "
        << chrono::duration_cast<chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;
#endif

  delete[] xpos;
  delete[] ypos;
  delete[] zpos;
  delete[] v;
  delete[] ro;
  delete[] vx;
  delete[] vy;
  delete[] vz;

  MPI_Finalize();
}

int main1(int argc, char** argv)
{
  MPI_Init(&argc, &argv);

  try_fft();

  MPI_Finalize();

  return 0;
}