#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Gridding.hpp"
//#include "curlv.hpp"

int main(){

int n= 100;
double xmin=0.0,xmax=1.0;
double halfbox=0.5*(xmax-xmin);           // Define parameters
bool centered_at_0=true;
double partperpixel=8.; // Doesn't work approximately lower than 6

double Lbox=xmax-xmin;
int npart=std::pow(n,3);
double mass=1./npart;
int npixels=2*n;
int npixels3=npixels*npixels*npixels;
double xpos[npart],ypos[npart],zpos[npart],vx[npart],vy[npart],vz[npart],v[npart],
       h[npart],u[npart],ro[npart],radius[npart],curlv[npart],p[npart],cs[npart],
       divv[npart], ax[npart], ay[npart], az[npart];
double Gvx[npixels3],Gvy[npixels3],Gvz[npixels3],Gv[npixels3];
double xcenters[npixels3],ycenters[npixels3],zcenters[npixels3],centers[npixels3];
int neighbours[npart];

std::ifstream infile;
std::ofstream v_xyz;
std::ofstream vfile;

infile.open("sphexa_100_18s.txt"); // N=1024000
//v_xyz.open("v_xyz.txt",std::ios::trunc);
vfile.open("v_sphexa_100_18s.txt",std::ios::trunc);

if(!infile.is_open()){
        std::cout << "Input file not found" << std::endl;
        return -1;
    }
std::cout << "Reading Data from file..." << abs(1.535) << ' ' << std::abs(1.535) << std::endl;
for(int k=0;k<npart;k++){
  //infile >> n >> xpos[k] >> ypos[k] >> zpos[k] >> h[k] >> u[k] >> ro[k] >> vx[k] >> vy[k] >> vz[k];
  //infile.ignore(13*29); // SPHYNX
  infile >>  xpos[k] >> ypos[k] >> zpos[k] >> h[k] >>  ro[k] >> vx[k] >> vy[k] >> vz[k] >> cs[k] >> p[k] >> u[k] >> neighbours[k]
         >> divv[k] >> curlv[k] >> ax[k] >> ay[k] >> az[k]; //sphexa
  //infile >>  xpos[k] >> ypos[k] >> zpos[k] >> u[k] >> ro[k] >> vx[k] >> vy[k] >> vz[k] ; //sphexa
  v[k]=std::sqrt(vx[k]*vx[k]+vy[k]*vy[k]+vz[k]*vz[k]);

  if(centered_at_0){
    xpos[k]=xpos[k]+halfbox;
    ypos[k]=ypos[k]+halfbox;  // Griding is done in box starting at xmin=0;
    zpos[k]=zpos[k]+halfbox;

  }
}

  gridding(npart,Lbox,xpos,ypos,zpos,v,ro,mass,partperpixel,npixels,centers,Gv);
  //gridding3(npart,Lbox,xpos,ypos,zpos,vx,vy,vz,promro,mass,partperpixel,npixels,xcenters,ycenters,zcenters,Gvx,Gvy,Gvz);

for(int k=0;k<npixels3;k++){

  //v_xyz << std::setprecision(16) << std::scientific <<  Gvx[k] << "  " << Gvy[k] << "   " << Gvz[k] << "   " << std::endl;
  vfile << std::setprecision(8) << std::scientific << Gv[k] << std::endl;
}

}
