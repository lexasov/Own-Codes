#include <cmath>
#include <iostream>
#include <chrono>

double cubickernel(double r, double h)
{

  double u;
  const double pi = 3.141592653589793;
  double W;
  u = r / h;
  if (u >= 0.0 && u <= 1.0)
  {
    W = 1.0 / (pi * h * h * h) * (1 - 1.5e0 * u * u * (1.0e0 - 0.5e0 * u));
  }
  else if (u > 1.0e0 && u < 2.0e0)
  {
    W = 1.0 / (pi * h * h * h) * 0.25e0 * (2.0e0 - u) * (2.0e0 - u) * (2.0e0 - u);
  }
  else
  {
    W = 0.0;
  }

  return W;
}

void gridding3(int npart, double Lbox, double xpos[], double ypos[], double zpos[],
               double vx[], double vy[], double vz[], double ro[], double mass, double partperpixel,
               int npixels, double Gx_1D[], double Gy_1D[], double Gz_1D[])
{

  int i, iii, j, k, n, xindex, yindex, zindex;
  long long int counts;
  int min_intx, max_intx, min_inty, max_inty, min_intz, max_intz;
  double dx, dx2, ydx, x, y, z, r, W, weight, h1, h2, h2_2, s1, s1_2, s2;
  int npixels3 = npixels * npixels * npixels;
  double D[npixels], norm[npixels][npixels][npixels];
  double Gx[npixels][npixels][npixels], Gy[npixels][npixels][npixels], Gz[npixels][npixels][npixels];

  dx = Lbox / npixels;
  ydx = 1.e0 / dx;
  dx2 = dx / 2;
  counts = 0;

  for (i = 0; i < npixels; i++)
  {
    D[i] = (2 * i + 1) * dx2;
  }

  for (i = 0; i < npixels; i++)
  {
    for (j = 0; j < npixels; j++)
    {
      for (k = 0; k < npixels; k++)
      {
        Gx[i][j][k] = 0;
        Gy[i][j][k] = 0;
        Gz[i][j][k] = 0;
        norm[i][j][k] = 0;
      }
    }
  }
  // h1=std::cbrt(partperpixel/npart)*Lbox/4;
  h1 = 1.003 * std::cbrt(3. / 4. / M_PI * partperpixel / npart) * Lbox / 2;
  h2 = 2 * h1;
  h2_2 = h2 * h2;
  for (n = 0; n < npart; n++)
  {
    weight = mass / ro[n];
    if (n % 100000 == 0)
    {
      std::cout << "particle " << n << " of " << npart << std::endl;
    }
    // h1=std::max(h[i],dx2);
    // h2=2*h1;
    // h1=dx2;
    // h2_2=h2*h2;
    max_intz = std::floor((zpos[n] + h2) * ydx - 0.5e0);
    min_intz = std::ceil((zpos[n] - h2) * ydx - 0.5e0);
    for (i = min_intz; i <= max_intz; i++)
    {
      if (i < 0)
      {
        zindex = i + npixels;
      }
      else if (i >= npixels)
      {
        zindex = i - npixels;
      }
      else
      {
        zindex = i;
      }
      z = std::min(std::abs(zpos[n] - D[zindex]), std::abs(zpos[n] + Lbox - D[zindex]));
      z = z * z;
      if (z > h2_2)
        continue;
      s1_2 = h2_2 - z;
      s1 = std::sqrt(s1_2);
      min_intx = std::ceil((xpos[n] - s1) * ydx - 0.5e0);
      max_intx = std::floor((xpos[n] + s1) * ydx - 0.5e0);
      for (j = min_intx; j <= max_intx; j++)
      {
        if (j < 0)
        {
          xindex = j + npixels;
        }
        else if (j >= npixels)
        {
          xindex = j - npixels;
        }
        else
        {
          xindex = j;
        }
        x = std::min(std::abs(xpos[n] - D[xindex]), std::abs(xpos[n] + Lbox - D[xindex]));
        x = x * x;
        if (x > s1_2)
          continue;
        s2 = std::sqrt(s1_2 - x);
        min_inty = std::ceil((ypos[n] - s2) * ydx - 0.5e0);
        max_inty = std::floor((ypos[n] + s2) * ydx - 0.5e0);

        for (k = min_inty; k <= max_inty; k++)
        {
          if (k < 0)
          {
            yindex = k + npixels;
          }
          else if (k >= npixels)
          {
            yindex = k - npixels;
          }
          else
          {
            yindex = k;
          }

          y = std::min(std::abs(ypos[n] - D[yindex]), std::abs(ypos[n] + Lbox - D[yindex]));
          r = std::sqrt(x + y * y + z);

          if (r < h2)
          {
            W = weight * cubickernel(r, h1);
            Gx[xindex][yindex][zindex] = Gx[xindex][yindex][zindex] + W * vx[n];
            Gy[xindex][yindex][zindex] = Gy[xindex][yindex][zindex] + W * vy[n];
            Gz[xindex][yindex][zindex] = Gz[xindex][yindex][zindex] + W * vz[n];
            norm[xindex][yindex][zindex] = norm[xindex][yindex][zindex] + W;
            counts++;
          }
        }
      }
    }
  }
  std::cout << "number of particles contributing to a pixel: " << counts * 1.0 / npixels / npixels / npixels << std::endl;
  for (i = 0; i < npixels; i++)
  {
    for (j = 0; j < npixels; j++)
    {
      for (k = 0; k < npixels; k++)
      {
        iii = i + npixels * (j + npixels * k);
        if (norm[i][j][k] == 0)
        {
          std::cout << "Failed Rasterization" << std::endl;
          std::terminate();
        }
        Gx_1D[iii] = Gx[i][j][k] / norm[i][j][k];
        Gy_1D[iii] = Gy[i][j][k] / norm[i][j][k];
        Gz_1D[iii] = Gz[i][j][k] / norm[i][j][k];
      }
    }
  }
}
