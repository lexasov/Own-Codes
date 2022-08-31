#include <cmath>
#include <iostream>
#include <chrono>
#include "include/heffte.h"
#include <vector>

void try_fft()
{
  // Create a 4x4x4 box
  heffte::box3d<> inbox = {{0, 0, 0}, {3, 3, 3}};
  heffte::box3d<> outbox = {{0, 0, 0}, {3, 3, 3}};

  // Create 8 2x2x2 boxes
  // heffte::box3d<> inbox1  = {{0, 0, 0}, {1, 1, 1}};
  // heffte::box3d<> inbox2  = {{0, 0, 2}, {1, 1, 3}};
  // heffte::box3d<> inbox3  = {{0, 2, 0}, {1, 3, 1}};
  // heffte::box3d<> inbox4  = {{2, 0, 0}, {3, 1, 1}};
  // heffte::box3d<> inbox5  = {{0, 2, 2}, {1, 3, 3}};
  // heffte::box3d<> inbox6  = {{2, 0, 2}, {3, 1, 3}};
  // heffte::box3d<> inbox7  = {{2, 2, 0}, {3, 3, 1}};
  // heffte::box3d<> inbox8  = {{2, 2, 2}, {3, 3, 3}};
  std::array<int, 3> grid = {2, 2, 2};
  std::vector<heffte::box3d<>> all_boxes = heffte::split_world(inbox, grid);
  std::vector<heffte::box3d<>> reordered = heffte::reorder(all_boxes, {1, 2, 0});

  heffte::fft3d<heffte::backend::fftw> fft(inbox, outbox, MPI_COMM_WORLD);

  // vectors with the correct sizes to store the input and output data
  // taking the size of the input and output boxes
  std::vector<double> input(fft.size_inbox());
  std::vector<std::complex<double>> output(fft.size_outbox());

  // fill the input vector with data that looks like 0, 1, 2, ...
  std::iota(input.begin(), input.end(), 1); // put some data in the input
  // for(int i = 0; i < 4; i++) {
  //   input.at(i) = i+1;
  // }

  // input.at(0) = 1;
  // input.at(1) = 2;
  // input.at(2) = 3;
  // input.at(3) = 4;
  // input.at(4) = 1;
  // input.at(5) = 2;
  // input.at(6) = 3;
  // input.at(7) = 4;

  // for(int i = 0; i < 4; i++)
  // {
  //   std::cout << input.at(i) << std::endl;
  // }

  // perform a forward DFT
  fft.forward(input.data(), output.data());

  auto real_inverse = fft.backward_real(output);

  for (int i = 0; i < 64; i++)
  {
    std::cout << input.at(i) << std::endl;
  }
  std::cout << "================" << std::endl;
  for (int i = 0; i < 64; i++)
  {
    std::cout << output.at(i) << std::endl;
  }
  std::cout << "================" << std::endl;
  for (int i = 0; i < 64; i++)
  {
    std::cout << real_inverse.at(i) / 64 << std::endl;
  }
}

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

void fft3D(double G_1D[], int npixels)
{
  uint64_t npixels3 = npixels * npixels * npixels;
  heffte::box3d<> inbox = {{0, 0, 0}, {npixels - 1, npixels - 1, npixels - 1}};
  heffte::box3d<> outbox = {{0, 0, 0}, {npixels - 1, npixels - 1, npixels - 1}};

  // define the heffte class and the input and output geometry
  heffte::fft3d<heffte::backend::fftw> fft(inbox, outbox, MPI_COMM_WORLD);

  // vectors with the correct sizes to store the input and output data
  // taking the size of the input and output boxes
  std::vector<double> input(fft.size_inbox());
  std::vector<std::complex<double>> output(fft.size_outbox());

  // fill the input vector with data that looks like 0, 1, 2, ...
  // std::iota(input.begin(), input.end(), 0); // put some data in the input
  for (uint64_t i = 0; i < npixels3; i++)
  {
    input.at(i) = G_1D[i];
  }

  // perform a forward DFT
  fft.forward(input.data(), output.data());

  for (uint64_t i = 0; i < npixels3; i++)
  {
    G_1D[i] = abs(output.at(i));
  }
}

void gridding(int npart, double Lbox, double xpos[], double ypos[], double zpos[], double v[], double ro[],
              double mass, double partperpixel, int npixels, double G_1D[])
{
  int64_t counts[npart] = {0};
  int64_t counts_reduced = 0;
  double dx, dx2, ydx, W, weight, h1, h2, h2_2, s1, s1_2, s2;
  double D[npixels], norm[npixels][npixels][npixels];
  double G[npixels][npixels][npixels];

  dx = Lbox / npixels;
  ydx = 1.e0 / dx;
  dx2 = dx / 2;
  // counts = 0;

  auto start = std::chrono::steady_clock::now();
#pragma omp parallel for
  for (int i = 0; i < npixels; i++)
  {
    D[i] = (2 * i + 1) * dx2;
  }

#pragma omp parallel for collapse(2)
  for (int i = 0; i < npixels; i++)
  {
    for (int j = 0; j < npixels; j++)
    {
      for (int k = 0; k < npixels; k++)
      {
        G[i][j][k] = 0;
        norm[i][j][k] = 0;
      }
    }
  }
  auto end = std::chrono::steady_clock::now();

  // std::cout << "Initialization took: "
  //       << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
  //       << " ms" << std::endl;
  // h1=std::cbrt(partperpixel/npart)*Lbox/4;
  // h1=std::sqrt(3./2.)/std::cbrt(partperpixel/npart)*Lbox/4;
  h1 = 1.003 * std::cbrt(3. / 4. / M_PI * partperpixel / npart) * Lbox / 2;
  h2 = 2 * h1;
  h2_2 = h2 * h2;

#pragma omp parallel for firstprivate (s1_2, s1, s2, weight, W)
  for (int n = 0; n < npart; n++)
  {
    weight = mass / ro[n];
    // if (n % 100000 == 0)
    // {
    //   std::cout << "particle " << n << " of " << npart << std::endl;
    // }

    int max_intz = std::floor((zpos[n] + h2) * ydx - 0.5e0);
    int min_intz = std::ceil((zpos[n] - h2) * ydx - 0.5e0);
    for (int i = min_intz; i <= max_intz; i++)
    {
      int zindex;
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
      double z = std::min(std::abs(zpos[n] - D[zindex]), std::abs(zpos[n] + Lbox - D[zindex]));
      z = z * z;
      if (z > h2_2)
        continue;

      s1_2 = h2_2 - z;
      s1 = std::sqrt(s1_2);
      int min_intx = std::ceil((xpos[n] - s1) * ydx - 0.5e0);
      int max_intx = std::floor((xpos[n] + s1) * ydx - 0.5e0);

      for (int j = min_intx; j <= max_intx; j++)
      {
        int xindex;
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
        double x = std::min(std::abs(xpos[n] - D[xindex]), std::abs(xpos[n] + Lbox - D[xindex]));
        x = x * x;
        if (x > s1_2)
          continue;
        s2 = std::sqrt(s1_2 - x);
        int min_inty = std::ceil((ypos[n] - s2) * ydx - 0.5e0);
        int max_inty = std::floor((ypos[n] + s2) * ydx - 0.5e0);

        for (int k = min_inty; k <= max_inty; k++)
        {
          int yindex;
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

          double y = std::min(std::abs(ypos[n] - D[yindex]), std::abs(ypos[n] + Lbox - D[yindex]));
          double r = std::sqrt(x + y * y + z);

          if (r < h2)
          {
            W = weight * cubickernel(r, h1);
            G[xindex][yindex][zindex] = G[xindex][yindex][zindex] + W * v[n];
            norm[xindex][yindex][zindex] = norm[xindex][yindex][zindex] + W;
            counts[n]++;
            // counts_reduced++;
          }
        }
      }
    }
  }


#pragma omp parallel for reduction(+:counts_reduced)
  for (int i = 0; i < npart; i++)
  {
    counts_reduced += counts[i];
  }
  std::cout << "counts: " << counts_reduced << std::endl;
  std::cout << "number of particles contributing to a pixel: " << counts_reduced * 1.0 / npixels / npixels / npixels << std::endl;
  start = std::chrono::steady_clock::now();

  double average = 0.0;

#pragma omp parallel for collapse(2)
  for (int i = 0; i < npixels; i++)
  {
    for (int j = 0; j < npixels; j++)
    {
      for (int k = 0; k < npixels; k++)
      {
        uint64_t iii = i + npixels * (j + npixels * k);
        if (norm[i][j][k] == 0)
        {
          std::cout << "Failed Rasterization in iteration i = " << i << ", j = " << j << ", k = " << k << std::endl;
          std::terminate();
        }
        G_1D[iii] = G[i][j][k] / norm[i][j][k];
        average += G_1D[iii]*G_1D[iii];
      }
    }
  }

  std::cout << "root mean square: " << sqrt((average/(npixels*npixels*npixels))) << std::endl;

  // fft3D(G_1D, npixels);
 
  end = std::chrono::steady_clock::now();

  std::cout << "Rasterization took: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        << " ms" << std::endl;
}

