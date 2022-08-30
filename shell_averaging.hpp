#include <vector>
#include <complex>
#include <iostream>

void shells(double w[], int npixels, double E[], double k_center[])
{
  int halfnpixels = std::floor(npixels / 2);
  int Kmax = std::ceil(std::sqrt(3.0) * (0.5 * npixels));
  // double E[Kmax],k_center[Kmax];
  int Kx, Ky, Kz, K, iii;
  int npixels3 = npixels * npixels * npixels;
  double y_N3 = 1. / npixels3;
  int counts[Kmax];

  for (int i = 0; i < Kmax; i++)
  {
    E[i] = 0;
    counts[i] = 0;
    k_center[i] = 2 * M_PI * std::cbrt(0.5 * ((i + 1) * (i + 1) * (i + 1) + i * i * i));
  }

  for (int i = 0; i < npixels; i++)
  {
    Kx = i;
    if (i > halfnpixels)
    {
      Kx = npixels - i;
    }
    std::cout << i << " of " << npixels << std::endl;
    for (int j = 0; j < npixels; j++)
    {
      Ky = j;
      if (j > halfnpixels)
      {
        Ky = npixels - j;
      }
      for (int k = 0; k < npixels; k++)
      {
        Kz = k;
        if (k > halfnpixels)
        {
          Kz = npixels - k;
        }
        iii = i + npixels * (j + npixels * k);
        K = std::floor(std::sqrt(Kx * Kx + Ky * Ky + Kz * Kz));

        // std::cout << K << std::endl; getchar();
        E[K] = E[K] + w[iii] * w[iii];

        counts[K]++;
      }
    }
  }

  for (int i = 0; i < Kmax; i++)
  {
    E[i] = 4 * M_PI * y_N3 * k_center[i] * k_center[i] * E[i] / counts[i];
  }
}
