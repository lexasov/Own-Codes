#!/usr/bin/env python3

import h5py
import numpy as np

import sys


def printSteps(fname):
    """ Display contents of HDF5 file: step, iteration and time """
    ifile = h5py.File(fname, "r")
    print(fname, "contains the following steps:")
    print("hdf5 step number".rjust(15), "sph iteration".rjust(15), "time".rjust(15))
    for i in range(len(list(ifile["/"]))):
        h5step = ifile["Step#%d" % i]
        print("%5d".rjust(14) % i, "%5d".rjust(14) % h5step.attrs["step"][0],
              "%5f".rjust(14) % h5step.attrs["time"][0])


def readStep(fname, step):
    ifile = h5py.File(fname, "r")
    try:
        h5step = ifile["Step#%s" % step]
        return h5step
    except KeyError:
        print(fname, "step %s not found" % step)
        printSteps(fname)
        sys.exit(1)

def gridding(npart, chunkSize, x, y, z, vx, vy, vz, n):
    xmin = 0.0
    xmax = 1.0
    halfbox = 0.5 * (xmax-xmin)
    partperpixel = 8.0
    Lbox = xmax - xmin
    npixels = 2*n
    mass = 1.0/npart
    Kmax = np.ceil(npixels * 0.5 * np.sqrt(3))

    xpos = [i + halfbox for i in x]
    ypos = [i + halfbox for i in y]
    zpos = [i + halfbox for i in z]


def plotSlice(fname, step):
    """ Plot a 2D xy-cross section with particles e.g. abs(z) < 0.1, using density as color """

    h5step = readStep(fname, step)

    # x = np.array(h5step["x"])
    # y = np.array(h5step["y"])
    # z = np.array(h5step["z"])
    # h = np.array(h5step["h"])
    # vx = np.array(h5step["vx"])
    # vy = np.array(h5step["vy"])
    # vz = np.array(h5step["vz"])

    # rho = np.array(h5step["rho"])
    # divv = np.array(h5step["divv"])
    # curlv = np.array(h5step["curlv"])

    f = open("turb_3000_output2.txt", "w")

    npart = len(h5step["z"])
    n = int(np.cbrt(npart))
    chunkNo = 10
    chunkSize = int(npart/chunkNo)
    for ind in range(chunkNo):
        x = np.array(h5step["x"][ind*chunkSize:(ind+1)*chunkSize])
        y = np.array(h5step["y"][ind*chunkSize:(ind+1)*chunkSize])
        z = np.array(h5step["z"][ind*chunkSize:(ind+1)*chunkSize])
        # h = np.array(h5step["h"][ind*chunkSize:(ind+1)*chunkSize])
        vx = np.array(h5step["vx"][ind*chunkSize:(ind+1)*chunkSize])
        vy = np.array(h5step["vy"][ind*chunkSize:(ind+1)*chunkSize])
        vz = np.array(h5step["vz"][ind*chunkSize:(ind+1)*chunkSize])
        for index in range(chunkSize):
            if (abs(z[index]/h[index]) < 2):
                v = np.sqrt (vx[index]*vx[index] + vy[index]*vy[index] + vz[index]*vz[index])
                str = "%lf %lf %lf %lf\n" % (x[index], y[index], z[index], v)
                f.write(str)
        print("chunk = ", ind)

    f.close()


if __name__ == "__main__":
    # first cmdline argument: hdf5 file name to plot
    fname = sys.argv[1]

    # second cmdline argument: hdf5 step number to plot or print (-p) and exit
    step = sys.argv[2]
    if step == "-p":
        printSteps(fname)
        sys.exit(1)

    plotSlice(fname, step)
