import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
filename="v_sphexa_100_18s_4.txt"
input = np.loadtxt(filename, dtype='d', delimiter=' ') #Load Rasterized image to perform FFT

npixel = round(len(input)**(1/3));

v=np.zeros((npixel,npixel,npixel));
for i in range(npixel):
    for j in range(npixel):
        for k in range(npixel):                     #Reconstruct 3D data from 1D file
            iii=i+npixel*((j-1)+npixel*(k-1));
            v[i,j,k]=input[iii];

n=int(npixel/2);
plt.figure(1)
plt.figure(figsize=(10,10))
plt.imshow(v[n,:,:],extent=[0,1,0,1],cmap='afmhot')
plt.title('Velocity Field (v) N = ' + str(int(npixel/4)) +'^3' ) #Image of transversal cut
plt.colorbar()

w=np.fft.fftn(v)/npixel**3; #FFT from numpy

Ew=2*np.pi*abs(w)**2; #Compute Energy

Ewx=np.sum(np.sum(Ew,axis=2),axis=1);
print(Ewx.shape)
Ewy=np.sum(np.sum(Ew,axis=2),axis=0);   #Perform averaging over te dfferent axis
print(Ewy.shape)
Ewz=np.sum(np.sum(Ew,axis=1),axis=0);
print(Ewz.shape)

Ew_mean=(Ewx*Ewy*Ewz)**(1/3);   # Geometrical Averaging
k=np.arange(npixel);
np.savetxt("PowerSpectra_"+filename,Ew_mean, delimiter=' ') #Save Spectra in file

## Plotting different Power Spectra

N= 6 #number of spectra to be plotted
k0=3 #Starting k value
Nyquist=2 #Arriving up to half the maximum k
plt.figure(3)
plt.figure(figsize=(15,10))
for i in range(N):
    E=np.loadtxt("PowerSpectra_v_sphexa_100_18s_"+ str(i+1) + ".txt")
    length=len(E)
    k=np.arange(length)
    k=k[k0:round(length/Nyquist)]; E=E[k0:round(length/Nyquist)];
    plt.loglog(k,k*E);

exp1=-5/3; exp2=-2;
plt.loglog(k,(k[0]**(-exp1)*E[0])*k**(exp1+1),'k',linestyle=(0,(5,10)));    #Kolmogorov -5/3 slope
plt.loglog(k,(k[0]**(-exp2)*E[0])*k**(exp2+1.),'k',linestyle=(0,(5,10)));   #-2 Slope

plt.xlabel('k'); plt.ylabel('k*E')
plt.title("Comparison Power Spectra")
