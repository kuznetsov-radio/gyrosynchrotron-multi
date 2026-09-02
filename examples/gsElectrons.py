import numpy as np
import matplotlib.pyplot as plt
import GScodes # initialization library - located either in the current directory or in the system path

libname='./GSmulti_64.dll' # name of the executable library - located where Python can find it

GET_MW=GScodes.initGET_MW(libname) # load the library

Nf=1280          # number of frequencies
L=725e5*15       # source depth, cm
S=(725e5*20)**2  # source area, cm^2
 
Lparms=np.zeros(4, dtype='int32') # array of dimensions etc.
Lparms[0]=Nf
Lparms[1]=1  # number of components
Lparms[2]=32 # number of integration nodes
Lparms[3]=1  # free-free off
 
Rparms=np.zeros(4, dtype='double') # array of global floating-point parameters
Rparms[0]=S      # source area, cm^2
Rparms[1]=3e8    # starting frequency to calculate spectrum, Hz
Rparms[2]=0.0025 # logarithmic step in frequency
Rparms[3]=1e-4   # mu_0 computation accuracy
 
Sparms=np.zeros(5, dtype='double') # array of global source parameters
Sparms[0]=L   # voxel depth, cm
Sparms[1]=5e9 # n_0 - thermal electron density, cm^{-3}
Sparms[2]=2e7 # T_0, K
Sparms[3]=270 # B - magnetic field, G
Sparms[4]=75  # theta - viewing angle, /degrees
 
Bparms=np.zeros((14, 1), dtype='double', order='F') # array of voxel parameters
Bparms[0, 0]=1    # q/e
Bparms[1, 0]=1    # m/e
Bparms[2, 0]=1e6  # n_b
Bparms[3, 0]=3    # energy distribution (PLW)
Bparms[4, 0]=0.1  # E_min
Bparms[5, 0]=10   # E_max
Bparms[6, 0]=1.0  # E_break
Bparms[7, 0]=4    # \delta_1
Bparms[8, 0]=6    # \delta_2
Bparms[9, 0]=0    # pitch-angle distribution (ISO)
Bparms[10, 0]=60  # \alpha_c
Bparms[11, 0]=0.1 # \Delta\mu
 
Flux=np.zeros((3, Nf), dtype='double', order='F') # input/output array

# calculating the emission
res=GET_MW(Lparms, Rparms, Sparms, Bparms, Flux)
 
# retrieving the results
f=Flux[0]
IX=Flux[1]
IO=Flux[2]

# plotting the results
plt.figure(1)
plt.plot(f, IX)
plt.plot(f, IO)
plt.xscale('log')
plt.yscale('log')
plt.title('Emission intensity')
plt.xlabel('Frequency, GHz')
plt.ylabel('Intensity, sfu')

plt.show()
