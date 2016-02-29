#!/bin/python3

#
# Reads the Power Spectral Density (PSD) and Cross Spectral Density (CSD)
# produced by psd.cpl and plots the results
#
# The spectra and co-spectra are averaged in the two channel halves
#
# This program is NOT parallel
#


# Inputs
py = 89
iy = 0
Ret = 1000
scale = 'log'
path = './nx-1024-nz-1024/'

case  = ['00'     , '03']

#U     = [ 17.4555, 20.9478]    # U(y+=150)
#U     = [ 19.9961 ,  23.3137]   # Ub
U     = [  1.0, 1.0 ]           # No normalization

#M     = [ 1.08, 1.32]            # U(y+=150)
#M     = [ 1.18, 1.18]            # Ub normalization
M     = [ 27.58, 27.58 ]         # No noramalization

lineStyle = [['-', '--'],
	     ['solid', 'dashed']]
# ----------------
from   dnsdata import *
# ----------------

# Define type
izd = lambda i: i+nz
iyd = lambda i: i+1

# Write to screen what you are looking at
print('y(py) = ',y[iyd(py)],)
print('y(iy) = ',y[iyd(iy)],)

# Load data
for i in range(len(case)):
  clev = np.array([-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, 0.2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])*M[i]
  psd = np.load(path+case[i]+'/xpsd_iyP_'+str(py)+'.npy','c')
  
# Show minima and maxima
  print(np.max(psd[iyd(iy)])/U[i],np.min(psd[iyd(iy)])/U[i])
 
# Plot data
  cax = plt.contour(2*pi/kx[1:nz+1],2*pi/kz[1:nx+1],psd[iyd(iy),1:nx+1,1:nz+1]/U[i],clev,linestyles=lineStyle[1][i])
  #plt.gcf().colorbar(cax)
  plt.gca().set_yscale('log'); plt.gca().set_xscale('log')
  
# Save the data
  np.savetxt('xpsd_xz_'+case[i]+'.dat', psd[iyd(iy),1:nx+1,1:nz+1]/U[i])


plt.show()
