#!/bin/python3

#
# Reads the Power Spectral Density (PSD) and Cross Spectral Density (CSD)
# produced by psd.cpl and plots the results
#
# The spectra and co-spectra are averaged in the two channel halves and in the two
# halves (-nz..-1) and (1..nz) 
#
# This program is NOT parallel
#

# Inputs
pl = 'uu'
py = 10;
utau = 1;
nu = 200;
# ----------------
from   dnsdata import *
# ----------------

iz = lambda i: i+nz

# Define type
SPECTRUM = np.dtype([('uu', np.float64), ('vv', np.float64), ('ww', np.float64),
                     ('uv', np.float64), ('uw', np.float64), ('vw', np.float64)])

# Read spectra
with open ('psd.bin','r') as psd_file:
  psd = np.fromfile(psd_file,SPECTRUM)
psd = np.reshape(psd,(ny+3,nx+1,2*nz+1))
psd = psd[1:ny+2] # throw away iy=-1 and iy=ny+1

# Average top and bottom plane
for iy in range(int(ny/2)):
  for uiuj in SPECTRUM.names:
    psd[iy,:,:][uiuj] = 0.5*(psd[iy,:,:][uiuj]+psd[ny-iy,:,:][uiuj])
psd = psd[0:int(ny/2)+1,:,:]    
    
# Average the two nz half-planes
for i in range(-nz,0):
  for uiuj in SPECTRUM.names:
    psd[:,:,iz(nz)-iz(i)][uiuj] = 0.5*(psd[:,:,iz(i)][uiuj]+psd[:,:,iz(nz)-iz(i)][uiuj])
psd = psd[:,:,iz(0):iz(nz)+1]
  
# Compute 1D spectra
psd_x = np.zeros((int(ny/2)+1,nx+1),SPECTRUM)
psd_z = np.zeros((int(ny/2)+1,nz+1),SPECTRUM)
for uiuj in SPECTRUM.names:
  psd_x[:,:][uiuj] = beta0*np.sum(psd[:,:,:][uiuj],2)
  psd_z[:,:][uiuj] = alfa0*np.sum(psd[:,:,:][uiuj],1)

# Premultiply spectra
for uiuj in SPECTRUM.names:
  for ix in range(0,nx+1):
    psd_x[:,ix][uiuj] = kx[ix]*psd_x[:,ix][uiuj]
    psd[:,ix,:][uiuj] = kx[ix]*psd[:,ix,:][uiuj]
  for iz in range(0,nz+1):
    psd_z[:,iz][uiuj] = kz[iz]*psd_z[:,iz][uiuj]
    psd[:,:,iz][uiuj] = kz[iz]*psd[:,:,iz][uiuj]
      
# Plot 2D spectra
plt.figure()
plt.contourf(2*pi/kz[1:nz+1],2*pi/kx[1:nx+1],psd[py,1:nx+1,1:nz+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.xlabel('$k_x h$'); plt.ylabel('$k_z h$'); plt.title(pl+' power spectrum at height '+'{:.2}'.format(y[py])+' $y/h$')

# Plot 1D spectra
plt.figure()
plt.contourf(2*pi/kx[1:nx+1],y[1:int(ny/2)+1],psd_x[1:int(ny/2)+1,1:nx+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.xlabel('$y/ h$'); plt.ylabel('$k_x h$'); plt.title(pl+' 1D power spectrum')

plt.figure()
plt.contourf(2*pi/kz[1:nz+1],y[1:int(ny/2)+1],psd_z[1:int(ny/2)+1,1:nz+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.xlabel('$y/ h$'); plt.ylabel('$k_z h$'); plt.title(pl+' 1D power spectrum')
plt.show()