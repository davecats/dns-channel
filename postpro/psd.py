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
pl = 'uv'
py = 10;
utau = 1;
nu = 200;
# ----------------
from   dnsdata import *
# ----------------

izd = lambda i: i+nz
iyd = lambda i: i+1

# Define type
SPECTRUM = np.dtype([('uu', np.float64), ('vv', np.float64), ('ww', np.float64),
                     ('uv', np.float64), ('uw', np.float64), ('vw', np.float64)])

# Load postprocessed PSD, if present
try:
    psd = np.load('psd.npy','r')
    psd_x = np.load('psd_x.npy','r')
    psd_z = np.load('psd_z.npy','r')
    psd_y = np.load('psd_y.npy','r')
# Otherwise read results
except IOError:
    with open ('psd.bin','r') as psd_file:
      psd = np.fromfile(psd_file,SPECTRUM)
    psd = np.reshape(psd,(ny+3,nx+1,2*nz+1))
    #psd = psd[1:ny+2] # throw away iy=-1 and iy=ny+1

# Average top and bottom plane
    for iy in range(iyd(int(ny/2))):
      for uiuj in SPECTRUM.names:
        psd[iy,:,:][uiuj] = 0.5*(psd[iy,:,:][uiuj]+psd[iyd(ny-iy),:,:][uiuj])
    psd = psd[0:iyd(int(ny/2))+1]    
    
# Average the two nz half-planes
    for iz in range(-nz,0):
      for uiuj in SPECTRUM.names:
        psd[:,:,izd(nz)-izd(iz)][uiuj] = 0.5*(psd[:,:,izd(iz)][uiuj]+psd[:,:,izd(nz)-izd(iz)][uiuj])
    psd = psd[:,:,izd(0):izd(nz)+1]
  
# Compute 1D spectra
    psd_x = np.zeros((int(ny/2)+2,nx+1),SPECTRUM)
    psd_z = np.zeros((int(ny/2)+2,nz+1),SPECTRUM)
    psd_y = np.zeros((nx+1,nz+1),SPECTRUM)
    for uiuj in SPECTRUM.names:
      psd_x[:,:][uiuj] = beta0*np.sum(psd[:,:,:][uiuj],2)
      psd_z[:,:][uiuj] = alfa0*np.sum(psd[:,:,:][uiuj],1)
      for ix in range(0,nx+1):
        for iz in range(0,nz+1):
          psd_y[ix,iz][uiuj] = np.trapz(psd[iyd(0):,ix,iz][uiuj],y[iyd(0):iyd(int(ny/2))+1])

# Premultiply spectra
    for uiuj in SPECTRUM.names:
      for ix in range(0,nx+1):
        psd_x[:,ix][uiuj] = kx[ix]*psd_x[:,ix][uiuj]
        psd[:,ix,:][uiuj] = kx[ix]*psd[:,ix,:][uiuj]
        psd_y[ix,:][uiuj] = kx[ix]*psd_y[ix,:][uiuj]
      for iz in range(0,nz+1):
        psd_z[:,iz][uiuj] = kz[iz]*psd_z[:,iz][uiuj]
        psd[:,:,iz][uiuj] = kz[iz]*psd[:,:,iz][uiuj]
        psd_y[:,iz][uiuj] = kz[iz]*psd_y[:,iz][uiuj]

# Save the spectra
    for outdata in ('psd','psd_x','psd_z','psd_y'):
      np.save(outdata,eval(outdata)) 
      
# Plot 2D spectra
plt.figure()
plt.contourf(2*pi/kz[1:nz+1],2*pi/kx[1:nx+1],psd[py,1:nx+1,1:nz+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.xlabel('$k_x h$'); plt.ylabel('$k_z h$'); plt.title(pl+' power spectrum at height '+'{:.2}'.format(y[py])+' $y/h$')

# Plot 2D spectra (y-integrated)
plt.figure()
plt.contourf(2*pi/kz[1:nz+1]*1000,2*pi/kx[1:nx+1]*1000,psd_y[1:nx+1,1:nz+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.xlabel('$k_x h$'); plt.ylabel('$k_z h$'); plt.title('Power spectrum of the y-integral of '+pl)

# Plot 1D spectra
plt.figure()
plt.contourf(2*pi/kx[1:nx+1],y[iyd(1):iyd(int(ny/2))+1],psd_x[iyd(1):,1:nx+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.ylabel('$y/ h$'); plt.xlabel('$k_x h$'); plt.title(pl+' 1D power spectrum')

plt.figure()
plt.contourf(2*pi/kz[1:nz+1],y[iyd(1):iyd(int(ny/2))+1],psd_z[iyd(1):,1:nz+1][pl])
plt.xscale('log'); plt.yscale('log')
plt.ylabel('$y/ h$'); plt.xlabel('$k_z h$'); plt.title(pl+' 1D power spectrum')
plt.show()
