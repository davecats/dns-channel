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
py = 15;
premultiply = True
y_normalize = True
scale = 'log'
antisymm = True
U=15.8956
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
    psd_1D_yx = np.load('psd_1D_yx.npy','r')
    psd_1D_yz = np.load('psd_1D_yz.npy','r')
    psd_1D_z = np.load('psd_1D_z.npy','r')
    psd_1D_x = np.load('psd_1D_x.npy','r')
    
# Otherwise read results
except IOError:
    with open ('psd.bin','r') as psd_file:
      psd = np.fromfile(psd_file,SPECTRUM)
    psd = psd.reshape(ny+3,nx+1,2*nz+1)

# Average top and bottom plane
    for uiuj in SPECTRUM.names:
      if antisymm:
        psd[1:iyd(int(ny/2)),:,:][uiuj] -= psd[ny+1:iyd(int(ny/2)):-1,:,:][uiuj]
      else:
        psd[1:iyd(int(ny/2)),:,:][uiuj] += psd[ny+1:iyd(int(ny/2)):-1,:,:][uiuj]
      psd[1:iyd(int(ny/2)),:,:][uiuj] *= 0.5*(-1 if antisymm else 1)*alfa0*beta0
    psd = psd[0:iyd(int(ny/2))+1]
    
# Sum the two nz half-planes
    for uiuj in SPECTRUM.names:
      psd[:,:,:][uiuj] += psd[:,:,::-1][uiuj]
      psd[:,:,izd(nz)][uiuj] *= 0.5
    psd = psd[:,:,izd(0):izd(nz)+1]
    
# Normalize (1-y/H)
    if y_normalize:
      for uiuj in SPECTRUM.names:
        psd[:,:,:][uiuj] *= (1-y[0:iyd(int(ny/2))+1]).reshape(int(ny/2)+2,1,1)

# Compute 1D spectra
    psd_x = np.zeros((int(ny/2)+2,nx+1),SPECTRUM)
    psd_z = np.zeros((int(ny/2)+2,nz+1),SPECTRUM)
    psd_y = np.zeros((nx+1,nz+1),SPECTRUM)
    psd_1D_x = np.zeros(nx+1,SPECTRUM)
    psd_1D_z = np.zeros(nz+1,SPECTRUM)
    psd_1D_yx= np.zeros((int(ny/2)+2),SPECTRUM)
    psd_1D_yz= np.zeros((int(ny/2)+2),SPECTRUM)
    for uiuj in SPECTRUM.names:
      psd_z[:,:][uiuj] = psd[:,0,:][uiuj]+2*np.sum(psd[:,1:,:][uiuj],1)
      psd_x[:,:][uiuj] = np.sum(psd[:,:,:][uiuj],2)
      psd_1D_yx[:][uiuj]= psd_x[:,0][uiuj]+2*np.sum(psd_x[:,1:][uiuj],1)
      psd_1D_yz[:][uiuj]= np.sum(psd_z[:,:][uiuj],1)
      psd_1D_z[:][uiuj] = np.trapz(psd_z[iyd(0):,:][uiuj],y[iyd(0):iyd(int(ny/2))+1],axis=0)
      psd_1D_x[:][uiuj] = np.trapz(psd_x[iyd(0):,:][uiuj],y[iyd(0):iyd(int(ny/2))+1],axis=0)
      psd_y[:,:][uiuj] = np.trapz(psd[iyd(0):,:,:][uiuj],y[iyd(0):iyd(int(ny/2))+1],axis=0)
      
# Premultiply spectra
    if premultiply:
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
    for outdata in ('psd','psd_x','psd_z','psd_y','psd_1D_yx','psd_1D_yz','psd_1D_z','psd_1D_x'):
      np.save(outdata,eval(outdata)) 
      
# uv - y
f, ax = plt.subplots(2, 3)      
ax[0,0].plot(psd_1D_yx[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U,2),y[iyd(0):iyd(int(ny/2))+1])
ax[0,0].plot(psd_1D_yz[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U,2),y[iyd(0):iyd(int(ny/2))+1])
ax[0,0].set_yscale(scale)
# uv - lambdax
ax[0,1].contourf(2*pi/kx[1:nx+1],y[iyd(1):iyd(int(ny/2))+1],psd_x[iyd(1):,1:nx+1][pl]/np.power(U,2))
ax[0,1].set_xscale(scale); ax[0,1].set_yscale(scale)
# uv - lambdaz
ax[0,2].contourf(2*pi/kz[1:nz+1],y[iyd(1):iyd(int(ny/2))+1],psd_z[iyd(1):,1:nz+1][pl]/np.power(U,2))
ax[0,2].set_xscale(scale); ax[0,1].set_yscale(scale)
# uv - x
ax[1,1].plot(2*pi/kx[1:nx+1],psd_1D_x[1:nx+1][pl]/np.power(U,2))
ax[1,1].set_xscale(scale)
# uv - z
ax[1,2].plot(2*pi/kz[1:nx+1],psd_1D_z[1:nx+1][pl]/np.power(U,2))
ax[1,2].set_xscale(scale)

print(6*(psd_1D_x[0][pl]+2*np.sum(psd_1D_x[1:nx+1][pl]))/np.power(U,2))
print(6*np.sum(psd_1D_z[0:nz+1][pl]/np.power(U,2)))
print(6*np.trapz(psd_1D_yx[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U,2),y[iyd(0):iyd(int(ny/2))+1]))
print(6*np.trapz(psd_1D_yz[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U,2),y[iyd(0):iyd(int(ny/2))+1]))

plt.show()      
