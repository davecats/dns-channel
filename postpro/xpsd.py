#!/bin/python3

#
# Reads the Cross Spectral Density (CSD)
# produced by xpsd.cpl and plots the results
#
# The co-spectra are averaged in the two channel halves and in the two
# halves (-nz..-1) and (1..nz) 
#
# This program is NOT parallel
#

# Inputs
py = 89
iy = 0
premultiply = True
y_normalize = False
scale = 'log'
U=1
# ----------------
from   dnsdata import *
# ----------------

#clevels = np.array([-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2,  0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95])*18.93

izd = lambda i: i+nz
iyd = lambda i: i+1

# Define type
SPECTRUM = np.dtype([('uu', np.float64), ('vv', np.float64), ('ww', np.float64),
                     ('uv', np.float64), ('uw', np.float64), ('vw', np.float64)])

# Load postprocessed PSD, if present
try:
    psd = np.load('xpsd_iyP_'+str(py)+'.npy','r')
    
# Otherwise read results
except IOError:
    with open ('xpsd_iyP_'+str(py)+'.bin','r') as psd_file:
      psd = np.fromfile(psd_file,np.float64)
    psd = psd.reshape(py+2,nx+1,2*nz+1)
    
# Sum the two nz half-planes
    psd[:,:,:] += psd[:,:,::-1]
    psd[:,:,izd(nz)] *= 0.5
    psd = psd[:,:,izd(0):izd(nz)+1]
    
# Normalize (1-y/H)
    if y_normalize:
      psd[:,:,:][uiuj] *= (1-y[0:iyd(py)+1]).reshape(py+2,1,1)
      
# Premultiply spectra
    if premultiply:
      for ix in range(0,nx+1):
        psd[:,ix,:] = kx[ix]*psd[:,ix,:]
      for iz in range(0,nz+1):
        psd[:,:,iz] = kz[iz]*psd[:,:,iz]

# Save the spectra
    np.save('xpsd_iyP_'+str(py)+'.npy',psd) 
      
print(np.max(psd[iyd(iy)]))

cax = plt.contourf(2*pi/kx[1:nz+1],2*pi/kz[1:nx+1],psd[iyd(iy),1:nx+1,1:nz+1])
plt.gcf().colorbar(cax)
plt.gca().set_yscale(scale); plt.gca().set_xscale(scale)

plt.show()      
