#/usr/bin/env python3
# ----------------
#%matplotlib inline
from dnsdata import *
# MATPLOTLIB 3D
from skimage import measure
# VISPY
from vispy import app,scene
# ----------------

# Inputs
# ----------------
path = '../../'
nxc = nxd
nzc = nzd//2
# ----------------
rx = np.arange(0,nxc)*(pi/alfa0)/nxd
rz = np.arange(0,nzc)*(2*pi/beta0)/nzd

# Define type
FLUXVEC  = np.dtype([('xx', np.float64), ('yy', np.float64),   ('zz', np.float64)])
GKETERMS = np.dtype([('rTURB', FLUXVEC), ('rMEAN', FLUXVEC),   ('rVISC', FLUXVEC),
                     ('cTURB', 'float64'), ('cPRES', 'float64'), ('cVISC', 'float64'),
                     ('scaleENER', 'float64'), ('scalePROD', 'float64')])

# Load data
startpos = np.zeros(ny//2+3,dtype=int); startpos[1:] = np.cumsum(ny-2*np.arange(-1,ny//2+1)+1)
gke = np.memmap(path+'gke-complete.bin',dtype=GKETERMS,mode='r',shape=(startpos[-1],nxc,nzc))

# Extract the case ry=0
gkei = np.zeros([ny//2+2,nxc,nzc],dtype=GKETERMS)
for iy in range(0,ny//2+2): gkei[iy]=gke[startpos[iy]]

# Plot something
plt.contourf(rz,y[1:ny//2+2],gkei[1:,0,:]['scalePROD']-gkei[1:,0,:]['rMEAN']['yy'])
plt.colorbar()

# Plot the plane iz=16, ix=0
del gkei; gkei = np.zeros(startpos[-1],dtype=GKETERMS); gkei = gke[:,0,16]; Y=np.zeros([startpos[-1],2]);
for y1 in range(-1,ny//2+1):
    for y2 in range(y1,ny-y1+1):
        Y[startpos[iyd(y1)]+y2-y1,0]=0.5*(y[iyd(y1)]+y[iyd(y2)])
        Y[startpos[iyd(y1)]+y2-y1,1]=0.5*(y[iyd(y2)]-y[iyd(y1)])
grid_Yc, grid_ry = np.mgrid[0:1:200j, 0:1:200j]
gkeI = sp.interpolate.griddata(Y, gkei[...]['scaleENER'], (grid_Yc, grid_ry), method='linear')
plt.figure()
plt.contourf(grid_ry,grid_Yc,gkeI)
plt.show()
