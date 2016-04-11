#!/bin/python3
# ----------------
from  dnsdata           import *
from  scipy.interpolate import interp2d,interpn
# MATPLOTLIB 3D
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# VISPY
from vispy import app,scene
# ----------------

#
# Reads the data produced by gke.cpl and plots the results
#
# This program is NOT parallel but can use storedstructures to 
# reduce memory requirements.
#


# Inputs
# ----------------
path = './dat/'
averagey = True
isoplot = True
streamplot = False
U=1
# ----------------
rx = np.linspace(0,2.*pi/alfa0,2*nxd)
rz = np.linspace(-2*pi/beta0,2.*pi/beta0,nzd)
# --- helper function to plot streamlines ---
def plotstreamlines(x,y,u,v,density=1,nx=100,ny=100):
  xx,yy = np.linspace(x.min(),x.max(),nx), np.linspace(y.min(),y.max(),ny)
  gu,gv = interp2d(x,y,u), interp2d(x,y,v)
  plt.streamplot(xx,yy,gu(xx,yy),gv(xx,yy), density)
def eqspcdgrd(x,y,z,F,nx=100,ny=100,nz=100):
  xx,yy,zz = np.linspace(x.min(),x.max(),nx), np.linspace(y.min(),y.max(),ny), np.linspace(z.min(),z.max(),nz)
  FF = interpn((x,y,z),F,[(ix,iy,iz) for ix in xx for iy in yy for iz in zz]).reshape(nx,ny,nz)
  return xx,yy,zz,FF

# Define type
FLUXVEC   = np.dtype([('xx', np.float64), ('yy', np.float64), ('zz', np.float64)])
SCALEFLUX = np.dtype([('TURB', FLUXVEC), ('MEAN', FLUXVEC), ('VISC', FLUXVEC)])
SPACEFLUX = np.dtype([('TURB', np.float64), ('PRES', np.float64), ('VISC', np.float64)])

# Load already postprocessed data, if present
try:
    psd = np.load('XXXXXXX.npy','r')
    
# Otherwise read results
except IOError:
    # Read source/sink terms
    with open (path+'source.dat','r') as source_file:
      source = np.fromfile(source_file,np.float64,(ny+3)*(2*nxd)*(nzd)).reshape(ny+3,2*nxd,nzd)
      sink   = np.fromfile(source_file,np.float64,(ny+3)).reshape(ny+3,1,1)
    # Read space fluxes
    with open (path+'phiC.dat','r') as space_file:
      phiC = np.fromfile(space_file,SPACEFLUX,(ny+3)*(2*nxd)*(nzd)).reshape(ny+3,2*nxd,nzd)
    # Read scale fluxes
    with open (path+'phiR.dat','r') as scale_file:
      phiR = np.fromfile(scale_file,SCALEFLUX,(ny+3)*(2*nxd)*(nzd)).reshape(ny+3,2*nxd,nzd)
    # Average top and bottom plane
    if averagey:
      source[1:iyd(int(ny/2)),:,:] += source[ny+1:iyd(int(ny/2)):-1,:,:]
      sink[1:iyd(int(ny/2)),:,:] += sink[ny+1:iyd(int(ny/2)):-1,:,:]
      source[1:iyd(int(ny/2)),:,:] *= 0.5
      sink[1:iyd(int(ny/2)),:,:] *= 0.5
      for term in SPACEFLUX.names:
        phiC[1:iyd(int(ny/2)),:,:][term] -= phiC[ny+1:iyd(int(ny/2)):-1,:,:][term]
        phiC[1:iyd(int(ny/2)),:,:][term] *= 0.5
      for term in SCALEFLUX.names:
        for cord in FLUXVEC.names:
          phiR[1:iyd(int(ny/2)),:,:][term][cord] += (-1 if cord=='yy' else 1)*phiR[ny+1:iyd(int(ny/2)):-1,:,:][term][cord]
          phiR[1:iyd(int(ny/2)),:,:][term][cord] *= 0.5           
    # Remove top half
      source, sink = source[0:iyd(int(ny/2))+1], sink[0:iyd(int(ny/2))+1]
      phiC, phiR, y = phiC[0:iyd(int(ny/2))+1], phiR[0:iyd(int(ny/2))+1], y[0:iyd(ny/2)+1]
    # Compute total flux
    phi = np.zeros(np.shape(source),FLUXVEC)
    for term in SCALEFLUX.names:
      for cord in ['xx','zz']:
        phi[...][cord] += phiR[...][term][cord]
    for term in SPACEFLUX.names:
      phi[...]['yy'] += phiC[...][term]
    # Saving results
    

# Contour plot (rz,Yx)-space
plt.figure()
cax = plt.contourf(rz[0:nzd//2],y,(source+sink)[:,0,0:nzd//2])
Q = plt.quiver(rz[0:nzd//2],y,phi[:,0,0:nzd//2]['zz'], phi[:,0,0:nzd//2]['yy'])
plt.colorbar(cax)
if streamplot: plotstreamlines(rz[0:nzd//2],y,phi[:,0,0:nzd//2]['zz'],phi[:,0,0:nzd//2]['yy'],nx=300,ny=300,density=5)
plt.xlabel(r'$r_z$')
plt.ylabel(r'$Y_c$')


# Contour plot (rx,Yx)-space
plt.figure()
cax = plt.contourf(rx[0:nxd],y,(source+sink)[:,0:nxd,0])
Q = plt.quiver(rx[0:nxd],y,phi[:,0:nxd,0]['xx'], phi[:,0:nxd,0]['yy'])
plt.colorbar(cax)
if streamplot: plotstreamlines(rx[0:nxd],y,phi[:,0:nxd,0]['xx'],phi[:,0:nxd,0]['yy'],nx=300,ny=300,density=5)
plt.xlabel(r'$r_x$')
plt.ylabel(r'$Y_c$')


# Surface PLOTS
# ---------------------------
if isoplot:
  RX,RY,RZ,SOURCE = eqspcdgrd(y,rx[0:nxd],rz[0:nzd//2],(source+sink)[:,0:nxd,0:nzd//2],nx=200,ny=200,nz=200)
  print(np.shape(SOURCE))
  levels=[0.5*((source+sink).min()), 0.7*((source+sink).max())]

# Surface plot with Matplotlib
# ---------------------------
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for l in levels:
    verts, faces = measure.marching_cubes(SOURCE, l, spacing=(RX[1]-RX[0], RY[1]-RY[0], RZ[1]-RZ[0]))
    coll= ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
                          cmap='jet', lw=0, shade=True)
 
# Surface plot with VisPy
# --------------------------
  canvas=scene.SceneCanvas(keys='interactive',bgcolor='white')
  view=canvas.central_widget.add_view()
  for l in levels:
    surf = scene.visuals.Isosurface(SOURCE,level=l,shading='smooth',parent=view.scene)
    surf.transform = scene.transforms.STTransform(translate=(0, 0, 0),scale=(RX[1]-RX[0], RY[1]-RY[0], RZ[1]-RZ[0]))
    axis=scene.visuals.XYZAxis(parent=view.scene)
    cam=scene.TurntableCamera(elevation=30., azimuth=60.,distance=50.0)
    view.camera=cam


# Show all plots
# --------------------------
  canvas.show()
plt.show()
