#!/bin/python3
# ----------------
from  dnsdata           import *
from  scipy.interpolate import interp2d
# VISVIS
import visvis as vv
from visvis.utils.iso import isosurface, isocontour 
from visvis.utils import iso
# MATPLOTLIB 3D
from skimage import measure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
# VISPY
from vispy import app,scene
# ----------------

#
# Reads the data produced by kolmogorov.cpl and plots the results
#
# This program is NOT parallel
#


# Inputs
# ----------------
us = 2
averagey = True
averagez = False
isoplot = False
streamplot = True
U=1
# ----------------
nz_us = int(nz/us)
nx_us = int(nx/us)
deltaX = 2.*pi/(2*nx_us)
deltaZ = 2.*pi/(2*nz_us)
rx = np.linspace(0,2.*pi/alfa0,nx_us+1)
rz = np.linspace(-2*pi/beta0,2.*pi/beta0,2*nz_us+1)
izd = lambda i: i+nz_us
# --- helper function to plot streamlines ---
def plotstreamlines(x,y,u,v,density=1.0,nx=100,ny=100):
  xx,yy = np.linspace(x.min(),x.max(),nx), np.linspace(y.min(),y.max(),ny)
  gu,gv = interp2d(x,y,u), interp2d(x,y,v)
  plt.streamplot(xx,yy,gu(xx,yy),gv(xx,yy), density)

# Define type
FLUXVEC   = np.dtype([('xx', np.float64), ('yy', np.float64), ('zz', np.float64)])
SCALEFLUX = np.dtype([('TURB', FLUXVEC), ('MEAN', FLUXVEC), ('VISC', FLUXVEC)])
DIVSCALEFLUX = np.dtype([('TURB', np.float64), ('MEAN', np.float64), ('VISC', np.float64)])
SPACEFLUX = np.dtype([('TURB', np.float64), ('PRES', np.float64), ('VISC', np.float64)])

# Load already postprocessed data, if present
try:
    psd = np.load('XXXXXXX.npy','r')
    
# Otherwise read results
except IOError:
    # Read source/sink terms
    with open ('source.dat','r') as source_file:
      source = np.fromfile(source_file,np.float64,(ny+3)*(nx_us+1)*(2*nz_us+1)).reshape(ny+3,nx_us+1,2*nz_us+1)
      sink   = np.fromfile(source_file,np.float64,(ny+3)).reshape(ny+3,1,1)
    # Read space fluxes
    with open ('phiC.dat','r') as space_file:
      phiC = np.fromfile(space_file,SPACEFLUX,(ny+3)*(nx_us+1)*(2*nz_us+1)).reshape(ny+3,nx_us+1,2*nz_us+1)
    # Read scale fluxes
    with open ('phiR.dat','r') as scale_file:
      phiR = np.fromfile(scale_file,SCALEFLUX,(ny+3)*(nx_us+1)*(2*nz_us+1)).reshape(ny+3,nx_us+1,2*nz_us+1)
      #DIVphiR = np.fromfile(scale_file,DIVSCALEFLUX,(ny+3)*(nx_us+1)*(2*nz_us+1)).reshape(ny+3,nx_us+1,2*nz_us+1)
      DIVphiR=np.zeros([ny+3,nx_us+1,2*nz_us+1],DIVSCALEFLUX)
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
        DIVphiR[1:iyd(int(ny/2)),:,:][term] += DIVphiR[ny+1:iyd(int(ny/2)):-1,:,:][term] 
        DIVphiR[1:iyd(int(ny/2)),:,:][term] *= 0.5  
        for cord in FLUXVEC.names:
          phiR[1:iyd(int(ny/2)),:,:][term][cord] += (-1 if cord=='yy' else 1)*phiR[ny+1:iyd(int(ny/2)):-1,:,:][term][cord]
          phiR[1:iyd(int(ny/2)),:,:][term][cord] *= 0.5           
    # Average the two nz half-planes 
    if averagez:
      source += source[:,:,::-1]; source *= 0.5
      for term in SPACEFLUX.names:
        phiC[...][term] += phiC[...,::-1][term]
        phiC[...][term] *= 0.5
      for term in SCALEFLUX.names:
        DIVphiR[...][term] += DIVphiR[...,::-1][term]
        DIVphiR[...][term] *= 0.5
        for cord in FLUXVEC.names:
          phiR[...,izd(1):][term][cord] += (-1 if cord=='zz' else 1)*phiR[...,izd(-1)::-1][term][cord]
          phiR[...,izd(1):][term][cord] *= 0.5
    #phiC, phiR, DIVphiR = phiC[:,:,izd(0):izd(nz_us)+1], phiR[:,:,izd(0):izd(nz_us)+1], DIVphiR[:,:,izd(0):izd(nz_us)+1]
    #source = source[:,:,izd(0):izd(nz_us)+1]
    # Derive the space flux
    DIVphiC = np.zeros([ny+3,nx_us+1,nz_us+1],SPACEFLUX)
    for ix,iz in ( (ix,iz) for ix in range(0,nx_us+1) for iz in range(0,nz_us+1) ):
      for term in SPACEFLUX.names:
        DIVphiC[:,ix,iz][term] = deriv(d1,phiC[:,ix,iz][term])  
    # Remove top half
    source, sink = source[0:iyd(int(ny/2))+1], sink[0:iyd(int(ny/2))+1]
    phiC, phiR, DIVphiR, DIVphiC = phiC[0:iyd(int(ny/2))+1], phiR[0:iyd(int(ny/2))+1], DIVphiR[0:iyd(int(ny/2))+1], \
                                   DIVphiC[0:iyd(int(ny/2))+1]
    y=y[0:iyd(ny/2)+1]
    # Compute the residual
    tot = source.copy()+sink.copy()
    #for term in SPACEFLUX.names:
    #  tot += DIVphiC[...][term]
    #for term in DIVSCALEFLUX.names:
    #  tot += DIVphiR[...][term]
    # Total Flux
    #phi = np.zeros([int(ny/2)+2,nx_us+1,nz_us+1],FLUXVEC)
    phi = np.zeros([int(ny/2)+2,nx_us+1,2*nz_us+1],FLUXVEC)
    for term in SCALEFLUX.names:
      for cord in ['xx','zz']:
        phi[...][cord] += phiR[...][term][cord]
    for term in SPACEFLUX.names:
      phi[...]['yy'] += phiC[...][term]

# Contour plot (rz,Yx)-space
plt.figure()
cax = plt.contourf(rz,y,tot[:,0,:])
Q = plt.quiver(rz,y,phi[:,0,:]['zz'], phi[:,0,:]['yy'])
plt.colorbar(cax)
if streamplot: plotstreamlines(rz,y,phi[:,0,:]['zz'],phi[:,0,:]['yy'],nx=300,ny=300,density=5)
plt.xlabel(r'$r_z$')
plt.ylabel(r'$Y_c$')


# Contour plot (rx,Yx)-space
plt.figure()
cax = plt.contourf(tot[:,:,0])
Q = plt.quiver(phi[:,:,0]['xx'], phi[:,:,0]['yy'])
plt.colorbar(cax)
plt.xlabel(r'$r_x$')
plt.ylabel(r'$Y_c$')


# Surface PLOTS
# ---------------------------
if isoplot:
  levels=[-0.001, 0.003]

# Surface plot with VISVIS
# ---------------------------
  try:
    vv.figure(1); vv.clf()
    for l in levels:
      bm1 = isosurface(tot, l, useValues=True)
      m1=vv.mesh(bm1,colormap=vv.CM_JET); #t=vv.volshow(tot)
    vv.gca().SetLimits(rangeX=[0.,nx_us],rangeY=[0.,nz_us],rangeZ=[0.,100.])
  except RuntimeError:
    print("VisVis: No isosurface with given levels found!")

# Surface plot with Matplotlib
# ---------------------------
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  for l in levels:
    verts, faces = measure.marching_cubes(tot, l, spacing=(1, 1, 1))
    coll= ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
                          cmap='jet', lw=0, shade=True)
 
# Surface plot with VisPy
# --------------------------
  canvas=scene.SceneCanvas(keys='interactive')
  view=canvas.central_widget.add_view()
  for l in levels:
    surf = scene.visuals.Isosurface(tot,level=l,color=plt.cm.jet(l-tot.min()/tot.max()),shading='smooth',parent=view.scene)
    axis=scene.visuals.XYZAxis(parent=view.scene)
    cam=scene.TurntableCamera(elevation=30., azimuth=30.)
    view.camera=cam


# Show all plots
# --------------------------
  canvas.show()
plt.show()





    
