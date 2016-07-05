#/usr/bin/env python3
%matplotlib inline
from dnsdata import *
from matplotlib.pyplot import pcolormesh

ref = [1001,2001,3001,4001,5001]
ntot = 5

# Initialize velocity
v=velocity()
# Initialize variables
mean  = np.zeros([4,ny+3])
rms = np.zeros([6,ny+3])

for i in ref:
    # Reset velocity field
    v.reset()
    # Read field
    fname = './Field'+str(i)+'.fld'
    dns,U,W,V = readfield(fname)
    # Convert from (v,eta) to (u,v,w)
    veta2uvw(V,v)
    # Compute
    mean[0,:]+=U; mean[1,:]+=W; mean[2,:]+=deriv('d1',U); mean[3,:]+=deriv('d1',W)
    # Compute rms
    for iy in range(ny+3):
        rms[0,iy]+=np.real(2.0*np.vdot(v.u[iy,1:nx+1,:],v.u[iy,1:nx+1,:])+np.vdot(v.u[iy,0,1:],v.u[iy,0,1:]))
        rms[3,iy]+=np.real(2.0*np.vdot(v.u[iy,1:nx+1,:],v.v[iy,1:nx+1,:])+np.vdot(v.u[iy,0,1:],v.v[iy,0,1:]))
        rms[4,iy]+=np.real(2.0*np.vdot(v.u[iy,1:nx+1,:],v.w[iy,1:nx+1,:])+np.vdot(v.u[iy,0,1:],v.w[iy,0,1:]))
        rms[1,iy]+=np.real(2.0*np.vdot(v.v[iy,1:nx+1,:],v.v[iy,1:nx+1,:])+np.vdot(v.v[iy,0,1:],v.v[iy,0,1:]))
        rms[5,iy]+=np.real(2.0*np.vdot(v.v[iy,1:nx+1,:],v.w[iy,1:nx+1,:])+np.vdot(v.v[iy,0,1:],v.w[iy,0,1:]))
        rms[2,iy]+=np.real(2.0*np.vdot(v.w[iy,1:nx+1,:],v.w[iy,1:nx+1,:])+np.vdot(v.w[iy,0,1:],v.w[iy,0,1:]))
# Average
mean*=1/ntot; rms*=1/ntot
# Compute inner units
ut = np.sqrt(0.5*(mean[2,1]-mean[2,-2])/dns['ni'])
Retau = ut*dns['ni']
