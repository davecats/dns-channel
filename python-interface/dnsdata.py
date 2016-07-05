#!/bin/python3

#
# This Python3 modules declares variables and functions
# useful to postprocess DNS data obtained with dns-channel CPL
# program.
#

# ----------
import numpy             as np
from   scipy.linalg      import solve,solve_banded
import matplotlib        as mp
mp.use("Qt4Agg")
import scipy             as sp
import pyfftw            as fftw
import matplotlib.pyplot as plt
from   numpy   import pi
# ----------

# Define some data tiypes
VETA = np.dtype([('v', np.complex128), ('eta', np.complex128)])
VELOCITY = np.dtype([('u', np.complex128), ('v', np.complex128), ('w', np.complex128)])
REALVELOCITY = np.dtype([('u', np.float64), ('v', np.float64), ('w', np.float64)])
DERIVS = np.dtype([('d0', np.float64), ('d1', np.float64), ('d2', np.float64), ('d4', np.float64)])

# Read dns.in
with open ('dns.in','r') as in_data:
   data=in_data.read().replace('\t','\n').replace(' ','\n')\
          .replace('yes','True').replace('=no','=False')\
          .replace("./Dati.cart.out","""str('./dati.cart.out')""")
dns_in={}; exec(data, globals(), dns_in);
nx = int(dns_in['nx']); nz = int(dns_in['nz']); ny = dns_in['ny']
alfa0=dns_in['alfa0'];  beta0=dns_in['beta0'];  a = dns_in['a']

# Index conversion
izd = lambda i: i+nz   # converts (-nz..nz)  to (0..2*nz)
iyd = lambda i: i+1    # converts (-1..ny+1) to (0..ny+2)

# Find nxd,nzd
nxd=3*nx//2-1; nzd=3*nz-1;
def fftfit(x):
  while not x & 1:
    x = x >> 1
  return (x==1 or x==3)
while not fftfit(nxd): nxd += 1
while not fftfit(nzd): nzd += 1

# Set grid
y = np.tanh(a*(2*np.arange(-1,ny+2)/ny-1))/np.tanh(a)+1
kx = alfa0*np.arange(0,nx+1); kz = beta0*np.arange(0,nz+1)

# Put data in banded-matrix form
def banded_form(A):
  m=np.shape(A); D=np.zeros([m[1],m[0]]); u=m[1]//2;
  for i in range(m[0]):
    for j in range(m[1]):
      J=j+i-u;
      if J>-1 and J<m[0]: D[u+i-J,J] = A[i,j]
  return D

# Define y-compact derivatives
der=np.zeros([ny+3,5],DERIVS); matder=np.zeros([5,5],np.float64); rhsder=np.zeros([5],np.float64); cder=np.zeros([5],np.float64)
for iy in range(2,ny+1):
  matder=np.fromfunction(lambda i,j: (y[iy-2+j]-y[iy])**(4.0-np.float64(i)), (5,5), dtype=int)
  rhsder*=0; rhsder[0]=np.float64(24); der[iy,:]['d4']=solve(matder,rhsder)
  matder=np.fromfunction(lambda i,j: (5.0-i)*(6.0-i)*(7.0-i)*(8.0-i)*((y[iy-2+j]-y[iy])**(4.0-np.float64(i))), (5,5), dtype=int)
  rhsder*=0
  for i in range(5): rhsder[i]=np.sum(der[iy,:]['d4']*(y[iy-2:iy+3]-y[iy])**np.float64(8-i))
  der[iy,:]['d0']=solve(matder,rhsder);
  matder=np.fromfunction(lambda i,j: (y[iy-2+j]-y[iy])**(4-np.float64(i)), (5,5), dtype=int)
  rhsder*=0
  for i in range(3): rhsder[i]=np.sum((4.0-i)*(3.0-i)*der[iy,:]['d0']*(y[iy-2:iy+3]-y[iy])**np.float64(2-i))
  der[iy,:]['d2']=solve(matder,rhsder)
  rhsder*=0
  for i in range(4): rhsder[i]=np.sum((4.0-i)*der[iy,:]['d0']*(y[iy-2:iy+3]-y[iy])**np.float64(3-i))
  der[iy,:]['d1']=solve(matder,rhsder)
  # Bottom wall
der[1,1]['d0']=np.float64(1.0);
matder=np.fromfunction(lambda i,j: (y[0+j]-y[1])**(4-np.float64(i)), (5,5), dtype=int)
rhsder=np.zeros([5],np.float64); rhsder[3]=np.float64(1.0); der[1,:]['d1']=solve(matder,rhsder)
rhsder=np.zeros([5],np.float64); rhsder[2]=np.float64(2.0); der[1,:]['d2']=solve(matder,rhsder)
matder=np.fromfunction(lambda i,j: (y[0+j]-y[0])**(4-np.float64(i)), (5,5), dtype=int)
rhsder=np.zeros([5],np.float64); rhsder[3]=np.float64(1.0); der[0,:]['d1']=solve(matder,rhsder)
rhsder=np.zeros([5],np.float64); rhsder[2]=np.float64(2.0); der[0,:]['d2']=solve(matder,rhsder)
  # Top wall
der[ny+1,3]['d0']=np.float64(1.0);
matder=np.fromfunction(lambda i,j: (y[ny-2+j]-y[ny+1])**(4-np.float64(i)), (5,5), dtype=int)
rhsder=np.zeros([5],np.float64); rhsder[3]=np.float64(1.0); der[ny+1,:]['d1']=solve(matder,rhsder)
rhsder=np.zeros([5],np.float64); rhsder[2]=np.float64(2.0); der[ny+1,:]['d2']=solve(matder,rhsder)
matder=np.fromfunction(lambda i,j: (y[ny-2+j]-y[ny+2])**(4-np.float64(i)), (5,5), dtype=int)
rhsder=np.zeros([5],np.float64); rhsder[3]=np.float64(1.0); der[ny+2,:]['d1']=solve(matder,rhsder)
rhsder=np.zeros([5],np.float64); rhsder[2]=np.float64(2.0); der[ny+2,:]['d2']=solve(matder,rhsder)
d0mat=banded_form(der[2:ny+1,:]['d0'])

# Derivate
def deriv(d,f,df):
  bcast = np.ones(len(np.shape(f)),int); bcast[0]=5
  for iy in range(0,ny+3):
    iyc = iy*(iy>1 and iy < ny+1) + ny*(iy>=ny+1) + 2*(iy<=1)
    df[iy,...] += np.sum(der[iy,:][d].reshape(bcast)*f[iyc-2:iyc+3,...],axis=0)
  df[2]-=der[2,1]['d0']*df[1]+der[2,0]['d0']*df[0];             df[3]-=der[3,0]['d0']*df[1]
  df[ny]-=der[ny,3]['d0']*df[ny+1]+der[ny,4]['d0']*df[ny+2];    df[ny-1]-=der[ny-1,4]['d0']*df[ny+1]
  df[2:ny+1]=solve_banded((2,2),d0mat,df[2:ny+1])

# Read a field file
def readfield(fname):
  DNSINFO = np.dtype([('ny', np.int32), ('nx', np.int32), ('nz', np.int32), ('dummy', np.int32),
                      ('time', np.float64), ('ymin', np.float64), ('ymax', np.float64),
                      ('a', np.float64), ('alfa0', np.float64), ('beta0', np.float64),
                      ('ni', np.float64)])
  with open(fname,'r') as diskfield:
    dns=np.fromfile(diskfield,DNSINFO,1)
  U=np.memmap(fname,dtype=np.float64,mode='r',shape=(ny+3),offset=(8*9))
  W=np.memmap(fname,dtype=np.float64,mode='r',shape=(ny+3),offset=(8*9+8*(ny+3)))
  V=np.memmap(fname,dtype=VETA,mode='r',shape=(ny+3,nx+1,2*nz+1),offset=(8*9+8*(ny+3)*2))
  return dns,U,W,V

# Read scalar field file
def readscalarfield(fname):
  return np.memmap(fname,dtype=np.complex128,mode='r',shape=(ny+3,nx+1,2*nz+1),offset=(8*9+8*(ny+3)*2+16*2*((ny+3)*(nx+1)*(2*nz+1))))

# Define a velocity class, with methods
class velocity:
    cdata=fftw.empty_aligned((3,ny+3,nx+1,2*nz+1),dtype='complex128')
    rdata=cdata.view(dtype=np.float64).reshape([3,ny+3,2*nx+2,2*nz+1])
    u=cdata[0,...].view(); v=cdata[1,...].view(); w=cdata[2,...].view()
    ru=cdata[0,...].view(dtype=np.float64).reshape([ny+3,2*nx+2,2*nz+1]);
    rv=cdata[1,...].view(dtype=np.float64).reshape([ny+3,2*nx+2,2*nz+1]);
    rw=cdata[2,...].view(dtype=np.float64).reshape([ny+3,2*nx+2,2*nz+1]);
    fftin=fftw.empty_aligned((nx+2,2*nz+1),dtype='complex128')
    fftout=fftw.empty_aligned((2*nx+2,2*nz+1),dtype='float64')
    fft=fftw.FFTW(fftin,fftout, axes=(-1,-2), direction='FFTW_BACKWARD')
    def ift(self): # consider weather adding dealiasing
      for iV,iy in [ (iV,iy) for iV in range(3) for iy in range(ny+3)]:
        self.fftin[:nx+1,0:nz+1]=self.cdata[iV,iy,:,nz:]; self.fftin[:nx+1,nz+1:]=self.cdata[iV,iy,:,0:nz]; self.fftin[nx+1,:]*=0
        self.fft.execute()
        self.rdata[iV,iy,...]=self.fftout
    def reset(self):
      self.cdata*=0

# Convert veta2uvw, Remember: eta=+I*beta*u-I*alfa*
def veta2uvw(veta,uvw):
  alfa=alfa0*np.arange(nx+1).reshape(1,nx+1,1); beta=beta0*np.arange(-nz,nz+1).reshape(1,1,2*nz+1); k2=alfa*alfa+beta*beta; k2[0,0,izd(0)]=1.0
  uvw.u+=veta[...]['eta']; uvw.u*=1j/k2; uvw.v+=veta['v']; deriv('d1',uvw.v,uvw.w); uvw.w*=1j/k2;
  for iy in range(ny+3):  vy=uvw.w[iy,...].copy(); uvw.w[iy,...]=beta[0,...]*vy+alfa[0,...]*uvw.u[iy,...]; uvw.u[iy,...]*=-beta[0,...]; uvw.u[iy,...]+=alfa[0,...]*vy
  uvw.u[:,0,izd(0)]*=0; uvw.v[:,0,izd(0)]*=0; uvw.w[:,0,izd(0)]*=0
