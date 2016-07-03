#!/bin/python3

#
# This Python3 modules declares variables and functions
# useful to postprocess DNS data obtained with dns-channel CPL
# program.
#

# ----------
import numpy             as np
import matplotlib        as mp
mp.use("Qt4Agg")
import matplotlib.pyplot as plt
from   numpy   import pi
# ----------

# Define some data tiypes
VETA = np.dtype([('v', np.complex128), ('eta', np.complex128)])

# Read dns.in
with open ('dns.in','r') as in_data:
   data=in_data.read().replace('\t','\n').replace(' ','\n')\
          .replace('yes','True').replace('=no','=False')\
          .replace("./Dati.cart.out","""str('./dati.cart.out')""")
dns_in={}; exec(data, globals(), dns_in);
nx = int(dns_in['nx']); nz = int(dns_in['nz']); ny = dns_in['ny']
alfa0=dns_in['alfa0'];  beta0=dns_in['beta0'];  a = dns_in['a']

# Index conversion
izd = lambda i: i+nz
iyd = lambda i: i+1

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

# Define y-compact derivatives
d1 = np.zeros([ny-1,5],np.float64); matder=np.zeros([5,5],np.float64)
d2 = np.zeros([ny-1,5],np.float64); tnder=np.zeros([5],np.float64)
for iy in range(1,ny):
  for ir,ic in ( (ir,ic) for ir in range(0,5) for ic in range(0,5) ):
    matder[ir,ic] = (y[iyd(iy-2+ic)] - y[iyd(iy)])**(4.0-np.float64(ir))
  tnder*=0; tnder[3]=1;     d1[iy-1,:]=np.linalg.solve(matder,tnder)
  tnder*=0; tnder[2]=2;     d2[iy-1,:]=np.linalg.solve(matder,tnder)
def deriv(d,f):
  df = np.zeros([ny+3],np.float64)
  for iy in range(1,ny):
    df[iyd(iy)] = np.sum(d[iy-1,:]*f[iyd(iy)-2:iyd(iy)+3])
  return df

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
