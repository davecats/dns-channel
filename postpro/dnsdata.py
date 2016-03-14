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

# Read dns.in
with open ('dns.in','r') as in_data:
   data=in_data.read().replace('\t','\n').replace(' ','\n')\
          .replace('yes','True').replace('=no','=False')\
          .replace("./Dati.cart.out","""str('./dati.cart.out')""")
dns_in={}; exec(data, globals(), dns_in);
nx = dns_in['nx']; nz = dns_in['nz']; ny = dns_in['ny']
alfa0=dns_in['alfa0']; beta0=dns_in['beta0']; a = dns_in['a']

# Index conversion
izd = lambda i: i+nz
iyd = lambda i: i+1

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
