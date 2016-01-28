#!/bin/python3

#
# This Python3 modules declares variables and functions
# useful to postprocess DNS data obtained with dns-channel CPL
# program. 
#

# ----------
import numpy             as np
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

# Set grid
y = np.tanh(a*(2*np.arange(ny+1)/ny-1))/np.tanh(a)+1
kx = alfa0*np.arange(0,nx+1); kz = beta0*np.arange(0,nz+1)