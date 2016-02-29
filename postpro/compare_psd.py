#!/bin/python3

#
# Reads the Power Spectral Density (PSD) and Cross Spectral Density (CSD)
# produced by psd.cpl and plots the results
#
# The spectra and co-spectra are averaged in the two channel halves
#
# This program is NOT parallel
#


# Inputs
pl = 'uv'
premultiply = False
y_normalize = True
path = './nx-1024-nz-1024/'


case  = ['00'     , '03']
U     = [ 19.9961 ,  23.3137]
M     = [ 0.321, 0.321 ]
lineStyle = [['-', '--'],
	     ['solid', 'dashed']]
# ----------------
from   dnsdata import *
# ----------------



# Define type
SPECTRUM = np.dtype([('uu', np.float64), ('vv', np.float64), ('ww', np.float64),
                     ('uv', np.float64), ('uw', np.float64), ('vw', np.float64)])
izd = lambda i: i+nz
iyd = lambda i: i+1

# Define some figure properties
f, ax = plt.subplots(2, 3)

f2 = plt.figure(figsize=(8,8))
plt.set_cmap('gnuplot2')

# Load data
for i in range(len(case)):
  clev = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])*M[i]
  psd_x = np.load(path+case[i]+'/psd_x.npy','c')
  psd_z = np.load(path+case[i]+'/psd_z.npy','c')
  psd_1D_yx = np.load(path+case[i]+'/psd_1D_yx.npy','c')
  psd_1D_z = np.load(path+case[i]+'/psd_1D_z.npy','c')
  psd_1D_x = np.load(path+case[i]+'/psd_1D_x.npy','c')
  
# Compute CfT
  CfT=6*(np.sum(2*psd_1D_x[1:nx+1][pl])+psd_1D_x[0][pl])/np.power(U[i],2)
  
# Show minima and maxima
  print(np.max(6*2*psd_x[iyd(1):,1:nx+1][pl]/np.power(U[i],2))/CfT,np.max(6*psd_z[iyd(1):,1:nz+1][pl]/np.power(U[i],2))/CfT)
 
# Plot data
  plt.figure(f.number)
  # uv - y    
  ax[0,0].plot(6*psd_1D_yx[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U[i],2)/CfT,y[iyd(0):iyd(int(ny/2))+1],lineStyle[0][i])
  np.savetxt('uv_y_'+case[i]+'.dat', [6*psd_1D_yx[iyd(0):iyd(int(ny/2))+1][pl]/np.power(U[i],2)/CfT,y[iyd(0):iyd(int(ny/2))+1]])
  # uv - lambdax
  ax[0,1].contour(2*pi/kx[1:nx+1],y[iyd(1):iyd(int(ny/2))+1],
                  6*2*psd_x[iyd(1):,1:nx+1][pl]/np.power(U[i],2)/CfT,clev,linestyles=lineStyle[1][i])
  np.savetxt('uv_yx_'+case[i]+'.dat', 6*2*psd_x[iyd(1):,1:nx+1][pl]/np.power(U[i],2)/CfT)
  # uv - lambdaz
  ax[0,2].contour(2*pi/kz[1:nz+1],y[iyd(1):iyd(int(ny/2))+1],
                  6*psd_z[iyd(1):,1:nz+1][pl]/np.power(U[i],2)/CfT,clev,linestyles=lineStyle[1][i])
  np.savetxt('uv_yz_'+case[i]+'.dat', 6*psd_z[iyd(1):,1:nz+1][pl]/np.power(U[i],2)/CfT)
  # uv - x
  ax[1,1].plot(2*pi/kx[1:nx+1],2*6*psd_1D_x[1:nx+1][pl]/np.power(U[i],2)/CfT,lineStyle[0][i])
  print(6*(np.sum(2*psd_1D_x[1:nx+1][pl])+psd_1D_x[0][pl])/np.power(U[i],2))
  np.savetxt('uv_x_'+case[i]+'.dat', [2*pi/kx[1:nx+1],2*6*psd_1D_x[1:nx+1][pl]/np.power(U[i],2)]/CfT)
  # uv - z
  ax[1,2].plot(2*pi/kz[1:nx+1],6*psd_1D_z[1:nz+1][pl]/np.power(U[i],2)/CfT,lineStyle[0][i])
  print(np.sum(6*psd_1D_z[0:nz+1][pl]/np.power(U[i],2)))
  np.savetxt('uv_z_'+case[i]+'.dat', [2*pi/kz[1:nx+1],6*psd_1D_z[1:nz+1][pl]/np.power(U[i],2)]/CfT)
  
  plt.figure(f2.number)
  plt.contour(2*pi/kz[1:nz+1],y[iyd(1):iyd(int(ny/2))+1],
              6*psd_z[iyd(1):,1:nz+1][pl]/np.power(U[i],2)/CfT,clev,linestyles=lineStyle[1][i]) 
  

plt.show()
