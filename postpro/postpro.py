#--------------------------------------------------#
#   						   #
#              Elaborate DNS results               #
#						   #
#--------------------------------------------------#
#
#
#  Author: Davide Gatti
#  Rev   : 0.2 20/Jul/2015
#  Email : davide.gatti@kit.edu
#

#!========================================
#! Input data
#!========================================
CUT         = 0;            # Minimum amount of initial transient to be cut in time units
correct     = 'False';         # Clean data files from unreadable lines? (Not implemented yet)
elaborate   = 'True';	    # Compute uncertainties
istransient = 'True'          # Cut the transient?
path       = './'
control    = 'TW'
fname_run  = ['Run', 'Runtimedata']
fname_pow  = ['Pow', 'Powerdata']
#!========================================
from re                   import findall
from sys     		  import stdout, exit
from os.path 		  import isfile
from os                   import walk
from numpy   		  import loadtxt, ones, zeros
from numpy   		  import concatenate, arange, argmin, mean
from numpy   		  import sign, power, sqrt, round, linspace
from numpy   	          import std, amax, amin, transpose, sum, subtract
from numpy   	          import where, empty
from scipy.interpolate    import interp1d
from scipy.signal         import lfilter, lfiltic
from scipy.constants      import pi
from ar                   import arsel
from hashlib              import md5
from itertools            import count
from openpyxl.compat      import range
from openpyxl             import load_workbook     
import matplotlib.pyplot  as plt
import openpyxl
#=========================================

# Read Runtimedata and Powerdata
for fname in fname_run:
  try: 
    data = loadtxt(path+fname)
    break
  except: 
    continue
else: 
  exit("ERROR: Neither Runtimedata nor Run found in path")

iscontrolled=True
for fname in fname_pow:
  try: 
    pdata = loadtxt(path+fname); 
    break
  except: continue
else: 
  iscontrolled=False

# Read dns.in
with open (path+'dns.in','r') as in_data:
   dns_in=in_data.read().replace('\t','\n').replace(' ','\n')\
          .replace('yes','True').replace('=no','=False')\
          .replace("./Dati.cart.out","""str('./dati.cart.out')""")
d1={}; exec(dns_in, globals(), d1); Re=d1["ni"];

# Clean Runtimedata from false lines if re-run
indx = where(data[:,1]==0)[0]
while indx.size>0:
  if indx[0]==data.shape[0]-1:
    data = data[0:indx[0]-2,:];
  elif indx[0]==0:
    data = data[1:data.shape[0]-1,:]
  else:
    data = concatenate((data[0:indx[0]-1,:],data[indx[0]+1:data.shape[0]-1,:]))
  indx = where(data[:,1]==0)[0]

k = 1; [m,n] = data.shape
while k<m-1:         
    if (data[k,0]<=data[k-1,0]):
        q = k+1;
        while (data[q,0]<data[k-1,0]) and (q<m-2):
            q = q+1;
        data = concatenate((data[0:k-1,:],data[q+1:m-1,:]))
        k = 1;
        [m,n] = data.shape
    k = k + 1;

if iscontrolled:
  k = 1; [h,l] = pdata.shape
  while k<h-1:         
      if pdata[k,0]<pdata[k-1,0]:
          q = k+1;
          while (pdata[q,0]<pdata[k-1,0]) and (q<h-2):
              q = q+1;
          pdata = concatenate((pdata[0:k-1,:],pdata[q+1:m-1,:]))
          k = 1;
          [h,l] = pdata.shape
      k = k + 1;

# Cutting the initial transient
CPI=False
if std(data[:,5])<1e-6:  # CFR
  dat = 0.5*(data[:,1]+data[:,2])
  print('CFR')
else:                    # CPI or CPG
  dat = data[:,5]
  print('CPI/CPG'); CPI=True

t = 3 if istransient else 1
s=empty([t,2],dtype=float); k=0; indx=3; ds_o = 0; d_max=0;
while (indx<m-2) and (k<t):
    indx=indx+5; sp1=std(dat[indx+1:m-1]); sm1=std(dat[indx:m-1]); 
    d=abs(dat[indx]-mean(dat[indx:m-1]))
    ds_n=sp1-sm1; d_max=(d_max if d<d_max else d)
    if (ds_n*ds_o) < 0:
      s[k,:] = [indx, sm1]; k=k+1
    ds_o=ds_n
indx=(s[argmin(s[:,1]),0] if d_max>6*amin(s[:,1]) else 1 )

if data[indx,0]<CUT:
    indx = where(data[:,0]>CUT)[0][0]

# Computing averages
print('\n\n  Analysis Results: ')
print('  ================================== ')
Tsim   = data[data.shape[0]-1,0] - data[indx,0];             print('  \tTsim\t:\t', Tsim)
qx     = sum(data[indx:m-1,5]*data[indx:m-1,10])/sum(data[indx:m-1,10]); 
Ub     = qx/2;                                               print('  \tUb\t:\t', Ub)
dudy   = 0.5*sum((data[indx:m-1,1]+data[indx:m-1,2])*
                  data[indx:m-1,10])/sum(data[indx:m-1,10]); print('  \tdudy\t:\t', dudy)
ni     = 1/Re;					       print('  \tni\t:\t', ni)
tw     = dudy*ni;
P0     = tw*Ub;                                              print('  \tP0\t:\t', P0)
utau   = sqrt(tw);				               print('  \tutau\t:\t', utau)
Ret    = utau/ni;					       print('  \tRetau\t:\t', Ret)
Reb    = Ub/ni;					       print('  \tReb\t:\t', Reb)
Cf     = 2*tw/power(Ub,2);				       print('  \tCf\t:\t', Cf)
Cf_dean= 0.073*power(Reb,-0.25);			       print('  \tCf\t:\t', \
                                                                    Cf_dean,'  (Dean)')
Cf_pope= 0.0336*power(Ret,-0.273);                           print('  \tCf\t:\t',
                                                                    Cf_pope,'  (Pope)')
meandt = sum(power(data[indx:m-1,10],2))/sum(data[indx:m-1,10])
if iscontrolled:
  if (control=='TW') or (control=='OW'):
    wdwdy = 0.5*sum((pdata[indx:h-2,1]-pdata[indx:h-2,4])*
                    (pdata[indx+1:h-1,0]-pdata[indx:h-2,0]))/ \
                sum(pdata[indx+1:h-1,0]-pdata[indx:h-2,0]);print('  \twdwdy\t:\t', wdwdy)
    wdwdyN= 0.5*sum((pdata[indx:h-2,3]-pdata[indx:h-2,5])*
                    (pdata[indx+1:h-1,0]-pdata[indx:h-2,0]))/ \
                sum(pdata[indx+1:h-1,0]-pdata[indx:h-2,0]);print('  \twdwdyN\t:\t', wdwdyN)
    wdwdyP= 0.5*sum((pdata[indx:h-2,2]-pdata[indx:h-2,6])*
                    (pdata[indx+1:h-1,0]-pdata[indx:h-2,0]))/ \
                sum(pdata[indx+1:h-1,0]-pdata[indx:h-2,0]);print('  \twdwdyN\t:\t', wdwdyP)
    Pin = ni*wdwdy;                                        print('  \tPin\t:\t', Pin)
    PinN = ni*wdwdyN;                                      print('  \tPinN\t:\t', PinN)
    PinP = ni*wdwdyP;                                      print('  \tPinP\t:\t', PinP)
    
  elif control=='VC':
    Pin = 0.5*sum((pdata[indx:h-2,1]-pdata[indx:h-2,2])*
                    (pdata[indx+1:h-1,0]-pdata[indx:h-2,0]))/ \
              sum(pdata[indx+1:h-1,0]-pdata[indx:h-2,0]);  print('  \twdwdy\t:\t', Pin)
    PinN = Pin; PinP = Pin
  gamma = Pin/P0;                                          print('  \tgamma\t:\t', gamma)

# Computing uncertainties
sbar_dudy=0.0; sbar_Ub=0.0; T0_dudy=0.0; T0_Ub=0.0
if elaborate:
  t = linspace(data[indx,0],data[m-2,0],m-indx-1); dt=t[1]-t[0]
  dudyI0 = zeros((1,m-indx-1)); dudyIn = zeros((1,m-indx-1));
  dudyI0f = interp1d(data[indx-1:m-1,0], data[indx-1:m-1,1]); dudyI0[0,:] = dudyI0f(t)
  dudyInf = interp1d(data[indx-1:m-1,0], data[indx-1:m-1,2]); dudyIn[0,:] = dudyInf(t)
  ar_dudyI0 = arsel(dudyI0); ar_dudyIn = arsel(dudyIn)
  zi0  = lfiltic([1], ar_dudyI0.AR[0], ar_dudyI0.autocor[0])
  zin  = lfiltic([1], ar_dudyIn.AR[0], ar_dudyIn.autocor[0])
  rho0 = ones(m-indx-1); rho0[1:] = lfilter([1], ar_dudyI0.AR[0], \
              zeros(m-indx-2), zi=zi0)[0]
  rhon = ones(m-indx-1); rhon[1:] = lfilter([1], ar_dudyIn.AR[0], \
              zeros(m-indx-2), zi=zin)[0]
  sbar_dudy = std(concatenate((data[indx:m-1,1], data[indx:m-1,2])))/ \
              sqrt(ar_dudyI0.eff_N[0]+ar_dudyIn.eff_N[0])
  T0_dudy =  0.5*dt*(ar_dudyI0.T0[0]+ar_dudyIn.T0[0])
  if CPI:
    UbI = zeros((1,m-indx-1))
    UbIf = interp1d(data[indx-1:m-1,0], data[indx-1:m-1,5]); UbI[0,:] = 0.5*UbIf(t)
    ar_UbI  = arsel(UbI); ziUb  = lfiltic([1], ar_UbI.AR[0], ar_UbI.autocor[0])
    rhoU = ones(m-indx-1); rhoU[1:] = lfilter([1], ar_UbI.AR[0],  \
                zeros(m-indx-2), zi=ziUb)[0]
    sbar_Ub = std(data[indx:m-1,5])/sqrt(ar_UbI.eff_N[0])
    T0_Ub = dt*(ar_UbI.T0[0])

  print('  sigma\tdudy \t:\t', sbar_dudy)
  print('  Teff\tdudy \t:\t', T0_dudy)
  print('  sigma\tUb \t:\t', sbar_Ub)
  print('  Teff\tUb \t:\t', T0_Ub)
  
print('  ================================== \n\n')

# Plotting wall shear
plt.figure(1)
plt.plot(data[:,0], data[:,1], 'm-')
plt.plot(data[:,0], data[:,2], 'c-')
plt.plot([data[indx,0], data[indx,0]], transpose([plt.ylim()]), 'k-')
plt.plot(transpose([plt.xlim()]), [dudy, dudy], 'k-')
plt.xlabel('$t U_{ref}/H$')
plt.ylabel('$d\overline{u}/dy$')

# Plotting bulk velocity
plt.figure(2)
plt.plot(data[:,0], 0.5*data[:,5], 'm-')
plt.plot([data[indx,0], data[indx,0]], transpose([plt.ylim()]), 'k-')
plt.plot(transpose([plt.xlim()]), [Ub, Ub], 'k-')
plt.xlabel('$t U_{ref}/H$')
plt.ylabel('$\overline{U_b}$')
plt.show()


