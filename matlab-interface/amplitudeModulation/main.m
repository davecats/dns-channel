close all
clear all

re=8000;
iy=20;
path='../';


% Load data & coords
load(strcat(path,'emd_scales_u1_',num2str(re),'_',num2str(iy)));
load(strcat(path,'coords_',num2str(re),'_','116-338.mat'));
t=t-t(1);

% Space filtering
dz=mean(diff(zz/ll)); freq=1/dz; cutoff=1/400;

figure(1); hold on
plot(zz,u_Splus(1,:),'k')
plot(zz,envl(u_Splus(1,:),1,cutoff,freq),'linewidth',2)
plot(zz,envl(u_Splus(1,:),-1,cutoff,freq),'linewidth',2)
plot(zz,envl_mathis(u_Splus(1,:),cutoff,freq),'linewidth',2)
legend('SS signal','E_{SS,up}, Agostini et al. (2016)','E_{SS,low}, Agostini et al. (2016)','E_{SS} Mathis et al. (2009)')
xlabel('z^+')

% Time filtering
dt=mean(diff(t)); freq=1/dt; cutoff=1/400;

figure(2); hold on
plot(t,u_Splus(:,1),'k')
plot(t,envl(u_Splus(:,1),1,cutoff,freq),'linewidth',2)
plot(t,envl(u_Splus(:,1),-1,cutoff,freq),'linewidth',2)
plot(t,envl_mathis(u_Splus(:,1),cutoff,freq),'linewidth',2)
legend('SS signal','E_{SS,up}, Agostini et al. (2016)','E_{SS,low}, Agostini et al. (2016)','E_{SS} Mathis et al. (2009)')
xlabel('t^+')
