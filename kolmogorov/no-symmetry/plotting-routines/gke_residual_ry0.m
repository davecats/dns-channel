clear variables
close all

ny=128; a=2.1;	         
nx=14; nxc=48;
nz=14; nzc=48;
beta0=8.9760;
alfa0=4;
path = '/home/davide/Desktop/no-symm/';

f=fopen(strcat(path,'residual.bin')); res=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'energy.bin'));   ener=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'phiry.bin'));    phiry=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'dphiry.bin'));   dphiry=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'phiC.bin'));    phiC =reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'dphiC.bin'));    dphiC =reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'dphirx.bin'));   dphirx=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen(strcat(path,'dphirz.bin'));   dphirz=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);


figure(); contourf(0.5*squeeze(res(:,10,1:end)+res(:,10,end:-1:1))'); title('res'); ylim([1 66])
figure(); contourf(0.5*squeeze(ener(:,10,1:end)+ener(:,10,end:-1:1))'); title('energy'); ylim([1 66])
figure(); contourf(0.5*squeeze(phiry(:,10,1:end)-phiry(:,10,end:-1:1))'); title('phiry'); ylim([1 66])
figure(); contourf(0.5*squeeze(dphiry(:,10,1:end)+dphiry(:,10,end:-1:1))'); title('dphiry'); ylim([1 66])
figure(); contourf(0.5*squeeze(dphirx(:,10,1:end)+dphirx(:,10,end:-1:1))'); title('dphirx'); ylim([1 66])
figure(); contourf(0.5*squeeze(dphirz(:,10,1:end)+dphirz(:,10,end:-1:1))'); title('dphirz'); ylim([1 66])
figure(); contourf(0.5*squeeze(phiC(:,10,1:end)-phiC(:,10,end:-1:1))'); title('phiC'); ylim([1 66])
figure(); contourf(0.5*squeeze(dphiC(:,10,1:end)+dphiC(:,10,end:-1:1))'); title('dphiC'); ylim([1 66])