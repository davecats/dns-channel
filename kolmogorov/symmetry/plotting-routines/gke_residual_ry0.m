
clear variables
close all

ny=128; a=2.1;	         
nx=14; nxc=48;
nz=14; nzc=48;
beta0=8.9760;
alfa0=4;
path = '/home/davide/Desktop/symm/';

f=fopen(strcat(path,'residual.bin'));  res=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'energy.bin'));  ener=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'phiry.bin')); phiry=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'dphiry.bin')); dphiry=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'phiC.bin'));  phiC =reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'dphiC.bin'));  dphiC =reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'dphirx.bin')); dphirx=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);
f=fopen(strcat(path,'dphirz.bin')); dphirz=reshape(fread(f,'double'),nzc,nxc,ny/2+2); fclose(f);


figure(); contourf(squeeze(res(:,10,1:end))'); title('res')
figure(); contourf(squeeze(ener(:,10,1:end))'); title('energy')
figure(); contourf(squeeze(phiry(:,10,1:end))'); title('phiry')
figure(); contourf(squeeze(dphiry(:,10,1:end))'); title('dphiry')
figure(); contourf(squeeze(dphirx(:,10,1:end))'); title('dphirx')
figure(); contourf(squeeze(dphirz(:,10,1:end))'); title('dphirz')
figure(); contourf(squeeze(phiC(:,10,1:end))'); title('phiC')
figure(); contourf(squeeze(dphiC(:,10,1:end))'); title('dphiC')
