ny=128; a=2.1;	         
nx=14; nxc=48;
nz=14; nzc=48;
beta0=8.9760;
alfa0=4;

f=fopen('residual.bin');  res=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen('dphiry.bin'); dphiry=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen('dphiC.bin');  dphiC =reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen('dphirx.bin'); dphirx=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);
f=fopen('dphirz.bin'); dphirz=reshape(fread(f,'double'),nzc,nxc,ny+3); fclose(f);


figure(); contourf(0.5*squeeze(res(:,1,1:end)+res(:,1,end:-1:1))'); title('res')
figure(); contourf(squeeze(dphiry(:,1,1:end))'); title('dphiry')
figure(); contourf(squeeze(dphirx(:,1,1:end))'); title('dphirx')
figure(); contourf(squeeze(dphirz(:,1,1:end))'); title('dphirz')
figure(); contourf(squeeze(dphiC(:,1,1:end))'); title('dphiC')