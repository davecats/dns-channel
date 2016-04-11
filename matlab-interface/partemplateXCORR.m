%
% This is a sample script
% that does stuff
%

% INPUT
% ------------------------------
EMD.path = '/net/istmuranus/volume1/shared/marion/Re1000/nx-1024-nz-1024/00/';
Vfield.path = '/net/istmuranus/volume1/data/hi221/BACKUP-WORK-HLR1/CONDITIONAL-AVERAGE/Re1000/nx-1024-nz-1024/00/postpro-fld/';
PSD.fname='./PSD.bin';
Vfield.nmin=1;
Vfield.nmax=23;
fftw('planner','patient')
nproc=8;
% ------------------------------

clc
addpath ./base
addpath ./emd
addpath ./psd

% Map field file to disk
[dns,Ubar,Wbar,VDisk]=open_fld(strcat(Vfield.path,'Field',num2str(Vfield.nmin),'.fld'));
% Create vector on disk for PSD and XCORR
PSDimage=storedstructure(PSD.fname,{'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'PSD';
                            'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'CORR'}, ...
                            true);
for iy=1:dns.ny+1
 PSDimage.data.PSD(:,:,iy)=zeros([2*dns.nz+1 2*dns.nx+1]);
 PSDimage.data.CORR(:,:,iy)=zeros([2*dns.nz+1 2*dns.nx+1]);
end
% Total number of fields
Vfield.ntot=Vfield.nmax-Vfield.nmin+1;
  
% Define y-coordinates
y=zeros(dns.ny+3,1);
for i=1:dns.ny+3
    y(i)=dns.ymin+0.5*(dns.ymax-dns.ymin)* ...
         (tanh(dns.a*(2*(i-2)/dns.ny-1))/tanh(dns.a)...
        +0.5*(dns.ymax-dns.ymin));
end

spmd(nproc)
for n=Vfield.nmin:Vfield.nmax   
disp(['Reading Field ' num2str(n)]);
% Map field file to disk
EMD.fname=strcat(EMD.path,'EMD',num2str(n),'.bin');
EMDimage=storedstructure(EMD.fname, {'double' [2*dns.nz+1 2*dns.nx+1 2+1 dns.ny+1] 'EMDC';
                    'double' 2 'EMD'}, true);
dUdy = FD(squeeze(EMDimage.data.EMDC(:,:,3,1:5)),y(2:6),4); s1=var(dUdy(:)); V1 = fft2(dUdy)/(2*dns.nx+1)/(2*dns.nz+1); V1(1,1)=0;
ny=round(dns.ny/2);
for i=0+floor((labindex-1)*(ny+1)/numlabs):floor(labindex*(ny+1)/numlabs)-1
    iy=i+2; IY=i+1;
    % Compute u-XCORR
    V2=squeeze(EMDimage.data.EMDC(:,:,3,IY)); s2=var(V2(:)); V2=fft2(V2)/(2*dns.nx+1)/(2*dns.nz+1); V2(1,1)=0; V = conj(V1).*V2;
    PSDimage.data.CORR(:,:,IY)= PSDimage.data.CORR(:,:,IY) +...
       ifftshift(real(ifft2(V)))*(2*dns.nx+1)*(2*dns.nz+1)/Vfield.ntot/sqrt(s1*s2);
    PSDimage.data.PSD(:,:,IY)=PSDimage.data.PSD(:,:,IY)+real(V)/Vfield.ntot;
end
labBarrier
end
end
% Plot correlation
disp('Plotting PSD')
PSDimage=PSDimage{1};
figure(3); hold on; box on
contourf(squeeze(PSDimage.data.CORR(dns.nz+5,:,:))')