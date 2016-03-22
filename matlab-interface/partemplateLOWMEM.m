%
% This is a sample script
% that does stuff
%

% INPUT
% ------------------------------
file = '/home/davide/MATTEST/matlab-interface/Field200.fld';
Vdisk_fname  = './Field.bin';
EMD.fname = './EMD.bin';
EMD.type='hD';   
EMD.n=1;         
EMD.MNAI=1;
PSD.fname = './PSD.bin';
PSD.y = 0.075;
fftw('planner','patient')
nproc=2;
% ------------------------------

clc
addpath ./base
addpath ./emd
addpath ./psd

% Map field file to disk
[dns,Ubar,Wbar,VDisk]=open_fld(file);
% Create vector on disk for velocity output
fh = fopen(Vdisk_fname,'w'); 
for iy=0:dns.ny
 fwrite(fh,zeros(3,2*dns.nz+1,2*dns.nx+1),'double'); 
end
fclose(fh);
Uimage=memmapfile(Vdisk_fname, ...
                 'Format', {'double' [3 2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'U'}, ...
                 'Writable',true);
% Create vector on disk for EMD
fh = fopen(EMD.fname,'w'); 
for iy=0:dns.ny
    fwrite(fh,zeros(2*dns.nz+1,2*dns.nx+1,EMD.n+1),'double'); 
end
fclose(fh);
EMDimage=memmapfile(EMD.fname, ...
                 'Format', {'double' [2*dns.nz+1 2*dns.nx+1 EMD.n+1 dns.ny+1] 'EMDC'}, ...
                 'Writable',true);
% Create vector on disk for PSD and XCORR
fh = fopen(PSD.fname,'w'); 
for iy=0:dns.ny
    fwrite(fh,zeros(2,2*dns.nz+1,2*dns.nx+1),'double'); 
end
fclose(fh);
PSDimage=memmapfile(PSD.fname, ...
                 'Format', {'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'PSD';
                            'double' [2*dns.nz+1 2*dns.nx+1 dns.ny+1] 'CORR'}, ...
                 'Writable',true);
             
% Define y-coordinates
y=zeros(dns.ny+3,1);
for i=1:dns.ny+3
    y(i)=dns.ymin+0.5*(dns.ymax-dns.ymin)* ...
         (tanh(dns.a*(2*(i-2)/dns.ny-1))/tanh(dns.a)...
        +0.5*(dns.ymax-dns.ymin));
end

spmd(nproc)
for i=0+floor((labindex-1)*(dns.ny+1)/numlabs):floor(labindex*(dns.ny+1)/numlabs)-1
    tic
    % Define helping indices
    iy=i+2; IY=i+1;
    % Extract a slice consisting of three planes
    U=complex(zeros(3,2*dns.nz+1,2*dns.nx+1,3));
    % Copy a slice to RAM memory
    U(2:3,:,dns.nx+1:2*dns.nx+1,:)=complex(VDisk.data.V(1,:,:,:,iy-1:iy+1),...
                                           VDisk.data.V(2,:,:,:,iy-1:iy+1));
    % Convert slice from (v,eta) to (u,v,w)
    U=veta2uvwIY(U,Ubar(iy-1:iy+1),Wbar(iy-1:iy+1),y(iy-1:iy+1),dns);
    % Convert to physical space 
    for iV=1:3
      U(iV,:,:,2)=plane_ifft(U(iV,:,:,2),dns.nx,dns.nz);
    end
    % Keep only central plane and convert to double
    U=real(U(:,:,:,2));
    % Empirical Mode Decomposition
    EMDC=zeros(2*dns.nz+1,2*dns.nx+1,EMD.n+1);
    EMDC(:,:,:)=FABEMD(squeeze(U(1,:,:,1)),EMD,0);
    EMDC(:,:,end)=EMDC(:,:,end)-Ubar(iy); 
    % Write the slice to disk
    Uimage.data.U(:,:,:,IY)=U;
    EMDimage.data.EMDC(:,:,:,IY)=EMDC;
    toc
end
% Compute PSD and XCORR of large scale
iy = find(y>=PSD.y,1,'first')-1; % reference plane
EMDCiy=fft2(EMDimage.data.EMDC(:,:,1,iy))/(2*dns.nx+1)/(2*dns.nz+1);
for i=0+floor((labindex-1)*(round(dns.ny/2)+1)/numlabs):floor(labindex*(round(dns.ny/2)+1)/numlabs)-1
    IY=i+1;
    EMDCIY=fft2(EMDimage.data.EMDC(:,:,1,IY))/(2*dns.nx+1)/(2*dns.nz+1);
    PSDimage.data.CORR(:,:,IY)=...
       ifftshift(real(ifft2(conj(EMDCiy).*EMDCIY)))*(2*dns.nx+1)*(2*dns.nz+1);
    PSDimage.data.PSD(:,:,IY)=...
       real(conj(EMDCiy).*EMDCIY);
end
end
% Plot some large scale isosurfaces from disk!
disp('Plotting some large-scale isosurfaces')
EMDimage=EMDimage{1}; PSDimage=PSDimage{1};
levels = [-1.2,0.8];
figure(1); hold on; box on
for lev=levels
    isosurface(squeeze(EMDimage.data.EMDC(:,:,end,:)),lev);
end