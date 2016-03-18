%
% This is a sample script
% that does stuff
%

clc
addpath ./base
addpath ./emd
addpath ./psd

% Read field file

% Convert (v,eta) to (u,v,w)
[U,Ubar,Wbar,dns]=veta2uvw('/home/davide/MATTEST/matlab-interface/Field200.fld');
% Compute the PSD
PSD=plane_psd(squeeze(U(1,:,:,:)),squeeze(U(1,:,:,:)));
% Convert to physical space 
U=IFFTU(U);
% Remove ghost planes
[U,Ubar,Wbar]=remove_ghost(U,Ubar,Wbar);
% Set some parameters for EMD
parameter.type='hD';
parameter.n=1;
parameter.MNAI=1;
% Start empirical mode decomposition
disp('Performing empirical mode decomposition')
EMDC=zeros(2*dns.nz+1,2*dns.nx+1,parameter.n+1,dns.ny+1);
h=waitbar(0,'EMD: please wait...');
parfor iy=1:dns.ny+1
  waitbar(iy/(dns.ny+1),h);
  EMDC(:,:,:,iy)=FABEMD(squeeze(U(1,:,:,iy)),parameter,0);
end
for iy=1:dns.ny+1
    EMDC(:,:,end,iy)=EMDC(:,:,end,iy)-Ubar(iy);
end
% Plot some large scale isosurfaces
disp('Plotting some large-scale isosurfaces')
levels = [-1.2,0.8];
figure(1); hold on; box on
for lev=levels
    isosurface(squeeze(EMDC(:,:,end,:)),lev);
end
