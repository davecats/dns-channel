%
% This is a sample script
% that does stuff
%

addpath ./base
addpath ./EMD

% Read field file
[U,Ubar,Wbar,dns]=veta2uvw('Field200.fld');
% Set some parameters for EMD
parameter.type='hD';
parameter.n=1;
parameter.MNAI=1;
% Start empirical mode decomposition
disp('Performing empirical mode decomposition')
EMDC=zeros(2*dns.nz+1,2*dns.nx+1,parameter.n+1,dns.ny+1);
parfor iy=1:dns.ny+1
  EMDC(:,:,:,iy)=FABEMD(squeeze(U(1,:,:,iy)),parameter,1,1);
end
for iy=1:dns.ny+1
    EMDC(:,:,end,iy)=EMDC(:,:,end,iy)-Ubar(iy);
end
% Plot some large scale isosurfaces
levels = [-1.2];
figure(1); hold on; box on
for lev=levels
    isosurface(squeeze(EMDC(:,:,2,:)),lev);
end