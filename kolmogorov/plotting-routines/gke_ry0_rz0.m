%% Initialization
%  ==============

clear variables
close all

ny=128; a=2.1;	         
nx=14; nxc=48;
nz=14; nzc=48;
beta0=8.9760;
alfa0=4;
file = '/run/media/davide/812f81ed-fed2-4bb6-af2d-2b64759fb87f/kolmo-mfu/gke.bin';

% Memory map to gke data (this may fail if data is too big)
gkedata = memmapfile(file, 'Format', {'double', [6,nzc,nxc,ny+3,ny+3], 'gke'});

% Define y-coordinates
y=zeros(ny+3,1);
for i=1:ny+3
    y(i)=tanh(a*(2*(i-2)/ny-1))/tanh(a)+1;
end
y_2 = y(1:floor(ny/2)+2);

% Load hyperplane rx=0
gke=squeeze(gkedata.Data.gke(:,1,1,:,:));  % i, rz, ry, y2, y1

% Compute (Yc,ry)
Yc = zeros(ny+3); ry = zeros(ny+3);
for y1=-1:ny+1
    Y1=y1+2;
    for y2=-1:ny+1
        Y2=y2+2;
        Yc(Y2,Y1)=0.5*(y(Y2)+y(Y1));
        ry(Y2,Y1)=(y(Y2)-y(Y1));
    end
end

% Define cartesian grid for (Yc,ry)
[RY,YC]=meshgrid(linspace(0,2,4*ny),linspace(0,1,2*ny));

% Interpolate data onto the cartesian grid
F = scatteredInterpolant(ry(:),Yc(:),squeeze(gke(2,:))','natural','none');
phiry=F(RY,YC);
F = scatteredInterpolant(ry(:),Yc(:),squeeze(gke(4,:))','natural','none');
phiC=F(RY,YC);
F = scatteredInterpolant(ry(:),Yc(:),squeeze(gke(5,:))','natural','none');
source=F(RY,YC);
% 
% % Remove duplicate data
% for yc=1:ny/2+2
%   for ry=1:ny+3
%     if y(ry)>2*y(yc)
%       phiry(yc,ry)=NaN;
%       phiC(yc,ry)=NaN;
%       source(yc,ry)=NaN;
%     end
%   end
% end

% Plot
figure(); box on; hold on
[~,c] = contourf(RY,YC,source,10); set(c,'LineStyle','None');
q = quiver(RY(1:6:end,1:6:end),YC(1:6:end,1:6:end),phiry(1:6:end,1:6:end),phiC(1:6:end,1:6:end));
set(q, 'LineWidth',1,'Color','k','AutoScaleFactor',1.5);
s = streamline(RY,YC,phiry,phiC,[0.2 0.3 0.4 0.6 0.8 1 ],[0.2 0.3 0.3 0.4 0.5 0.6]);
set(s,'Color','w','Linewidth',2)
xlabel('r_y')
ylabel('Y_c')