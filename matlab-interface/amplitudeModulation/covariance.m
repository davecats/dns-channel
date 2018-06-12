close all
clear all

% Input parameters
re=8000;
path='../';
iy=1:47;
dt=0.8; freq=1/dt; cutoff=1/400;
% ----------------

%% Declare variables
ny=numel(iy);
CAM=zeros(ny,ny);
load(strcat(path,'coords_',num2str(re),'_','116-338.mat'));
fname = @(re,iy) strcat(path,'emd_scales_u1_',num2str(re),'_',sprintf('%02d',iy));

%% Compute envelopes in time and load large scales
disp('Loading data and computing envelopes...')
for iy1=iy
    disp(['        ' num2str(iy1) ' of ' num2str(ny)])
    u1=load(fname(re,iy1));
    ESS{iy1}=envl_mathis(u1.u_Splus,cutoff,freq);
    LS{iy1}=u1.u_Lplus;
end

%% Start computing the covariance
disp('Computing covariance...')
for iy1=iy
    disp(['        ' num2str(iy1) ' of ' num2str(ny)])
    for iy2=iy
      c=cov(ESS{iy1},LS{iy2});
      CAM(iy1,iy2)=c(1,2);
    end
end

%% Save it
save('../cov_ESS_uLS.mat','CAM')

%% Plot it
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot, 'defaultTextInterpreter','latex');
figure
box on
pcolor(yy/ll,yy/ll,CAM); shading interp
cb=colorbar; cb.TickLabelInterpreter = 'latex';
set(gca(),'XScale','log')
set(gca(),'YScale','log')
xlabel('$y_1^+$')
ylabel('$y_2^+$')
