close all
clear all

%% Input parameters
re=8000;    % reynolds number
path='../'; % path
iy=1:47;    % array of iy indices
filt=3;     % filter type, see below
nproc=24;   % number of parallel processes
nSS=12;     % number of IMF to be considered small scales
% -----------------
switch filt
    case 1 % Time cutoff-filter
        dt=0.8; freq=1/dt; cutoff=1/400;
        filter.fun=@lowpass; filter.params={cutoff freq 'ImpulseResponse' 'iir' 'Steepness' 0.95 'StopbandAttenuation' 100};
    case 2 % Spanwise cutoff filter
        dz=(zz(2)-zz(1))/ll; freq=1/dz; cutoff=1/400;
        filter.fun=@lowpass; filter.params={cutoff freq 'ImpulseResponse' 'iir' 'Steepness' 0.95 'StopbandAttenuation' 100};
    case 3 % EMD filter
        filter.fun=@FABEMD_philipp; EMDparams.type='LD'; EMDparams.n=15; EMDparams.MNAI=1; filter.params={EMDparams};
end
% ----------------

%% Declare variables
ny=numel(iy);
CAM=zeros(ny,ny);
load(strcat(path,'coords_',num2str(re),'_','116-338.mat'));
fname = @(re,iy) strcat(path,'emd_scales_u1_',num2str(re),'_',sprintf('%02d',iy));

%% Compute envelopes and large scales
disp('Loading data and computing envelopes...')
m = storedstructure('AM.bin',...
                    {'double',[numel(t) numel(zz) EMDparams.n.*(filt==3)+1 ny],'E_SS';
                     'double',[numel(t) numel(zz) ny],            'u_LS'},...
                    true);
spmd(nproc)
    for iy1=iy(1)+floor((labindex-1)*iy(end)/numlabs):floor(labindex*iy(end)/numlabs)-1+iy(1)
        disp(['        ' num2str(iy1) ' of ' num2str(ny)])
        u1=load(fname(re,iy1));
        m.Data.E_SS(:,:,:,iy1)=envl_mathis(u1.u_Splus,filter.fun,filter.params);
        m.Data.u_LS(:,:,iy1)=u1.u_Lplus;
    end
end
m=m{1}; % from composite to standard workspace variable
m.Writable=false; % make storedstructure read-only

%% Load envelopes and large scales to memory
disp('Loading envelopes and large scales to memory...')
E_SS=zeros(numel(t),numel(zz),ny);
u_LS=zeros(numel(t),numel(zz),ny);
for iy1=iy
    E_SS(:,:,iy1)=squeeze(sum(m.Data.E_SS(:,:,nSS*(filt==3)+1:end,iy1),3));
    u_LS(:,:,iy1)=m.Data.u_LS(:,:,iy1);
end

%% Start computing the covariance
disp('Computing covariance...')
for iy1=iy
    disp(['        ' num2str(iy1) ' of ' num2str(ny)])
    for iy2=iy
      c=cov(squeeze(E_SS(:,:,iy1)),squeeze(u_LS(:,:,iy2)));
      CAM(iy1,iy2)=c(1,2);
    end
end

%% Save it
save(strcat('./cov_ESS_uLS_nSS_',num2str(nSS),'.mat'),'CAM')

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
set(gca(),'Layer','top')