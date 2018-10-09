clear all
close all
%clc

% Inputs
% ---------
iF = 1; %field number
params.nxd=192;
params.nzd=192;
params.ny=193;
params.a=1.6;
params.alfa0=0.8;
params.beta0=2;
params.y = tanh( params.a * ( 2 * ( (-1:params.ny+1) )/params.ny - 1) )/tanh( params.a ) + 1;
params.x = (0:params.nxd-1)*(2*pi/params.alfa0)/params.nxd;
params.z = (0:params.nzd-1)*(2*pi/params.beta0)/params.nzd;
params.size = [3,params.ny+3,params.nzd,params.nxd];
% ---------

addpath('./base');
%% load field (v,eta)
[field,params]=readfield(strcat('../rField',num2str(iF),'.bin'),params);
%% setup derivatives (run only once)
derivatives=compute_derivatives(params);
%% remove mean velocity component (if needed)
field.U=field.U-repmat(mean(mean(field.U(:,:,:,:),3),4),[1,1,params.size(3:4)]);
%% compute velocity gradient tensor
field.dU=zeros(3,params.ny+3,params.nzd,params.nxd,3);    % first index: velocity component, last index: derivative direction
% in x-direction derivative in Fourier
field.dU(:,:,:,:,1) = real(ifft(repmat(1j*reshape(params.alfa0*[0:params.nxd/2, -params.nxd/2+1:-1],[1,1,1,params.nxd]),[params.size(1:3),1]).*fft(field.U,[],4),[],4));
% in z-direction derivative in Fourier
field.dU(:,:,:,:,3) = real(ifft(repmat(1j*reshape(params.beta0*[0:params.nzd/2, -params.nzd/2+1:-1],[1,1,params.nzd,1]),[params.size(1:2),1,params.size(4)]).*fft(field.U,[],3),[],3));
% in y-direction compact finite differences
for iV=1:3
   field.dU(iV,:,:,:,2) = reshape(derivatives.d0\(derivatives.d1*reshape(field.U(iV,:,:,:),[params.size(2),params.size(3)*params.size(4)])),params.size(2:4));
end
%% vorticity vector
field.omega=zeros(3,params.ny+3,params.nzd,params.nxd);
field.omega(1,:,:,:)=field.dU(3,:,:,:,2)-field.dU(2,:,:,:,3);
field.omega(2,:,:,:)=field.dU(1,:,:,:,3)-field.dU(3,:,:,:,1);
field.omega(3,:,:,:)=field.dU(2,:,:,:,1)-field.dU(1,:,:,:,2);
%% divergence of vorticity
field.divomega=real(ifft(repmat(reshape(1j*(params.alfa0*[0:params.nxd/2, -params.nxd/2+1:-1]),[1,1,1,params.nxd]),[1,params.size(2:3),1]).*fft(field.omega(1,:,:,:),[],4),[],4))+ ... d( lapl(omegax) )/dx
               real(ifft(repmat(reshape(1j*(params.beta0*[0:params.nzd/2, -params.nzd/2+1:-1]),[1,1,params.nzd,1]),[1,params.size(2:2),1,params.size(4)]).*fft(field.omega(3,:,:,:,:),[],3),[],3))+... % d( lapl(omegaz) )/dz
               reshape(derivatives.d0\(derivatives.d1*reshape(field.omega(2,:,:,:,:),[params.size(2),params.size(3)*params.size(4)])),[1,params.size(2:4)]); % d( lapl(omegay))/dy
%% laplacian of the vorticity
field.d2omega=zeros(3,params.ny+3,params.nzd,params.nxd,3);    % first index: vorticity component, last index: derivative direction
% in x-direction derivative in Fourier
field.d2omega(:,:,:,:,1) = real(ifft(repmat(-reshape((params.alfa0*[0:params.nxd/2, -params.nxd/2+1:-1]).^2,[1,1,1,params.nxd]),[params.size(1:3),1]).*fft(field.omega,[],4),[],4));
% in z-direction derivative in Fourier
field.d2omega(:,:,:,:,3) = real(ifft(repmat(-reshape((params.beta0*[0:params.nzd/2, -params.nzd/2+1:-1]).^2,[1,1,params.nzd,1]),[params.size(1:2),1,params.size(4)]).*fft(field.omega,[],3),[],3));
% in y-direction compact finite differences
for iV=1:3
   field.d2omega(iV,:,:,:,2) = reshape(derivatives.d0\(derivatives.d2*reshape(field.omega(iV,:,:,:),[params.size(2),params.size(3)*params.size(4)])),params.size(2:4));
end
%% divergence of laplacian of vorticity
field.divd2omega=real(ifft(repmat(1j*reshape((params.alfa0*[0:params.nxd/2, -params.nxd/2+1:-1]),[1,1,1,params.nxd]),[1,params.size(2:3),1]).*fft(sum(field.d2omega(1,:,:,:,:),5),[],4),[],4))+ ... d( lapl(omegax) )/dx
                 real(ifft(repmat(1j*reshape((params.beta0*[0:params.nzd/2, -params.nzd/2+1:-1]),[1,1,params.nzd,1]),[1,params.size(2:2),1,params.size(4)]).*fft(sum(field.d2omega(3,:,:,:,:),5),[],3),[],3))+... % d( lapl(omegaz) )/dz
                 reshape(derivatives.d0\(derivatives.d1*reshape(sum(field.d2omega(2,:,:,:,:),5),[params.size(2),params.size(3)*params.size(4)])),[1,params.size(2:4)]); % d( lapl(omegay) )/dy
%% vorticity magnitude
om = squeeze(sqrt(sum(field.omega(:,:,:,:).^2,1)));
[X,Y,Z]=meshgrid(params.z,params.y,params.x);

%% Plot something
zslice=[0 1 2]; %linspace(0,params.x(end),10)
yslice=[];
xslice=[0.2, 1];
slice(X,Y,Z,squeeze(field.divd2omega),zslice,yslice,xslice)
shading interp
axis equal
xlabel('z')
ylabel('y')
zlabel('x')

