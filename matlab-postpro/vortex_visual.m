% Load Field
% --------------------
load Field10.mat

% Compute gradients
% --------------------
deltaX = (2*pi/dns.alfa0)/(2*dns.nx+1);
deltaZ = (2*pi/dns.beta0)/(2*dns.nz+1);
Lx = 2*pi/dns.alfa0; Ly = 2*pi/dns.beta0;
uy = zentral_diff(squeeze(data2(1,:,:,:)),y,dns.nx,dns.ny,dns.nz); uy=uy(2:end-1,2:end-1,:);
vy = zentral_diff(squeeze(data2(2,:,:,:)),y,dns.nx,dns.ny,dns.nz); vy=vy(2:end-1,2:end-1,:);
wy = zentral_diff(squeeze(data2(3,:,:,:)),y,dns.nx,dns.ny,dns.nz); wy=wy(2:end-1,2:end-1,:);
data2=data2(:,:,:,2:dns.ny+2);
uz = squeeze(data2(1,3:end,2:end-1,:)-data2(1,1:end-2,2:end-1,:))/(2*deltaZ);
vz = squeeze(data2(2,3:end,2:end-1,:)-data2(2,1:end-2,2:end-1,:))/(2*deltaZ);
wz = squeeze(data2(3,3:end,2:end-1,:)-data2(3,1:end-2,2:end-1,:))/(2*deltaZ);
ux = squeeze(data2(1,2:end-1,3:end,:)-data2(1,2:end-1,1:end-2,:))/(2*deltaX);
vx = squeeze(data2(2,2:end-1,3:end,:)-data2(2,2:end-1,1:end-2,:))/(2*deltaX);
wx = squeeze(data2(3,2:end-1,3:end,:)-data2(3,2:end-1,1:end-2,:))/(2*deltaX);

% Define grad(V)
% ---------------------
gradU = @(ix,iy,iz) [ux(ix,iz,iy),uy(ix,iz,iy),uz(ix,iz,iy); ...
                     vx(ix,iz,iy),vy(ix,iz,iy),vz(ix,iz,iy); ...
                     wx(ix,iz,iy),wy(ix,iz,iy),wz(ix,iz,iy)];
                 
Omega = @(ix,iy,iz) 0.5*(gradU(ix,iy,iz)-gradU(ix,iy,iz)');

Eps =   @(ix,iy,iz) 0.5*(gradU(ix,iy,iz)+gradU(ix,iy,iz)');

% Compute criteria
% ---------------------
[m,n,p] = size(ux);
Q=ux(:,:,1:round(p/2))*0; Delta=Q; lambda2=Q; lambdaci=Q;
for ix=1:m
    for iy=1:round(p/2)
       for iz=1:n
         O = Omega(ix,iy,iz); E=Eps(ix,iy,iz);
         Q(ix,iz,iy) = 0.5*(trace(O*O')-trace(E*E'));
         Delta(ix,iz,iy) = (0.5*det(gradU(ix,iy,iz)))^2+((1/3)*Q(ix,iz,iy))^3;
         lambda  = sort(eig(E*E+O*O),'descend');
         lambda2(ix,iz,iy) = lambda(2);
         lambda_gradU = eig(gradU(ix,iy,iz));
         for i=1:3
             if isreal(lambda_gradU(i))==0
                 lambdaci(ix,iz,iy)=abs(imag(lambda_gradU(i)));
             end
         end
       end
    end
end

% Q : 8000
% 

[X,Y,Z]=meshgrid((1:m)*deltaX,(1:n)*deltaZ,y(2:round(p/2)+1));

% Q
figure(1)
hold on
isosurface(X,Y,Z,Q(:,:,1:round(p/2)),0.05e4);
camlight; camlight(-80,-10); lighting phong;
h = slice(X,Y,Z,real(squeeze(data2(1,2:end-1,2:end-1,1:round(p/2)))),[Lx/2],[],[]);
set(h,'FaceColor','interp',...
 'EdgeColor','none')


% Delta
figure(3)
hold on
isosurface(X,Y,Z,Delta(:,:,1:round(p/2)),0.005e10);
camlight; camlight(-80,-10); lighting phong;
h = slice(X,Y,Z,real(squeeze(data2(1,2:end-1,2:end-1,1:round(p/2)))),[Lx/2],[],[]);
set(h,'FaceColor','interp',...
 'EdgeColor','none')

% lambdaCI
figure(4)
hold on
isosurface(X,Y,Z,lambdaci(:,:,1:round(p/2)),15);
camlight; camlight(-80,-10); lighting phong;
h = slice(X,Y,Z,real(squeeze(data2(1,2:end-1,2:end-1,1:round(p/2)))),[Lx/2],[],[]);
set(h,'FaceColor','interp',...
 'EdgeColor','none')

% lambda2
figure(2)
hold on
isosurface(X,Y,Z,lambda2(:,:,1:round(p/2)),-0.1e4);
camlight; camlight(-80,-10); lighting phong;
h = slice(X,Y,Z,real(squeeze(data2(1,2:end-1,2:end-1,1:round(p/2)))),[Lx/2],[],[]);
set(h,'FaceColor','interp',...
 'EdgeColor','none')