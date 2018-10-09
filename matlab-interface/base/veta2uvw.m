function field=veta2uvw(field,params,derivatives)
% Converts from (v,eta) -> (u,v,w)
% not memory efficient but fast, use only for small fields

ny=params.ny;
nx=params.nx;
nz=params.nz;

vy = reshape((derivatives.d0\(derivatives.d1*reshape(field.veta(1,:,:,:),[(nx+1)*(2*nz+1), ny+3])'))',[1,2*nz+1,nx+1,ny+3]);

field.uvw = zeros(3,2*nz+1,nx+1,ny+3);
alfa=repmat(reshape((0:nx),[1,1,nx+1,1])*params.alfa0,[1,2*nz+1,1,ny+3]); 
beta=repmat(reshape((-nz:nz),[1,2*nz+1,1,1])*params.alfa0,[1,1,nx+1,ny+3]); 
k2=alfa.^2 + beta.^2;
field.uvw(1,:,:,:) = 1j*(alfa.*vy-beta.*field.veta(2,:,:,:))./k2; % u
field.uvw(2,:,:,:) = field.veta(2,:,:,:); % v
field.uvw(3,:,:,:) = 1j*(beta.*vy+alfa.*field.veta(2,:,:,:))./k2; % w

field.uvw(1,nz+1,1,:)=field.U;
field.uvw(3,nz+1,1,:)=field.W;