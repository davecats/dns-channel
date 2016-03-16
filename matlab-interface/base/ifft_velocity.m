% function [data]=ifft_velocity(data,i,nx,ny,nz)

% --------------------------------- Input ---------------------------------
% data : matrix(3,2*nz+1,nx+1,ny+3), FFT der v-Komponente
% i = 1, 2 oder 3 : zeigt welche geschwindigkeitskomponente betrachtet wird
% 1=u, 2=v,3=w
% nx,ny,nz : number of points 

% --------------------------------- Output --------------------------------
% data : Matrix, so groß wie data(Input), enthält die reale Geschwindigkeit
% 8! komplexe Zahle >> imaginärer Teil vernachlässigbar (Rundungsfehler)

function [data]=ifft_velocity(data,i)

nx=(size(data,3)-1)/2;
ny=size(data,4)-3;
nz=(size(data,2)-1)/2;

parfor iy=1:ny+1
    
    data(i,:,:,iy+1)=plane_ifft(data(i,:,:,iy+1),nx,nz);
    
end