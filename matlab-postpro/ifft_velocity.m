% function [data]=ifft_velocity(data,i,nx,ny,nz)


% --------------------------------- Input ---------------------------------

% data : matrix(3,2*nz+1,nx+1,ny+3), FFT der v-Komponente
% i = 1, 2 oder 3 : ebene der ersten Komponente
% nx,ny,nz : number of points 


% --------------------------------- Output --------------------------------

% data : Matrix, so groß wie data(Input), enthält die reale Geschwindigkeit
% 8! komplexe Zahle >> imaginärer Teil vernachlässigbar (Rundungsfehler)


function [data]=ifft_velocity(data,i,nx,ny,nz)

for iy=1:ny+1
    datay(1:2*nz+1,nx+1:2*nx+1)=squeeze(data(i,1:2*nz+1,nx+1:2*nx+1,iy+1));
    for ix=nx+1:2*nx+1 datay(:,ix)=ifft(ifftshift(datay(:,ix)))*(2*nz+1); end
    for iz=1:2*nz+1
        datay(iz,1:nx)=conj(fliplr(datay(iz,nx+2:2*nx+1)));
        datay(iz,:)=ifft(ifftshift(datay(iz,:)))*(2*nx+1);
    end
    data(i,:,:,iy+1)=datay;
end