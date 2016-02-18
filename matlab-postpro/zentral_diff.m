% function dvdy=zentral_diff(data,y,nx,ny,nz)


% --------------------------------- Input ---------------------------------

% data : matrix(2*nz+1,nx+1,ny+3), FFT der v-Komponente
% y : Vektor(1:nx+3), position of the points in the y direction, y(2)=ymin
% and y(ny+2)=ymax
% nx,ny,nz : number of points 


% --------------------------------- Output --------------------------------

% dvdy : matrix(2*nz+1,nx+1,ny+2), dvdy(iy,ix,2:2*ny+2) contains the 
% central differences of v along y


function dvdy=zentral_diff(data,y,nx,ny,nz)

dvdy=zeros(2*nz+1,2*nx+1,ny+2);
for iy=2:ny+2
            dvdy(:,:,iy)=(data(:,:,iy+1)-data(:,:,iy-1))/(y(iy+1)-y(iy-1)); 
end

