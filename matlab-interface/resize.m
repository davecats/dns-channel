%
% This program resizes a Dati.cart field 
% to a given newer (smaller or larger) size
%
% The initial header of Dati.cart is not modified
% and must be adjusted manually (implies you can not
% use time_from_restart=YES) 
%

% This program is NOT parallel

% INPUT
% ------------------------------
infile = './base/Dati.cart.out';
outfile = './interp.out';
newdns.nx = 188;
newdns.ny = 150;
newdns.nz = 188;
newdns.a = 1.6;
Uscale = 2 / (4/3);
method = 'linear';
% ------------------------------

clc
addpath ./base

% Read input field
[olddns,oldV]=open_cart(infile);
% Create file for output field
newV=storedstructure(outfile, {'uint8'  1024 'header';
                               'double' [2 3 2*newdns.nz+1 newdns.nx+1 newdns.ny+3] 'V'}, true);
% Create y-grid for input field
oldy=zeros(olddns.ny+3,1);
for i=1:olddns.ny+3
    oldy(i)=olddns.ymin+0.5*(olddns.ymax-olddns.ymin)* ...
            (tanh(olddns.a*(2*(i-2)/olddns.ny-1))/tanh(olddns.a)...
            +0.5*(olddns.ymax-olddns.ymin));
end
% Create y-grid for output field
newdns.ymin=olddns.ymin;  newdns.ymax=olddns.ymax;
newy=zeros(newdns.ny+3,1);
for i=1:newdns.ny+3
    newy(i)=newdns.ymin+0.5*(newdns.ymax-newdns.ymin)* ...
            (tanh(newdns.a*(2*(i-2)/newdns.ny-1))/tanh(newdns.a)...
            +0.5*(newdns.ymax-newdns.ymin));
end
% Interpolate fields        
for ix=0:min(newdns.nx,olddns.nx)
   disp([num2str(round(100*ix/min(newdns.nx,olddns.nx))) '  % completed'])
   IX=ix+1;
   for iz=-min(newdns.nz,olddns.nz):min(newdns.nz,olddns.nz)
        oldIZ=iz+olddns.nz+1; newIZ=iz+newdns.nz+1;
        for iV=1:3
            for i=1:2     
              newV.data.V(i,iV,newIZ,IX,1:newdns.ny+3)= ...
                             interp1( oldy, ...
                                      squeeze(oldV.data.V(i,iV,oldIZ,IX,1:olddns.ny+3)), ...
                                      newy, ...
                                      method, ...
                                      'extrap')*Uscale;
            end
        end
    end
end
% Read header from old field
fileID=fopen(infile);  header=uint8(fread(fileID,1024,'uint8')); fclose(fileID);
% Write header to new field
newV.data.header=header;