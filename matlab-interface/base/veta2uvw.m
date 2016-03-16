% This function imports a field produced with dns-channel and 
% outputs the velocity field converted in (u,v,w). It can read the whole
% field or a single plane. 

% --------------- Programm zur Auswertung von *.fld Dateien ---------------
%
% Ziel :    Geschwindigkeitskomponenten (V,Eta) aus einer spectral-
%           Diskretisierung in reale Geschwindigkeiten (u,v,w) in einem 
%           spatialen Gitter umrechnen
%           und eine besondere Ebene y=const lesen (y=0:ny)
 
% Etapen des Programms :
%   1) die Rohdaten der binären Datei durch die Funktion 'importieren' 
%       speichern
%   2) (u,v,w) aus (v,Eta) berechnen (FFT-spectrum) :
%           u=1/k²*j(ix*kx*[dv/dy]-iz*kz*Eta)
%           w=1/k²*j(ix*kx*Eta+iz*kz*[dv/dy])
%       mit k²=ix²*kx²+iz²*kz² und dvdy die Ableitung von v in der y-Richt
%       (hier wird eine zentrale Differenzierung benutzt)
%        /\     Der erste Modus (z=0,x=0) ist mit dieser Berechnung nicht
%       /__\    beinhaltet, der muss als (U,W) dazu addiert werden.
%   3) die in FFT- representation Daten in realen daten umrechnen

% =========================================================================

function [data2,U,W,dns]=veta2uvw(filename)

%   1) Rohdaten speichern
% -------------------------------------------------------------------------
[dns,U,W,data2]=import_fld(filename);

%   2) (u,w) aus (v,Eta) berechnen
% -------------------------------------------------------------------------
y=zeros(dns.ny+3,1);
for i=1:dns.ny+3
    y(i)=dns.ymin+0.5*(dns.ymax-dns.ymin)*(tanh(dns.a*(2*(i-2)/dns.ny-1))/tanh(dns.a)+0.5*(dns.ymax-dns.ymin));
end
disp('Converting (v,eta) --> (u,v,w)')
for iy=2:dns.ny+2
  data2(1,:,:,iy)=(data2(2,:,:,iy+1)-data2(2,:,:,iy-1))./(y(iy+1)-y(iy-1)); %dvdy
end
for ix=1:dns.nx+1
    alfa=(ix-1)*dns.alfa0;
    for iz=1:2*dns.nz+1
        beta=(iz-1-dns.nz)*dns.beta0; k2=alfa^2+beta^2;
        for iy = 2:dns.ny+2
          tmp = 1j*(alfa*data2(1,iz,dns.nx+ix,iy)-beta*data2(3,iz,dns.nx+ix,iy))/k2;
          data2(3,iz,dns.nx+ix,iy)=1j*(beta*data2(1,iz,dns.nx+ix,iy)+alfa*data2(3,iz,dns.nx+ix,iy))/k2;
          data2(1,iz,dns.nx+ix,iy)=tmp;
        end
    end
end
data2(1,dns.nz+1,dns.nx+1,:)=U(:);
data2(3,dns.nz+1,dns.nx+1,:)=W(:);


%   3) inverse FFT : spectral >> spatial
% -------------------------------------------------------------------------
disp('Transforming to physical space')
data2=ifft_velocity(data2,1);
data2=ifft_velocity(data2,2);
data2=ifft_velocity(data2,3);

% Löschung der ersten und letzten y Daten
data2=real(squeeze(data2(:,:,:,2:size(data2,4)-1)));
U=U(2:end-1);
W=W(2:end-1);
