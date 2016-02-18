% --------------- Programm zur Auswertung von *.fld Dateien ---------------

% Ziel :    Geschwindigkeitskomponenten (V,Eta) aus einer spectral-
%           Diskretisierung in reale Geschwindigkeiten (u,v,w) in einem 
%           spatialen Gitter umrechnen
 
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

function [data2,dns]=auswertung_fld()

clear all;
% close all;

% beginn of a chronometer for the improvement of the code
tic;

%   1) Rohdaten speichern
% -------------------------------------------------------------------------

[dns,U,W,data2]=import_fld('Field10.fld');

%   2) (u,w) aus (v,Eta) berechnen
% -------------------------------------------------------------------------

y=zeros(dns.ny+3);
for i=1:dns.ny+3
    y(i)=dns.ymin+0.5*(dns.ymax-dns.ymin)*(tanh(dns.a*(2*(i-2)/dns.ny-1))/tanh(dns.a)+0.5*(dns.ymax-dns.ymin));
end

dvdy=zentral_diff(squeeze(data2(2,:,dns.nx+1:2*dns.nx+1,:)),y,dns.nx,dns.ny,dns.nz);


uy=zeros(2*dns.nz+1,dns.nx+1);
wy=zeros(2*dns.nz+1,dns.nx+1);
for iy=2:dns.ny+2
    for ix=1:dns.nx+1
        alfa=(ix-1)*dns.alfa0;
        for iz=1:2*dns.nz+1
            beta=(iz-1-dns.nz)*dns.beta0; k2=alfa^2+beta^2;
            uy(iz,ix)=1j*(alfa*dvdy(iz,ix,iy)-beta*data2(3,iz,dns.nx+ix,iy))/k2;
            wy(iz,ix)=1j*(beta*dvdy(iz,ix,iy)+alfa*data2(3,iz,dns.nx+ix,iy))/k2;
        end
    end
    data2(1,1:2*dns.nz+1,dns.nx+1:2*dns.nx+1,iy)=uy;
    data2(3,1:2*dns.nz+1,dns.nx+1:2*dns.nx+1,iy)=wy;
end
data2(1,dns.nz+1,dns.nx+1,2:dns.ny+2)=U(2:dns.ny+2);
data2(3,dns.nz+1,dns.nx+1,2:dns.ny+2)=W(2:dns.ny+2);
time_seconds=toc

%   3) inverse FFT : spectral >> spatial
% -------------------------------------------------------------------------

data2=ifft_velocity(data2,1,dns.nx,dns.ny,dns.nz);
data2=ifft_velocity(data2,2,dns.nx,dns.ny,dns.nz);
data2=ifft_velocity(data2,3,dns.nx,dns.ny,dns.nz);

%end of the chronometer, memory used at the end of the programm
time_seconds=toc