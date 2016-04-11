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

function [U]=veta2uvwIY(U,Ubar,Wbar,y,dns)

%   2) (u,w) aus (v,Eta) berechnen
% -------------------------------------------------------------------------
U(1,:,:,2)=(U(2,:,:,3)-U(2,:,:,1))./(y(3)-y(1)); %dvdy
for ix=1:dns.nx+1
    alfa=(ix-1)*dns.alfa0;
    for iz=1:2*dns.nz+1
        beta=(iz-1-dns.nz)*dns.beta0; k2=alfa^2+beta^2;
        tmp = 1j*(alfa*U(1,iz,dns.nx+ix,2)-beta*U(3,iz,dns.nx+ix,2))/k2;
        U(3,iz,dns.nx+ix,2)=1j*(beta*U(1,iz,dns.nx+ix,2)+alfa*U(3,iz,dns.nx+ix,2))/k2;
        U(1,iz,dns.nx+ix,2)=tmp;
    end
end
U=U(:,:,:,2);
U(1,dns.nz+1,dns.nx+1)=Ubar(2);
U(3,dns.nz+1,dns.nx+1)=Wbar(2);

