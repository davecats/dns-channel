% This function reads a field (v-eta) produced by dns-channel and saves it 
% into a Matlab array. It can read both a whole field or a single plane. 

% --------------------------------- Input ---------------------------------
% file : Name der zu lesenen Datei [char]
% y (optional) : Ebene zu lesen, y=0:ny [int]

% --------------------------------- Output --------------------------------
% dns : Parameter der Simulation
%     dns.ny : Anzahl der Knoten in y-Richtung [Integer]
%     dns.nx : Anzahl der Knoten in x-Richtung [Integer]
%     dns.nz : Anzahl der Knoten in z-Richtung [Integer]
%     dns.time : time instant at which the field was saved
%     dns.ymin : Koordinate der unteren Wand
%     dns.ymax : Koordinate der oberen Wand
%     dns.a : stretching parameter (used in the gittering in the y 
%           direction)
%     dns.alfa0 : base wavenumber in x. The length of the computational 
%           box is Lx=2*pi/alfa0
%     dns.beta0 : base wavenumber in z. The width of the computational 
%           box is Lz=2*pi/beta0
%     dns.ni : kinematic viscosity. If one works nondimensional this is 
%           1/Re, where Re is a Reynolds number
% U : gemittelte Geschwindigkeit in x-Richtung für y=-1:ny+1 (real)
% W : gemittelte Geschwindigkeit in z-Richtung für y=-1:ny+1 (real)
% VEta : 4D Matrix von double complex
%     1. Dim: nichts (null, will be used later to save the u component of the speed),Geschwindigkeit v, wandnormale Wirbelstärke Eta
%     2. Dim: Z=-nz:nz
%     3. Dim: X=0:nx
%     4. Dim: Y=-1:ny+1 OR Y=y-1:y+1

% -------------- Reihenfolge der Daten in der binären Datei ---------------
% dns : 3 Integer (int32), 7 double
%    /\    Datei in 64-bits Format >> nach den Integer ist der nächste Wert  
%   /__\   erst ab nächstem Mehrfach von 64 bits gespeichert
% U,W : jeweils ny+3 doubles
% VEta : 2*2*(2*nz+1)*(nx+1)*(ny+3) doubles
%     Zuordnung: Real / Imag / V / Eta / z / x / y

% =========================================================================

function[dns,VDisk]=open_cart(file)

% Open file
fileID=fopen(file);
% Read header
header=fread(fileID,1024,'*char')';
i=find(header=='=');
dns.ny=str2double(header(i(1)+1 : i(2)-3));
dns.nx=str2double(header(i(2)+1 : i(3)-3));
dns.nz=str2double(header(i(3)+1 : i(4)-6));
dns.alfa0=str2double(header(i(4)+1 : i(5)-6));
dns.beta0=str2double(header(i(5)+1 : i(6)-5));
dns.ymin=str2double(header(i(6)+1 : i(7)-5));
dns.ymax=str2double(header(i(7)+1 : i(8)-2));
dns.a=str2double(header(i(8)+1 : i(9)-3));
dns.Re=str2double(header(i(9)+1 : i(10)-5));
dns.time=str2double(header(i(10)+1 : end));
offset=ftell(fileID);
fclose(fileID);
VDisk=memmapfile(file, ...
                'Format', {'double' [2 3 2*dns.nz+1 dns.nx+1 dns.ny+3] 'V'}, ...
                'Offset', 1024);
