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

function[dns,U,W,VEta]=import_fld(file)

% Open file
fileID=fopen(file);
% Read header
dns=struct('ny',fread(fileID,1,'int32'),'nx',fread(fileID,1,'int32'),'nz',fread(fileID,1,'int32'));
fseek(fileID,32/8,'cof');
dns.time=fread(fileID,1,'double');
dns.ymin=fread(fileID,1,'double');
dns.ymax=fread(fileID,1,'double');
dns.a=fread(fileID,1,'double');
dns.alfa0=fread(fileID,1,'double');
dns.beta0=fread(fileID,1,'double');
dns.ni=fread(fileID,1,'double');
% Read mean velocity profiles
U=fread(fileID,(dns.ny+3),'double');
W=fread(fileID,(dns.ny+3),'double');
offset=ftell(fileID);
% Declare the velocity field
VEta=complex(zeros(3,2*dns.nz+1,2*dns.nx+1,dns.ny+3));
fclose(fileID);
% Map file to memory
disp(strcat('Reading velocity field : ', file))
VDisk=memmapfile(file, ...
                 'Format', {'double' [2 2 2*dns.nz+1 dns.nx+1 dns.ny+3] 'V'}, ...
                 'Offset', offset);
VEta(2:3,:,dns.nx+1:2*dns.nx+1,:)=complex(VDisk.data.V(1,:,:,:,:),...
                                           VDisk.data.V(2,:,:,:,:));


