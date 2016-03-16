% ----------function [emdc]=save_FABEMD(iy,parameter,save_folder)----------

% Purposes :
%   - load a plane y=const from a simulation file sim (*.fld)
%   - run FABEMD(sim) with user defines parameters for FABEMD (type, MNAI,
%   number of iterations to be computed) and register some parameters of
%   the simulation in a text file

% ------------------INPUT----------------------
% iy : node of the plane to be read // y_plus=adimensional distance from the wall
% parameter: structure, parameters of the DNS and the FABEMD
% save_folder: the path where the data should be saved

% ------------------OUTPUT--------------------
% emdc: matrix(2*nz+1,2*nx+1,n+1), contains the n modes + the residual

function emdc=save_FABEMD(iy,parameter,save_folder)

tic;
path='/net/istmuranus/volume1/data/hi221/BACKUP-WORK-HLR1/CONDITIONAL-AVERAGE/';

% -------------------------------------------------------------------------
% Go to the right folder 
% -------------------------------------------------------------------------

old_folder=cd(save_folder);

if exist(['Re' num2str(parameter.Re)],'dir')==7
    cd(['Re' num2str(parameter.Re)]);
else
    mkdir(['Re' num2str(parameter.Re)]);
    cd(['Re' num2str(parameter.Re)]);
end
if exist(parameter.box,'dir')==7
    cd(parameter.box);
else
    mkdir(parameter.box);
    cd(parameter.box);
    actual_folder=cd(old_folder);
    save_info(parameter);
    cd(actual_folder);
end
folder=struct('folder',{});
folder(1).folder=parameter.version;
folder(2).folder=['Field' num2str(parameter.fieldnb)];
folder(3).folder=parameter.type;
for i=1:3
    if exist(folder(i).folder,'dir')==7
        cd(folder(i).folder);
    else
        mkdir(folder(i).folder);
        cd(folder(i).folder);
    end
end

% -------------------------------------------------------------------------
% compute FABEMD
% -------------------------------------------------------------------------
if exist([parameter.velocity '_' num2str(iy) '.bin'],'file')==2
    cd(old_folder);
    emdc=read_emdc_bin(parameter,iy,save_folder);
else 
    if exist([parameter.velocity '.txt'],'file')==2
        fid=fopen([parameter.velocity '.txt'],'a');
    else fid=fopen([parameter.velocity '.txt'],'w');
        fprintf(fid,'iy\tmode\tMNAI\tw\tnb min\tnb max\ttime(s)\tmean\ttmin\ttdmin\ttfilter\ttsmin\t\n');
    end
    new_folder=cd([save_folder '/Re' num2str(parameter.Re) '/' parameter.box '/' parameter.version '/Field' num2str(parameter.fieldnb)]);
    if exist(['velocity_' num2str(iy) '.bin'],'file')==2
        cd(old_folder);
        [sim,dns]=read_vel_bin(parameter,iy,save_folder);
    else 
        file=[path '/Re' num2str(parameter.Re) '/' parameter.box '/' parameter.version '/postpro-fld/Field' num2str(parameter.fieldnb) '.fld'];
        actual_folder=cd(old_folder);
        [sim,dns]=veta2uvw(file,iy);
        cd(actual_folder);
        fid_vel_bin=fopen(['velocity_' num2str(iy) '.bin'],'w');
        fwrite(fid_vel_bin,[dns.nx dns.nz],'int32'); 
        fwrite(fid_vel_bin,reshape(sim,[1 3*(2*dns.nx+1)*(2*dns.nz+1)]),'double');
        fclose(fid_vel_bin);
        cd(old_folder);
    end
    switch parameter.velocity
        case 'u'
            vel=1;
        case 'v'
            vel=2;
        case 'w'
            vel=3;
    end
    emdc=FABEMD(squeeze(sim(vel,:,:)),parameter,fid,iy);
    fclose(fid);
    cd(new_folder);
    fid_bin=fopen([parameter.velocity '_' num2str(iy) '.bin'],'w');
    fwrite(fid_bin,[dns.nx dns.nz parameter.n],'int32'); 
    fwrite(fid_bin,reshape(emdc,[1 (parameter.n+1)*size(sim,2)*size(sim,3)]),'double');
    fclose(fid_bin);
    cd(old_folder);
end

computation_time_total=toc