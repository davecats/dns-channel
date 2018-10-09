function [field,params]=readfield(filename,params)
% reads a DNS field in the (u,v,w) format

f=fopen(filename);
% params.ny=fread(f,1,'uint32');
% params.nx=fread(f,1,'uint32');
% params.nz=fread(f,1,'uint32');
% fread(f,1,'uint32');
% params.time=fread(f,1,'double');
% params.ymin=fread(f,1,'double');
% params.ymax=fread(f,1,'double');
% params.a=fread(f,1,'double');
% params.alfa0=fread(f,1,'double');
% params.beta0=fread(f,1,'double');
% params.Re=fread(f,1,'double');
% field.U=fread(f,params.ny+3,'double');
% field.W=fread(f,params.ny+3,'double');
% params.y = tanh( params.a * ( 2 * ( (-1:params.ny+1) )/params.ny - 1) )/tanh( params.a ) + 1;
field.U=reshape(fread(f,'double'),[3,params.ny+3,params.nzd,params.nxd]);
% field.veta=complex(squeeze(field.veta(1,:,:,:,:)),squeeze(field.veta(2,:,:,:,:)));
fclose(f);
