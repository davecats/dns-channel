function [PSD]=plane_psd(U1,U2)

% Remove mean (is this really wanted?)
U1(floor(end/2)+1,floor(end/2)+1,:)=0;
U2(floor(end/2)+1,floor(end/2)+1,:)=0;

% Compute PSD
PSD=real(conj(U2).*U1);
