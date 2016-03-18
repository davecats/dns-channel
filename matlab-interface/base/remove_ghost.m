function [U,Ubar,Wbar]=remove_ghost(U,Ubar,Wbar)

% Löschung der ersten und letzten y Daten
U=real(squeeze(U(:,:,:,2:size(U,4)-1)));
Ubar=Ubar(2:end-1);
Wbar=Wbar(2:end-1);

