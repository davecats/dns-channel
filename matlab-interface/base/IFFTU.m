function V = IFFTU(V)

disp('Transforming to physical space')
V=ifft_velocity(V,1);
V=ifft_velocity(V,2);
V=ifft_velocity(V,3);