function f=envl(f,sig,filterfun,filterpars)
% Upper (sig=1) or lower (sig=-1) envelope of a small-scale real signal
% f, low-pass filtered by the custom filter @filterfun with parameters
% @filterpars. 
%
% If f is a matrix, then the function operates along the columns of f.
%
% Method as outlined in Eq. (24) and (25) of Agostini, Leschziner &
% Gaitonde, Phys. Fluids 28, 015110 (2016), also described in Figure 16 of
% the same paper. 


% Take only positive (negative) part of signal
f=sig*((sig*f)>0);
% Augmented signal
f=hilbert(f);
% Modulus of augmented signal
f=abs(f);
% Low-pass filter the signal
f=2*sig*filterfun(f,filterpars{:});
