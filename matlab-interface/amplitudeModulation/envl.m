function f=envl(f,sig,cutoff,freq)
% Upper (sig=1) or lower (sig=-1) envelope of a small-scale real signal
% f, sampled at a frequency freq and obtained by a 6th-order, delay-coorected 
% low-pass filter with cut-off frequency cutoff.
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
f=2*sig*lowpass(f,cutoff,freq,'ImpulseResponse','iir','Steepness',0.95,'StopbandAttenuation',100);
