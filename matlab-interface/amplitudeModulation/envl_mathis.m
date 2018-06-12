function f=envl_mathis(f,cutoff,freq)
% Envelope of a small-scale real signal
% f, sampled at a frequency freq and obtained by a 6th-order, delay-coorected 
% low-pass filter with cut-off frequency cutoff.
%
% If f is a matrix, then the function operates along the columns of f.
%
% Method as outlined in Figure 4 of Mathis, Hutchins & Marusic, J. Fluid
% Mech.,628, 311-337 (2009)


% Augmented signal
f=hilbert(f);
% Modulus of augmented signal
f=abs(f);
% Low-pass filter the signal
f=lowpass(f,cutoff,freq,'ImpulseResponse','iir','Steepness',0.95,'StopbandAttenuation',100);
