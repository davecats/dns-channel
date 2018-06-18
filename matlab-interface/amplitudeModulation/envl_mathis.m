function f=envl_mathis(f,filterfun,filterpars)
% Envelope of a small-scale real signal
% f, low-pass filtered by the custom filter @filterfun with parameters
% @filterpars. 
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
f=filterfun(f,filterpars{:});
