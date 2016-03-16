% ------function [signal_max,signal_min]=extrema_filter(signal,width)------

% Purpose :
% filtered_signal(x,y)=max(signal(s,t)) with (s,t)belongsto(Z), Z being a
% square of size width*width centered at (x,y)
% Input :
% signal: 2D matrix
% width: odd integer
% Output :
% filtered:signal: 2D matrix (same size as signal)

% to consider the periodicity of the input signal, a periodicity condition
% is added :
% if i>size(signal,1),  i=i-size(signal,1)
% if i<1,               i=i+size(signal,1)
% if j>size(signal,2),  j=j-size(signal,2)
% if j<1,               j=j+size(signal,2)

function [signal_max,signal_min]=extrema_filter(signal,width)

signal_max=zeros(size(signal));
signal_min=zeros(size(signal));
w=(width-1)/2;

[m,n]=size(signal);
signal(w+1:w+m,:)=signal(:,:);
signal(1:w,:)=signal(m+1:m+w,:);
signal(m+w+1:m+2*w,:)=signal(w+1:2*w,:);
signal(:,w+1:w+n)=signal(:,:);
signal(:,1:w)=signal(:,n+1:n+w);
signal(:,n+w+1:n+2*w)=signal(:,w+1:2*w);

for i=1:m
    for j=1:n
        signal_max(i,j)=max(max(signal(i:i+2*w,j:j+2*w)));
        signal_min(i,j)=min(min(signal(i:i+2*w,j:j+2*w)));
    end
end