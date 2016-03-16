% ------------ function signal_smoothed=smooth(signal, width) -------------
% Purpose :
% smooth a 2D signal (obtaiend from max/min_filter in FABEMD) over a indow
% of width w:the value at (x,y) is equal to the average of all the values 
% in the window of width w centered at (x,y)
% Input :
% signal : 2D matrix of doubles of size (m,n)
% width: integer, width of the windowing filter
% Output :
% signal_smoothed : 2D matrix of size (m,n)



function signal_smoothed=smooth(signal, width)

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
        signal_smoothed(i,j)=sum(sum(signal(i:i+2*w,j:j+2*w)))/(width*width);
    end
end