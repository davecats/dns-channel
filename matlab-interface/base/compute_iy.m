% ----------------function iy=compute_iy(y_plus,Re,ny)---------------------
% Purpose : 
% compute the iyth node corresponding to the y_plus value for the
% Reynolds number Re and ny+1 nodes

function iy=compute_iy(y_plus,Re,ny)

a=1.66;
ymax=2;
ymin=0;
y=y_plus/Re;
iy=round(ny/2*(atanh(tanh(a)*((y-ymin)/(0.5*(ymax-ymin))-0.5*(ymax-ymin)))/a+1));