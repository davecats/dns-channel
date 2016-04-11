function dudy=FD(EMDC,y,order)

% INPUT
% EMDC=(z,x,ymin:ymax)
% y=distance between the points for ymin:ymax
% order: order of the approximation, from 1 to 4, order=ymax-ymin !!

if order==1
    a=-1/(y(2)-y(1));
    b=1/(y(2)-y(1));
    dudy=a*squeeze(EMDC(:,:,1))+b*squeeze(EMDC(:,:,2));
elseif order==2
    a=-(y(3)+y(2))/(y(3)*y(2));
    b=y(3)/(y(2)*(y(3)-y(2)));
    c=y(2)/(y(3)*(y(2)-y(3)));
    dudy=a*squeeze(EMDC(:,:,1))+b*squeeze(EMDC(:,:,2))...
        +c*squeeze(EMDC(:,:,3));
elseif order==3
    b=y(3)*y(4)/(y(2)*(y(2)-y(3))*(y(2)-y(4)));
    c=y(2)*y(4)/(y(3)*(y(3)-y(2))*(y(3)-y(4)));
    d=y(2)*y(3)/(y(4)*(y(4)-y(3))*(y(4)-y(2)));
    a=-b-c-d;
%     eq1=a+b+c+d
%     eq2=b*y(2)+c*y(3)+d*y(4)
%     eq3=b*y(2)^2+c*y(3)^2+d*y(4)^2
%     eq4=b*y(2)^3+c*y(3)^3+d*y(4)^3
    dudy=a*squeeze(EMDC(:,:,1))+b*squeeze(EMDC(:,:,2))...
        +c*squeeze(EMDC(:,:,3))+d*squeeze(EMDC(:,:,4));
elseif order==4
    b=-y(3)*y(4)*y(5)/(y(2)*(y(2)-y(3))*(y(2)-y(4))*(y(2)-y(5)));
    c=-y(2)*y(4)*y(5)/(y(3)*(y(3)-y(2))*(y(3)-y(4))*(y(3)-y(5)));
    d=-y(2)*y(3)*y(5)/(y(4)*(y(4)-y(3))*(y(4)-y(2))*(y(4)-y(5)));
    e=-y(2)*y(3)*y(4)/(y(5)*(y(5)-y(3))*(y(5)-y(2))*(y(5)-y(4)));
    a=-b-c-d-e;
%     eq1=a+b+c+d+e
%     eq2=b*y(2)+c*y(3)+d*y(4)+e*y(5)
%     eq3=b*y(2)^2+c*y(3)^2+d*y(4)^2+e*y(5)^2
%     eq4=b*y(2)^3+c*y(3)^3+d*y(4)^3+e*y(5)^3
%     eq5=b*y(2)^4+c*y(3)^4+d*y(4)^4+e*y(5)^4
    dudy=a*squeeze(EMDC(:,:,1))+b*squeeze(EMDC(:,:,2))...
        +c*squeeze(EMDC(:,:,3))+d*squeeze(EMDC(:,:,4))...
        +e*squeeze(EMDC(:,:,5));
end