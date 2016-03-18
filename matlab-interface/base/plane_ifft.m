function [f]=plane_ifft(f,nx,nz)
    
    parfor ix=nx+1:2*nx+1
        f(1,:,ix)=ifft(ifftshift(f(1,:,ix)))*(2*nz+1); 
    end
    f(1,:,:) = cat(3,f(1,:,nx+1:end),conj(flip(f(1,:,nx+2:2*nx+1),3)));
    parfor iz=1:2*nz+1
        f(1,iz,:)=ifft( f(1,iz,:) )*(2*nx+1);
    end
    
end