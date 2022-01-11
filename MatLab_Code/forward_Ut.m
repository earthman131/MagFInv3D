function Ut = forward_Ut(m,u0,M,F,nx,ny,nz,dz,z,npts,k,kx,ky)
Utf=zeros(npts,npts,nz);anof=zeros(npts,npts,nz);
anoex=zeros(npts,npts);ano=zeros(ny,nx,nz);Ut=zeros(ny,nx);
ydiff=floor((npts-ny)/2); xdiff=floor((npts-nx)/2); 
mex=extend_copy2d(m,nx,ny,nz,npts);
for i = 1:nz
    Utf(:,:,i)=fftshift(fft2(mex(:,:,i)));
end
Cm=1i*M(1)*kx+1i*M(2)*ky+M(3)*k;
Cf=1i*F(1)*kx+1i*F(2)*ky+F(3)*k;
for K=1:nz
    h1=z(K); 
    h2=z(K)+dz;
    anof(:,:,K)=u0*Utf(:,:,K).*((Cm.*Cf+eps)./(2*k.^2+eps)).*(exp(-k*h1)-exp(-k*h2));
end
for K=1:nz
    anoex(:,:,K)=ifft2(ifftshift(anof(:,:,K)));
    ano(:,:,K)=real(anoex((ydiff+1):(ydiff+ny),(xdiff+1):(xdiff+nx),K))*1e9;
    Ut(:,:)=Ut(:,:)+(ano(:,:,K));
end
