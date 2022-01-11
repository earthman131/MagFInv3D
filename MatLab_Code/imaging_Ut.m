function mcal = imaging_Ut(Ut,u0,M,F,n,N,nx,ny,nz,dz,z,npts,k,kx,ky,pressfilter)
Ut=Ut*1e-9; 
ydiff=floor((npts-ny)/2); xdiff=floor((npts-nx)/2); 
Utex=extend_copy2d(Ut,nx,ny,nz,npts);
Utf=fftshift(fft2(Utex));
mf=zeros(npts,npts,nz);
mcal=zeros(ny,nx,nz);
Cm=1i*M(1)*kx+1i*M(2)*ky+M(3)*k;
Cf=1i*F(1)*kx+1i*F(2)*ky+F(3)*k;
for K=1:nz
    zo=z(K)+dz/2;
    imgf=((2*k+eps)./(u0*Cf.*Cm+eps))*((n+1)^(n+1)/(factorial(n)))*(zo/N)^n.*k.^(n+1).*exp(-k*n*zo/N);
    mf(:,:,K)=Utf.*imgf;
end
for K=1:nz
    mf(:,:,K)=ifft2(ifftshift(mf(:,:,K)));
    mcal(:,:,K)=real(mf((ydiff+1):(ydiff+ny),(xdiff+1):(xdiff+nx),K));
end