function fhat = ftransform(gridf,Nz,U)
    U(:,:,1:Nz)=gridf;
    U(:,:,Nz+1:2*Nz-2)=gridf(:,:,Nz-1:-1:2);
    U = fftn(U);
    fhat=zeros(size(gridf));
    fhat(:,:,1)=U(:,:,1);
    fhat(:,:,2:Nz-1)=U(:,:,2:Nz-1)+U(:,:,end:-1:Nz+1);
    fhat(:,:,Nz)=U(:,:,Nz);
    fhat=fhat/(2*Nz-2);
end
