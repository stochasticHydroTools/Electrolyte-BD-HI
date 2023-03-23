function fgrid = btransform(fhat,Nz,U,varargin)
    if nargin == 3
        U(:,:,1)=fhat(:,:,1);
        U(:,:,2:Nz-1)=fhat(:,:,2:Nz-1)/2;
        U(:,:,Nz)=fhat(:,:,Nz);
        U(:,:,Nz+1:2*Nz-2)=fhat(:,:,Nz-1:-1:2)/2;
        U=U*(2*Nz-2);
        U = real(ifftn(U));
        fgrid = U(:,:,1:Nz);
    else
        % go from cheb coeffs to function values on uniform grid in z
        Tu = varargin{1}; [Nzu,~] = size(Tu);
        [nxy,~,~] = size(fhat);
        fhatu = zeros(nxy,nxy,Nzu);
        for iX = 1:nxy
            for iY = 1:nxy
             cf = reshape(fhat(iY,iX,:),Nz,1);
             fhatu(iY,iX,:) = Tu*cf;
            end
        end 
        % backtransform in x,y
        fgrid = zeros(nxy,nxy,Nzu);
        for iZ = 1:Nzu
            fgrid(:,:,iZ) = real(ifft2(reshape(fhatu(:,:,iZ),nxy,nxy)));
        end
    end
end

