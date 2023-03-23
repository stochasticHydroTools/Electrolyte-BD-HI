% Computes the spread and interpolation matrices by fast Gaussian gridding.
% Vectorized and linear in the number of points. 
function [SWeights,IWeights] = SpreadWts(xEpts,yEpts,zEpts,IBpts,sup,gw,iwtsz)
    Nx=length(xEpts);
    Ny=length(yEpts);
    Nz=length(zEpts);
    [NIB,~]=size(IBpts);
    hex=xEpts(2)-xEpts(1);
    hey=yEpts(2)-yEpts(1);
    aex=min(xEpts);
    aey=min(yEpts);
    Lx=max(xEpts)-min(xEpts)+hex;
    Ly=max(yEpts)-min(yEpts)+hey;
    down=sup/2-1;
    up=sup/2;
    oddsup = false;
    if (mod(sup,2) == 1)
        down = floor(sup/2);
        up = down;
        oddsup = true;
    end
    mvals=-down:up;
    SWeights=sparse(Nx*Ny*Nz,NIB);
    IWeights=sparse(Nx*Ny*Nz,NIB);
    for ilam=1:NIB
        zpts=1:Nz;
        zpts=zpts(abs(zEpts-IBpts(ilam,3)) <= sup/2*hex);
        % Fast Gaussian gridding in x and y
        floory=mod(floor((IBpts(ilam,2)-aey)/hey+1e-10),Ny)+1;
        % correct mvals for when particle is closer to right face of cell i
        if oddsup && abs(yEpts(floory) - IBpts(ilam,2)) > hey/2
            floory = floory + 1;
        end  
        yclose=yEpts(floory)+roundm(IBpts(ilam,2)/Ly)*Ly;
        % Compute the y weights
        E1y = exp(-(IBpts(ilam,2)-yclose)^2/(2*gw^2));
        E2y = exp((IBpts(ilam,2)-yclose)*Ly/(Ny*gw^2));
        ywts = E1y.*E2y.^mvals.*exp(-(mvals.*Ly/Ny).^2/(2*gw^2));
        ypts=floory-down:floory+up;
        ypts=mod(ypts,Ny);
        ypts(ypts==0)=Ny;
        % Compute the x weights
        floorx=mod(floor((IBpts(ilam,1)-aex)/hex+1e-10),Nx)+1;
        % correct mvals for when particle is closer to right face of cell i
        if oddsup && abs(xEpts(floorx) - IBpts(ilam,1)) > hex/2
            floorx = floorx + 1;
        end      
        xclose=xEpts(floorx)+roundm(IBpts(ilam,1)/Lx)*Lx;
        E1x = exp(-(IBpts(ilam,1)-xclose)^2/(2*gw^2));
        E2x = exp((IBpts(ilam,1)-xclose)*Lx/(Nx*gw^2));
        xwts = E1x.*E2x.^mvals.*exp(-(mvals.*Lx/Nx).^2/(2*gw^2));
        xpts=floorx-down:floorx+up;
        xpts=mod(xpts,Nx);
        xpts(xpts==0)=Nx;
        xywts=1/(2*pi*gw^2)*xwts'*ywts;
        xywts=xywts(:);
        xyinds=xpts'+Nx*(ypts-1);
        xyinds=xyinds(:);
        zwts=1/sqrt(2*pi*gw^2)*exp(-(zEpts(zpts)-IBpts(ilam,3)).^2/(2*gw^2));
        ziwts=zwts.*iwtsz(zpts)';
        totwts=xywts*zwts';
        totiwts=xywts*ziwts'*hex*hey;
        totinds=xyinds+Nx*Ny*(zpts-1);
        SWeights=SWeights+sparse(totinds(:),ilam,totwts(:),Nx*Ny*Nz,NIB);
        IWeights=IWeights+sparse(totinds(:),ilam,totiwts(:),Nx*Ny*Nz,NIB);
    end
end

function val = roundm(x)
    val = round(x);
    if (mod(x,1)==0.5)
        if (x < 0)
            val = val+1;
        end
    end
end
        
    
