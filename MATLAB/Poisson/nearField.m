% Function to evaluate the near field for a given set of points and
% charges. 
% Inputs: points (nearpts), charges (nearcharges), number of charges in the
% system (nC - this does NOT include any images so it will in general be
% less than length(nearcharges)), width of Gaussian gw, Ewald parameter p,
% length of periodic domain L, and truncation of near field distance rcut. 
% Outpts: values of near field for phi and E at each of the nC charges. 
function [phiNear, ENear] = nearField(nearpts,nearcharges,nC,gw,p,L,rcut)
    % Figure out where to truncate
    phiNear=zeros(nC,1);
    ENear=zeros(nC,3);
    g0=Gnear(0,gw,p);
    for iC=1:nC
        phiNear(iC)=phiNear(iC)+nearcharges(iC)*g0;
        for jC=1:length(nearcharges)
            if (iC~=jC)
                pt2 = calcNearImg(nearpts(iC,:),nearpts(jC,:),L);
                r=norm(nearpts(iC,:)-pt2);
                if (r < rcut)
                    rhat = (nearpts(iC,:)-pt2)/r;
                    gr = Gnear(r,gw,p);
                    grE = GnearE(r,gw,p);
                    phiNear(iC)=phiNear(iC)+nearcharges(jC)*gr;
                    ENear(iC,:)=ENear(iC,:)+nearcharges(jC)*grE*rhat;
                end
            end
        end
    end
end

% Find the nearest periodic image 
function pt = calcNearImg(pt1,pt2,L)
    xar=[-L;0;L;-L;0;L;-L;0;L];
    yar=[-L;-L;-L;0; 0; 0; L; L; L];
    d=pt1(1:2)+[xar yar]-pt2(1:2);
    [~,index]=min(sum(d.*d,2));
    pt = pt2-[xar(index) yar(index) 0];
end