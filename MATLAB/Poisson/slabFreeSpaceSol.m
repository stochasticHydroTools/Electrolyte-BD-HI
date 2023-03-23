% Free space solution inside the slab. Uses infinitely many images (100
% total images for each charge - this might not be enough depending on how 
% fast the series decays)
function [phiFree, Efree] = slabFreeSpaceSol(pts,charges,gw,H,epsratiot,epsratiob) 
    nC = length(charges);
    bimgpts = [pts(:,1:2) -pts(:,3)];
    bimgch = -charges*(1/epsratiob-1)/(1/epsratiob+1);
    timgpts = [pts(:,1:2) 2*H-pts(:,3)];
    timgch = -charges*(1/epsratiot-1)/(1/epsratiot+1);
    pwithIm = [pts; bimgpts; timgpts];
    cwithIm = [charges; bimgch; timgch];
    for iImg=1:49
        % Images of images
        temppts=bimgpts;
        tempch=bimgch;
        bimgpts = [timgpts(:,1:2) -timgpts(:,3)];
        bimgch = -timgch*(1/epsratiob-1)/(1/epsratiob+1);
        timgpts = [temppts(:,1:2) 2*H-temppts(:,3)];
        timgch = -tempch*(1/epsratiot-1)/(1/epsratiot+1);
        pwithIm = [pwithIm; bimgpts; timgpts];
        cwithIm = [cwithIm; bimgch; timgch];
    end
    % Compare to the free space answer
    pwithIm=[pwithIm; 0 0 0];
    cwithIm=[cwithIm; 0];
    phiFree=zeros(length(cwithIm),1);
    Efree=zeros(length(cwithIm),3);
    for iC=length(cwithIm):-1:1
        phiFree(iC)=phiFree(iC)-cwithIm(iC)*(1/(4*gw*pi^(3/2))); % self potential
        for jC=iC-1:-1:1
            r=norm(pwithIm(iC,:)-pwithIm(jC,:));
            rhat = (pwithIm(iC,:)-pwithIm(jC,:))/r;
            gr = erf(r/(2*gw))./(4*pi*r);
            phiFree(iC)=phiFree(iC)-cwithIm(jC)*gr;
            phiFree(jC)=phiFree(jC)-cwithIm(iC)*gr;
            Efree(iC,:)=Efree(iC,:)+cwithIm(jC)*rhat.*(erf(r/(2*gw))./(4*pi*r.^2)-...
            exp(-r.^2/(4*gw^2))./(4*pi^(3/2)*gw*r));
            Efree(jC,:)=Efree(jC,:)+cwithIm(iC)*(-rhat).*(erf(r/(2*gw))./(4*pi*r.^2)-...
            exp(-r.^2/(4*gw^2))./(4*pi^(3/2)*gw*r));
        end
    end
    % The potential on z < 0 is the same as if it was a wall at z = 0 with all
    % of the charges and image charges above z = 0. For our algorithm this
    % means we will consider the charges and the top charges (and maybe their
    % first images)
    bpt=pwithIm(end,:);
    phibot=0;
    Ebot=[0 0 0];
    pbot=[pts; pwithIm(pwithIm(:,3) > H,:)];
    cbot=[charges(1:nC); cwithIm(pwithIm(:,3) > H)]*2/(1/epsratiob+1);
    for iC=1:length(cbot)
        r=norm(bpt-pbot(iC,:));
        rhat=(bpt-pbot(iC,:))/r;
        phibot=phibot-cbot(iC)*erf(r/(2*gw))./(4*pi*r);
        Ebot=Ebot+cbot(iC)*rhat*(erf(r/(2*gw))./(4*pi*r.^2)-...
            exp(-r.^2/(4*gw^2))./(4*pi^(3/2)*gw*r));
    end
    abs(phibot-phiFree(end))
    abs([Ebot(1:2) Ebot(3)/epsratiob]-Efree(end,:))
    phiFree=-(phiFree(1:nC)-phiFree(end));
    Efree=Efree(1:nC,:);
end