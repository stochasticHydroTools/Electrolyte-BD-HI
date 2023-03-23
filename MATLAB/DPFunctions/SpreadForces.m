% This function spreads forces on the grid
% while taking into account the images
 
% Input: 
%   S        - spread matrix
%   forces   - forces on particles (Np x 3)
%   imS      - image spread matrix
%   imForces - image forces
%             - these should be EQUAL to the forces on the corresponding particle
%             - we end up subtracting the spread of the image forces from that
%               of the forces on the actual particles
%   nxy,Nz   - number of points in x/y and z

% Output:
%   gridf,gridg,gridh - x,y and z component of spread forces on grid

function [gridf,gridg,gridh] = SpreadForces(S,forces,imS,imForces,nxy,Nz)
    % spread forces on the particles
    [gridf,gridg,gridh] = spread(S,forces,nxy,Nz);
    [nIm,~] = size(imForces);
    % add image correction to spread forces 
    if nIm ~= 0
        [imGridf, imGridg, imGridh] = spread(imS,imForces,nxy,Nz);
        gridf = gridf-imGridf; gridg = gridg-imGridg; gridh = gridh-imGridh;
    end
end
