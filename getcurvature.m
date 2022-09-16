function curvature = getcurvature( curr_x, curr_y )
%GETCURVATURE.m This function takes in matched x and y coordinates and
%spits out a curvature value for all points along the line (curvature =
%1/radius of the osculating circle along the body length);
Lines=[(1:(size(curr_x,1)-1))' (2:size(curr_y,1))'];
% Get left and right neighbor of each points
Na = zeros(size(curr_x)); Nb = zeros(size(curr_y));
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);

% Check for end of line points, without a left or right neighbor
checkNa = Na == 0; checkNb = Nb == 0;
Naa = Na; Nbb = Nb;
Naa(checkNa) = find(checkNa); Nbb(checkNb) = find(checkNb);

% If no left neighbor use two right neighbors, and the same for right...
Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)= Naa(Naa(checkNb));


slopeA = (curr_y-curr_y(Na))./(curr_x-curr_x(Na)); %ma
slopeB = (curr_y(Nb)-curr_y)./(curr_x(Nb)-curr_x); %mb
colin_inds = find(slopeA==slopeB); %to exclude from analysis
%points = [curr_x(Na) curr_y(Na) curr_x curr_y curr_x(Nb) curr_y(Nb)];
% points = points(non_colin_inds,:);
Ax = curr_x(Na); Ay = curr_y(Na); Bx = curr_x; By = curr_y; Cx = curr_x(Nb); Cy = curr_y(Nb);
% Ax = points(:,1); Ay = points(:,2); Bx = points(:,3); By = points(:,4); Cx = points(:,5); Cy = points(:,6);
% tic
% points_table = table(Ax,Ay,Bx,By,Cx,Cy);
% circles = rowfun(@threept_circ,points_table,'OutputVariableNames',{'centers' 'radii'});
% toc
% radii = (circles.radii)';
%alternate method simply using a for loop -- this is much faster!!!
% tic
curvature = zeros(size(Ax));
for k = 1:size(Ax,1)
    if sum(k==colin_inds)==1 % if 3 neighboring points are colinear, the curvature is 0
        curvature(k) = 0;
    else
        [~, radius] = threept_circ( Ax(k), Ay(k), Bx(k), By(k), Cx(k), Cy(k) );
        curvature(k) = 1/radius;
    end
end
% toc
%curvature = (1./radii)';

end

