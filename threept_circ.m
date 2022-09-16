function [center, radius] = threept_circ( Ax, Ay, Bx, By, Cx, Cy )
%THREEPT_CIRC takes in 3 points (A, B, C) and computes the center and
%radius of the circumcircle defined by this triangle. It takes advantage of
%the fact that perpendicular bisectors of two sides of a triangle meet at
%the circumcenter, which is the center of the circle we are looking for.
%full explanation can be found here: http://paulbourke.net/geometry/circlesphere/
A = [Ax Ay]; B = [Bx By]; C = [Cx Cy];

if isempty(A) || isempty(B) || isempty(C)
    error('require 3 points!');
end
if ~isequal(size(A),[1 2]) || ~isequal(size(B),[1,2]) || ~isequal(size(C),[1,2])
    error('require 3 [1 2] dimension points!');
end

slopeA = (B(2)-A(2))/(B(1)-A(1)); %ma
slopeB = (C(2)-B(2))/(C(1)-B(1)); %mb

if slopeA == slopeB
    error('Err:COLIN','3 points MUST not be collinear to calculate a circumcircle!');
%     center = [NaN NaN];
%     radius = NaN;
%     warning('WARN:COLIN','3 points MUST not be collinear to calculate a circumcircle!');
end
sideA_mid_x = (A(1)+B(1))/2;
sideA_mid_y = (A(2)+B(2))/2;

% x = ma*mb(y1-y3)+mb(x1+x2)-ma(x2+x3) / 2(mb-ma)
center_x = (slopeA*slopeB*(A(2)-C(2))+slopeB*(A(1)+B(1))-slopeA*(B(1)+C(1)))/(2*(slopeB-slopeA));
% subsitute center_x into equation of perpendicular bisector for side A,
% which is ya = -1/ma( x - (x1+x2)/2 ) + (y1+y2)/2
center_y = (-1/slopeA)*(center_x - sideA_mid_x) + sideA_mid_y;
center = [center_x center_y];

radius = norm([abs(A(1)-center_x) abs(A(2)-center_y)]);

end

