function [x,y] = Ellipse(c,Q,theta)
%
% Compute radius 
r2 = 1/sqrt(Q(1,1));
r1 = 1/sqrt(Q(2,2));
xc = c(1);
yc = c(2);
% compute points corresponding to axis-oriented ellipse
t = linspace(0, 2*pi, 200);
xt = r1 * cos(t);
yt = r2 * sin(t);
% aply rotation by angle theta
cot = cos(theta); sit = sin(theta);
x = xt * cot - yt * sit + xc;
y = xt * sit + yt * cot + yc;
end

