% µ÷ÕûÏß¿í
function y = myelipsnorm(m,cov,level,coline,lw)
%
%
% draws one contour plot of a bivariate
% gaussian density with mean "m" and covariance
% matrix "cov". 
% The level is controled by "level".
% line color and linestyle is controled by 'coline'
% linewidth is controled by 'lw'
% -----------------------------------------------------------------------
% ref Mario A. T. Figueiredo and Anil K. Jain

% ----------------------------------------------------------------------
if nargin<4 
   dashed=0;
end
[uu,ei,vv]=svd(cov);
a = sqrt(ei(1,1)*level*level);
b = sqrt(ei(2,2)*level*level);
theta = [0:0.01:2*pi];
xx = a*cos(theta);
yy = b*sin(theta);
cord = [xx' yy']';
cord = uu*cord;
plot(cord(1,:)+m(1),cord(2,:)+m(2),coline,'LineWidth',lw)
end