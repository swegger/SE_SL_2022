function [x,y] = ellipse(xradius,yradius,x0,y0,dradians)
%% ellipse
%
%   [x,y] = ellipse(xradius,yradius,x0,y0,dradians)
%      Gives x-y coordinates for an ellipse centered at x0,y0 with a radius
%      along the x-axis of xradius and along the y-axis of yradius.
%
%%

t=-pi:dradians:pi;
x=x0+xradius*cos(t);
y=y0+yradius*sin(t);
