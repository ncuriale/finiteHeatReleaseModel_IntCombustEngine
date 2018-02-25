function dVdt = dVdtheta(theta,Vd,R)

% This function provides the derivative of volume with respect to crank
% angle

dVdt =(Vd/2*(sind(theta)*(1+cosd(theta)*((R^2)-(sind(theta))^2)^(-1/2))))*pi/180;