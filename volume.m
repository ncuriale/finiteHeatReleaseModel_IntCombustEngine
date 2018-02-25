function V = volume(theta,Vd,r,R)

% This function provides the volume as a function of crank angle

V = Vd/(r-1)+ (Vd/2)*(R+1-cosd(theta)-((R^2)-(sind(theta))^2)^(0.5));