function dQdt = dQdtheta(theta,n,a,Qin,thetas,thetad)

% This function provides the derivative of heat release with respect to
% crank angle

if (theta < thetas)
    dQdt =0;
else
    xb= 1 - exp(-a*((theta-thetas)/thetad)^n);
    dQdt =n*a*(Qin/thetad)*(1-xb)*((theta-thetas)/thetad)^(n-1) ;
end