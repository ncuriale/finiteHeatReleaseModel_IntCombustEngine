function dQwdt = dQwdtheta(theta,P,b,N,Tw,V,s,m,l,viscg,cond)

% This function provides the derivative of heat transfer to the wall with 
% respect to crank angle

Rg = 287;
a = s/2;
y = a + l - ((l^2 - a^2*(sind(theta))^2)^0.5 + a*cosd(theta));
Aw = pi*b*y + pi/2*b^2;
Up = 2*N*s/60;
rho = m/V;
Re = rho*Up*b/viscg;
h = 0.49*(Re^0.7)*(cond/b);
Tg = P/(rho*Rg);

dQwdt = (h*Aw*(Tg-Tw)/(2*pi*N/60))*(pi/180);
%dQwdt=0;