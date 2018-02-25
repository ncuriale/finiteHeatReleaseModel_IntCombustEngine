%------------------------------------------------------------------%
% main.m - This is the main code for the finite heat release model %
% coupled with intake/exhaust stroke models                        %
%------------------------------------------------------------------%

clear all;

%------------------------------------------%
% Define physical parameters of the system %
%------------------------------------------%

% b .... bore [m]
% s .... stroke [m]
% l .... connecting rod length [m]
% r .... compression ratio
% Vd ... displaced volume (calculated) [m^3]
% R .... dimensionless group 2*l/s (calculated)

b = 0.089;
s = 0.063;
l = 0.11;
r = 8.5;
Vd = (pi*s*b^2)/4;
R = 2*l/s;

%---------------------------------------%
% Define finite heat release parameters %
%---------------------------------------%

% n ........ Weibe form factor
% a ........ Weibe efficiency factor
% Qin ...... total heat addition [J]
% thetad ... duration of heat release [degrees]
% thetas ... start of heat release [degrees]

n = 3;
a = 5;
As=15.03;
phi=1/.927;
Qc=25559000;
thetas = -10;
thetad = 40;

%-----------------------%
% Define gas properties %
%-----------------------%

% k ........... specific heat ratio (cp/cv)
% viscg ....... gas viscosity [Pa*s]
% cond ........ gas conductivity [W/mK]

k = 1.4;
viscg = 20e-6;
cond = 6e-2;

%-----------------------------%
% Define operating conditions %
%-----------------------------%

% Pi .......... inlet pressure [Pa]
% Pe .......... exhaust pressure [Pa]
% Ti .......... inlet temperature [K]
% Te .......... exhaust temperature (initial guess) [K]
% Tw .......... average cylinder wall temperature [K]
% f ........... residual fraction (initial guess)
% N ........... engine speed [RPM]

Pi = 101.3e3;
Pe = 140e3;
Ti = 299;
Te = 1029;
Tw = 400;
f  = 0.001;
N = 2816;

%-----------------------------%
% Define remaining parameters %
%-----------------------------%

% theta ....... array of angles [deg]
% ntheta ...... number of values for theta
% Rg .......... ideal gas constant [J/(kg K)]

theta = -180:1:180;
ntheta = length(theta);
Rg = 287;

%-----------------------------%
% Iterate to solve the system %
%-----------------------------%

res = 1e10;

iter = 1;
fprintf('Iterating ...\n')

while (res>1e-6),
    
   
    
    % Store Te and f from previous iteration
    
    Te0 = Te;
    f0 = f;

    % Temperature at end of intake

    T1 = (1-f)*Ti+f*Te*(1-(1-(Pi/Pe))*((k-1)/k));
    
    % Mass contained in cylinder after intake

    m = Pi*volume(180,Vd,r,R)/(Rg*T1); %ideal gas law
    Qin = 624.5/2; %(m*(1-f)/(As/phi))*Qc;
    %------------------------------------------------%
    % Solve the finite heat release ODE for pressure %
    %------------------------------------------------%

    % fun .745[........ function handle for dPdtheta that is a fn of theta and P only

    fun = @(thetavar,Pvar) dPdtheta(thetavar,Pvar,k,Vd,R,n,a,Qin,thetas,thetad,b,N,Tw,s,m,r,l,viscg,cond);
    [theta,P] = ode45(fun,theta,Pi);
    
    % Temperature during compression, heat addition, and expansion

    for i=1:ntheta
        T(i)= P(i)*volume(theta(i),Vd,r,R)/(m*Rg);
    end

    % Temperature and pressure after expansion

    T4 = T(ntheta);
    P4 = P(ntheta);

    % Temperature after blowdown

    T5 =T4*(Pe/P4)^((k-1)/k);

    % Temperature after exhaust stroke

    Te = T5 ;
    f =(1/r)*(Pe/P4)^(1/k);
    
    % Print convergence information
    
    res = max([abs(Te-Te0) abs(f-f0)]);
    fprintf(' Iteration %i ...\n',iter)
    fprintf('  Te residual: %9.5e\n',abs(Te-Te0));
    fprintf('  f residual:  %9.5e\n',abs(f-f0));
    
    iter = iter + 1;

end

fprintf('Converged.\n\n')
fprintf('Exhaust temperature: %9.5e\n',Te);
fprintf('Residual fraction:   %9.5e\n\n',f);

%----------------------------------------------%
% Add the intake and exhaust strokes to arrays %
%----------------------------------------------%

theta = -360:1:360;
nIntake = (length(theta)-ntheta)/2;
ntheta = length(theta);
nExhaust = nIntake;

for i=1:ntheta
    V(i)=volume(theta(i),Vd,r,R);
end

P=[Pi*ones(nIntake,1); P; Pe*ones(nExhaust,1)];
T=[T1*ones(nIntake,1); T'; Te*ones(nExhaust,1)];

%-----------------------%
% Compute the work done %
%-----------------------%

W = trapz(V,P);       %i.e. Z = trapz(Y), if Y is a vector, trapz(Y) is the 
                    %                   integral of Y.
                    %                   if Y is a matrix, trapz(Y) is a row
                    %                   vector with the integral over each
                    %                   column.
                    %i.e. Z = trapz(X,Y), computes the integral of Y with
                    %                     respect to X using trapezoidal 
                    %                     integration.

                   
%----------------------------------%
% Compute ideal thermal efficiency %
%----------------------------------%

NTeff = W/Qin;

%-----------------------%
%      Compute IMEP     %
%-----------------------%

IMEP = W/Vd ;

Cp=1450;
Eta_m=0.9;
EE=Cp*m*(1-f)*(Te-Ti)*N/60/2;
QL=Qin*N/60/2-EE-W*Eta_m*N/60/2;

%------%
% Plot %
%------%

figure(1)
plot(theta,P/1000,'b')
hold on
xlim([-360 360])
xlabel('Crank angle [deg]')
ylabel('Pressure [kPa]')

figure(2)
plot(theta,T,'b')
hold on
xlim([-360 360])
xlabel('Crank angle [deg]')
ylabel('Temperature [K]')

figure(3)
plot(V,P/1000, 'b')
hold on
plot([V(1) V(ntheta)],[P(1)/1000 P(ntheta)/1000])
hold on
xlabel('Volume [m^3]')
ylabel('Pressure [kPa]')

%---------------%
% Print summary %
%---------------%

fprintf('+------------------------------------------+\n')
fprintf('| Summary                                  |\n')
fprintf('+------------------------------------------+\n')
fprintf('| Heat input [J]                    | %9.5e |\n',Qin)
fprintf('| Max. pressure [kPa]               | %9.5e |\n',max(P)/1000)
fprintf('| Max. temperature [K]              | %9.5e |\n',max(T))
fprintf('| Indicated work [J]                | %9.5e |\n',W)
fprintf('| Indicated power [W]               | %9.5e |\n',W*N/60/2)
fprintf('| Indicated power [W] (4-cylinders) | %9.5e |\n',4*W*N/60/2)
fprintf('| Efficiency                        | %9.5e |\n',(W*N/60/2)/(0.0003044*Qc))
fprintf('| Ideal Thermal Efficiency          | %9.5e |\n',NTeff)
fprintf('| bsfc                              | %9.5e |\n',0.0003044/(W*N/60/2)*3600*1000000)
fprintf('| IMEP                              | %9.5e |\n',IMEP)
fprintf('| Exhaust Energy                    | %9.5e |\n',EE)
fprintf('| Energy Lost                       | %9.5e |\n',QL)
fprintf('+------------------------------------------+\n\n')