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

b = 0.0755;
b_iter = 0.0001; %incremental bore size for iteration
power_iter=0;%initialize iteration value
n_cyl=4;%number of cylinders
eff=0.34;%efficiency of brake power
power_wanted=50000;%power wanted by entire engine
while power_iter<(power_wanted/eff/n_cyl) % iteration loop to optimize bore size
    
clear T %normalize T array for concatenation later

s = b;
l = 1.5*s;
r = 10;
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
Qin = 3000.0;
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

Pi = 100e3;
Pe = 140e3;
Ti = 300;
Te = 400;
Tw = 400;
f  = 0.001;
N = 3000;

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
fprintf('| Heat input [J]                   | %9.5e |\n',Qin)
fprintf('| Max. pressure [kPa]              | %9.5e |\n',max(P)/1000)
fprintf('| Max. temperature [K]             | %9.5e |\n',max(T))
fprintf('| Indicated work [J]               | %9.5e |\n',W)
fprintf('| Indicated power [W]              | %9.5e |\n',4*W*N/60/2)
fprintf('| Ideal Thermal Efficiency         | %9.5e |\n',NTeff)
fprintf('| IMEP                             | %9.5e |\n',IMEP)
fprintf('| Bore Size [m]                    | %9.5e |\n',b)
fprintf('| Power Iteration Loop Error [W]   | %9.5e |\n',power_wanted/eff/n_cyl-(W*N/60/2))
fprintf('+------------------------------------------+\n\n')

power_iter=(W*N/60/2);%re-initialize power for loop
b=b+b_iter;%increase bore size
end
