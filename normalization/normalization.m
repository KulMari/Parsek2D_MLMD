% normalization 

%% program to calculate plasma parameters for simulation with PARSEK
clear all;close all

%% Quantities defined by user

format long
% Mass-ratio
miome=100;
% Reference magnetic field - gauss ( 1 nT = 10^-5 gauss)
B0=0.001015011757700;
%B0=0.0043;
% Reference density - cm^-3
n0=10;
% Reference temperature - kelvin
T=2.968950532338670e+06;

%% Physical constants - CGS system

%mi=1.6726e-24;me=mi/miome; % choose which one to keep at physical value
me=9.1094e-28;mi=me*miome;  % between mi and me
kb=1.3807e-16;
q=4.8032e-10;
c=3e10;

%%  Quantities calculated

T=T/1.1604e4; % T in eV;

Va=B0/sqrt(4*pi*n0*mi); % Alfven velocity

ld=7.43e2*sqrt(T/n0);  % Debye length - cm

om_pi=sqrt(4*pi*n0*q^2/mi); % ion plasma frequency - rad/s
om_pe=om_pi*sqrt(miome); % electron plasma frequency - rad/s

om_ci=q*B0/mi/c; % ion cyclotron frequency - rad/s
om_ce=om_ci*miome; % electron cyclotron frequency - rad/s
%vthi=9.79e5*sqrt(2*T); vthe=vthi*sqrt(miome);    %  Thermal velocities (physical mi)
vthe=4.19e7*sqrt(2*T);  vthi=vthe/sqrt(miome);    %  Thermal velocities (physical me)


%% Box size and Nx 
% Box size is N times the wavelength of a given perturbation
% Assign here N and the wavevector k;
% Assign also the number of the cell Ncell;
Nk=0.85/cosd(80);
k=Nk*om_pi/c;   % reference wavevector k
Nw=1;           % The box will be Nw wavelength long
L=2*pi*Nw/k;    % Length of the box         

Ncell=320;


%% Time and number of cycles;
% The total time of simulation is given by N/growth-rate of a given
% perturbation

% Assign here N and the growth-rate gamma 
N=3;
gamma=0.71*om_ci;
Time=2*pi*N/gamma;
%Time=N/gamma;
% Assign here how many time-steps for an electron gyromotion;
steps=13 ;
dt=1/steps/om_ce;

Ncycles=Time/dt;


% Assign :
Npcel= 10^2; % Number of particles per cell;
Nparticles=Npcel*Ncell^2;

dx=L/Ncell;
%% Check some constraints
if (vthe*dt/dx>1)
disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp(' vthe < dx\dt  is not satisfied --> increase steps or Nw , or decrease Ncell or Nk')
return
end

if (vthe*dt/dx<0.1)
disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp(' vthe *dt/dx > 0.1  is not satisfied --> finite grid instability !!')
%return
end

%% Decide normalization quantities
Norm_time = '1/om_pi'; % 
Norm_velocity = 'c';
Norm_space = ['((' Norm_time ')*' Norm_velocity ')'];

norm_time=eval(Norm_time);
norm_velocity=eval(Norm_velocity);
norm_space=norm_time*norm_velocity;
norm_qom=q/mi;  
norm_B0=1/norm_space*norm_velocity^2/norm_qom; 
norm_rho=norm_B0/norm_space;


DT=dt/norm_time;
DX=dx/norm_space;
L_normalized=L/norm_space;
Time_normalized=Time/norm_time;
qom_i_normalized=q/mi/norm_qom;
qom_e_normalized=-q/me/norm_qom;
B0_normalized=B0/norm_B0;
rho=q*n0;
rho_normalized=rho/norm_rho;
vthi_normalized=vthi/norm_velocity;
vthe_normalized=vthe/norm_velocity;
c_normalized=c/norm_velocity;
Va_normalized=Va/norm_velocity;
om_pi_normalized=om_pi*norm_time;
om_pe_normalized=om_pe*norm_time;
om_ci_normalized=om_ci*norm_time;
om_ce_normalized=om_ce*norm_time;
%
k_normalized=k*norm_space;
delta_E=om_ce_normalized^1.5/(vthe_normalized^0.5);
delta_rho=om_ce_normalized^2.5/(vthe_normalized^1.5);
delta_j=om_ce_normalized^2.5/(vthe_normalized^0.5);

%% Display
disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp (['Physical size of the box is ' num2str(L/1e5) ' km']);
disp (['Debye length is ' num2str(ld/1e2) ' m']);
disp (['Number of cell (dx = ' num2str(dx/ld) ' ld) is ' num2str(Ncell) ])
disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp (['Total physical time of the simulation is ' num2str(Time) ' sec']);
disp (['Physical dt is ' num2str(dt) ' sec']);
disp (['Number of cycles is ' num2str(Ncycles) ]);

disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp('Check some constrains on dt and dx :' )
disp(' Constrain          True ?')
disp(['vthe < dx\dt          ' num2str(vthe<dx/dt) ])
disp(['(om_pe*dt)^2 < Npc    ' num2str((om_pe*dt)^2 < Nparticles/Ncell^2) ])

disp(' ');disp('%%%%%%%%%%%%%%%%%%%%%')
disp(' Normalization parameters based on:')

disp(['Velocity--> ' Norm_velocity ' ; Time --> ' Norm_time])
disp(' ')
disp(['qom_i = ' num2str(qom_i_normalized) ' ; qom_e = ' num2str(qom_e_normalized) ])
disp(['L = ' num2str(L_normalized) ';    dx = ' num2str(DX) '; ncell = ' num2str(Ncell)])
disp(['Time = ' num2str(Time_normalized) '; dt = ' num2str(DT) '; ncycles = ' num2str(Ncycles)])
disp(['B0 = ' num2str(B0_normalized)]);
disp(['rho = ' num2str(rho_normalized)]);
disp(['c = ' num2str(c_normalized)]);
disp(['Va = ' num2str(Va_normalized)]);
disp(['vthi = ' num2str(vthi_normalized)]);
disp(['vthe = ' num2str(vthe_normalized)]);
disp(['om_pi = ' num2str(om_pi_normalized)]);
disp(['om_pe = ' num2str(om_pe_normalized)]);
disp(['om_ci = ' num2str(om_ci_normalized)]);
disp(['om_ce = ' num2str(om_ce_normalized)]);
disp(['delta_E = ' num2str(delta_E,'%10.12f')]);
disp(['delta_rho = ' num2str(delta_rho,'%10.12f')]);
disp(['delta_j = ' num2str(delta_j,'%10.12f')]);


disp(['1 electron gyroradius is ' num2str(0.5*vthe_normalized/om_ce_normalized/DX) ' cell size']);
disp(['vthe*dt/dx = ' num2str(vthe*dt/dx)]);