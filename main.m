% Spacecraft Attitude Dynamics and Control 
% A.Y. 2022-2023, prof. Franco Bernelli Zazzera
% Final project - 
% Cambielli Alessandro 

% Add folder & Data

%
%
%                                        ^
%                                        |
%                                        |   z
%
%                             +--------------------------+
%                            /                          /|
%                           /                          / |
%                          /                          /  |
%                         +--------------------------+   |  ---> y
%                         |                          |   +
%                         | b                        | h/
%                         |           a              | /
%                         +--------------------------+
%                                     /
%                                    /  x
%                                  \/
%

clear
clc

addpath('.\data')
addpath('.\scripts\')
addpath('functions\')

stepsize = 0.1; % set lower value for a more refined smulation

% DATA

%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 3.98600433e+5;  % Earth gravitational constant  [km^3/s^2]
c = 2.99792*10^8;  % Light velocity  [m/s] 

%%%%%%%%%%%%%%%%%%%%%%%%%% Orbit Initial Data %%%%%%%%%%%%%%%%%%%%%%%%%%

R_e = 6378;  % Earth's radius  [km]
a = 8152;  % Semi-major axis  [km]             
e = 0.1195;  % Eccentricity  [-]             
i = deg2rad(21.8583);  % Inclination  [rad]           
theta0 = deg2rad(30);  % true anomaly initial condition  [rad] 
n = sqrt(mu/a^3);  % mean angular velocity  [rad/s]         
T = 2*pi/n;  % orbit period  [s]   

%%%%%%%%%%%%%%%%%%%%%% Spacecraft Characteristics %%%%%%%%%%%%%%%%%%%%%%

MB = [1186; 1.8; 1.4; 1.5];% Main body  [kg; m; m; m] Mass, a, b, h
SP = [18.4; 1.78; 3.2; 0.1];% Solar panels  [kg; m; m; m] Mass, a, b, h

% Inertia Moments main body
Ix_mb = MB(1)/12 * (MB(3)^2 + MB(2)^2); %  [Kg*m^2]
Iy_mb = MB(1)/12 * (MB(3)^2 + MB(4)^2); %  [Kg*m^2]
Iz_mb = MB(1)/12 * (MB(2)^2 + MB(4)^2); %  [Kg*m^2]

% Inertia Moments solar panels
Ix_sp = (2*SP(1))/12 * (SP(3)^2 + SP(2)^2); %  [Kg*m^2]
Iy_sp = (2*SP(1))/12 * (SP(3)^2 + SP(4)^2); %  [Kg*m^2]
Iz_sp = (2*SP(1))/12 * (SP(2)^2 + SP(4)^2); %  [Kg*m^2]

% Terms for huygens steiner theorem
Ix_hs = Ix_sp + 2*SP(1)*(2.5)^2; %  [Kg*m^2]
Iy_hs = Iy_sp; % distance = 0  [Kg*m^2]
Iz_hs = Iz_sp + 2*SP(1)*(2.5)^2; %  [Kg*m^2]

% total inertia moments
Ix = Ix_mb + Ix_hs; %  [Kg*m^2]
Iy = Iy_mb + Iy_hs; %  [Kg*m^2]
Iz = Iz_mb + Iz_hs; %  [Kg*m^2]

Inertia_matrix = [Ix 0 0; 0 Iy 0; 0 0 Iz];  % inertia matrix  [Kg*m^2]
Inverse_inertia = inv(Inertia_matrix);  % inverse inertia matrix  [Kg*m^2]        

%%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%

wx0 = deg2rad(3);  %omega x  [rad/s]        
wy0 = deg2rad(11);  %omega y  [rad/s]                          
wz0 = deg2rad(14);  %omega z  [rad/s]                        
w0 = [wx0; wy0; wz0];  %omega vector  [rad/s] [3x1] 
A_BN0 = eye(3);

%%%%%%%%%%%%%%%%%%%%%%%%%% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravity Gradient
inertia1 = Iz-Iy;  % term to compute gravity gradient  [Kg*m^2]
inertia2 = Ix-Iz;  % term to compute gravity gradient  [Kg*m^2]
inertia3 = Iy-Ix;  % term to compute gravity gradient  [Kg*m^2]

% Magnetic Dsiturbance
m = [0.1; 0.1; 0.1];
[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
G = [gn, gm, gvali, gsvi];
H = [hn, hm, hvali, hsvi];

% Solar Radiation Pressure
rho_d_MB = 0.1;
rho_s_MB = 0.5;
rho_s_SP = 0.1;
rho_d_SP = 0.1;

%%%%%%%%%%%%%%%%%%%%%%% Sensors characteristics %%%%%%%%%%%%%%%%%%%%%%%

% Earth Horizon Sensor STD16
EHSampleRate = 1;  % Earth Horizon sensor sample rate  [Hz]  
EHAccuracy = 0.1020;  % Earth Horizon sensor accuracy HSNS  [deg]

% Gyroscope STIM300
GyroSampleRate = 10;  % Gyroscope sample rate  [Hz] 
ARW = 0.15;  % Angular Randon Walk gyroscope  [deg/sqrt(h)]
RRW = 0.3;  % Rate Randon Walk gyroscope  [deg/h]

% Magnetometer DTFM100S
MMAccuracy = 0.003;  % +-0.3%
MMSampleRate = 1;  % Magnetometer sample rate  [Hz] to speed up the simulation
MMnoise = 15e-9*[1; 1; 1];  % noise vector (bias)  [T]

% Extended State Observer
Lw = 0.8;
Ld = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%% Actuators characteristics %%%%%%%%%%%%%%%%%%%%%%%

% Magnetic Coils MT400-2-L model
maxdipole = 400; %  [Am^2]

% Reaction Wheels RW6000 astrofein
%wheel_axis = [1/sqrt(3); 1/sqrt(3); 1/sqrt(3)];
wheel_axis = [1; 0; 0];
Inertia_RW = 2.274e-1;  % [kg m^2]
h_RW_max = 100;  % [Nms]
torque_RW_max = 0.8;  % [Nm]
hr0_1RW = 0;
hr0 = [0;0;0;0];
A = 1/sqrt(3)*[-1, 1, 1, -1;...
    -1, -1, 1, 1;...
    1, 1, 1, 1];
A_star = sqrt(3)/4*[-1, -1, 1;...
    1, -1, 1;...
    1, 1, 1;...
    -1, 1, 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Detumbling %%%%%%%%%%%%%%%%%%%%%%%%%%%%

first_detumbling_time_1RW = 20000;  % first control law time with 1 RW
first_detumbling_time = 5000;  % first control law time with 4 RW
detumblingkb = 1e9;  
detumblingkp = 0.1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Slew %%%%%%%%%%%%%%%%%%%%%%%%%%%%

slewk1 = 5e-1; 
slewk2 = 5e-1; 

 
%% Uncontrolled Dynamics:

startTime = 0;
stopTime = T;

outUncontrolled = sim('uncontrolled');

% uncontrolled dynamics plots:

% Setup for the plots:
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

figure('Name','Uncontrolled - Angular Velocity'),
hold on, grid on, box on
plot(outUncontrolled.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([0, outUncontrolled.tout(end)])
ylim([-0.35, 0.35])
% saveFigAsPdf('uncont_omega',0.5,2)

figure('Name','Uncontrolled - Angular Velocity Detail'),
hold on, grid on, box on
plot(outUncontrolled.w.Time(1:(1000/stepsize)),squeeze(outUncontrolled.w.Data(:,:,1:(1000/stepsize))), 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([0, outUncontrolled.tout(1000/stepsize)])
ylim([-0.3, 0.3])
% saveFigAsPdf('uncont_omega',0.5,2)

figure('Name','Uncontrolled - Pointing Error Detail'),
hold on, grid on, box on
plot(outUncontrolled.point_error.Data(1:(1000/stepsize)), 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')
xlim([0, outUncontrolled.tout(1000/stepsize)])
ylim([0, 200])
% saveFigAsPdf('uncont_pointerror',0.5,2)

figure('Name','Uncontrolled Dynamics - Distubing Torques'),
hold on, grid on, box on
plot(outUncontrolled.M_GG.Time(1:end),squeeze(outUncontrolled.M_GG.Data), 'LineWidth', 1.5)
plot(outUncontrolled.M_SRP.Time(1:end),squeeze(vecnorm(outUncontrolled.M_SRP.Data(:,:,:))), 'LineWidth', 1.5)
plot(outUncontrolled.M_MAG.Time(1:end),squeeze(vecnorm(outUncontrolled.M_MAG.Data(:,:,:))), 'LineWidth', 1.5)
legend('GG','SRP','MAG')
xlabel('$t [s]$'), ylabel('Torque $[N m]$')
xlim([0, outUncontrolled.tout(end)])
% saveFigAsPdf('uncont_disturbances',1,2)


%% Detumbling 1 Reaction Wheel:

startTime = 0;
stopTime = 70000; 

outDetumbling_1RW = sim('detumbling1');

% Final Conditions:
finalDetumbling_A_BN_1RW = outDetumbling_1RW.A_BN(:,:,end); 
finalDetumbling_w_1RW = outDetumbling_1RW.w.Data(:,1,end);
finalDetumbling_theta_1RW = outDetumbling_1RW.theta(end);
finalDetumbling_hr_1RW = outDetumbling_1RW.hr(end)';

% plots detumbling with 1 reaction wheel

% Setup for the plots:
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

figure('Name','Detumbling 1RW - Angular Velocity'),
hold on, grid on, box on
plot(outDetumbling_1RW.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([0, outDetumbling_1RW.tout(end)])
% saveFigAsPdf('detum_1RW_omega',0.5,2)

figure('Name','Detumbling 1RW - Pointing Error'),
hold on, grid on, box on
plot(outDetumbling_1RW.point_error, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')
xlim([0, outDetumbling_1RW.tout(end)])

figure('Name','Detumbling 1RW - Magnetorquer Control Torque'),
hold on, grid on, box on
plot(outDetumbling_1RW.Mc_MT, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([0, outDetumbling_1RW.tout(20000/stepsize)])

figure('Name','Detumbling 1RW - Reaction Wheels Control Torque'),
hold on, grid on, box on
plot(outDetumbling_1RW.Mc_RW, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([outDetumbling_1RW.tout(20000/stepsize), outDetumbling_1RW.tout(end)])

figure('Name','Detumbling 1RW - total Control Torque'),
hold on, grid on, box on
plot(outDetumbling_1RW.Mc_tot, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([0, outDetumbling_1RW.tout(end)])
% saveFigAsPdf('detum_1RW_totaltorque',0.5,2)

figure('Name','Detumbling 1RW - delta omega'),
hold on, grid on, box on
plot(outDetumbling_1RW.delta_omega, 'linewidth',1.5);
xlabel('t $[s]$'), ylabel('$[rad/s]$')
legend('$(\omega_x - \bar{\omega}_x)$',...
    '$(\omega_y - \bar{\omega}_y)$',...
    '$(\omega_z - \bar{\omega}_z)$')
xlim([0, outDetumbling_1RW.tout(end)])
% saveFigAsPdf('detum_1RW_deltaomega',0.5,2)

figure('Name','Detumbling 1RW - delta Md'),
hold on, grid on, box on
plot(outDetumbling_1RW.delta_Md, 'linewidth',1.5);
legend('$(M_{d,x} - \bar{M}_{d,x})$',...
    '$(M_{d,y} - \bar{M}_{d,y})$',...
    '$(M_{d,z} - \bar{M}_{d,z})$')
xlabel('t $[s]$'), ylabel('$[N m]$')
xlim([0, outDetumbling_1RW.tout(end)])
ylim([-8e-4, 7e-4])
% saveFigAsPdf('detum_1RW_deltaMd',0.5,2)



%% Detumbling 4 Reaction Wheels:

startTime = 0;
stopTime = 70000; % 70000 con 4 reaction wheels

outDetumbling = sim('detumbling2');

% Final Conditions:
finalDetumbling_A_BN = outDetumbling.A_BN(:,:,end); 
finalDetumbling_w = outDetumbling.w.Data(:,1,end);
finalDetumbling_theta = outDetumbling.theta(end);
finalDetumbling_hr = outDetumbling.hr(end,:)';

% plots detumbling with 4 reaction wheels

% Setup for the plots:
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

figure('Name','Detumbling - Angular Velocity'),
hold on, grid on, box on
plot(outDetumbling.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([0, outDetumbling.tout(end)])
% saveFigAsPdf('detum_omega',0.5,2)

figure('Name','Detumbling - Pointing Error'),
hold on, grid on, box on
plot(outDetumbling.point_error, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')
xlim([0, outDetumbling.tout(end)])

figure('Name','Detumbling - Magnetorquer Control Dipole'),
hold on, grid on, box on
plot(outDetumbling.D, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$D [A m^2]$')
legend('$D_x$','$D_y$','$D_z$')
xlim([0, outDetumbling.tout(first_detumbling_time/stepsize)])
ylim([-600, 600])
% saveFigAsPdf('detum_MTdipole',0.5,2)

figure('Name','Detumbling - Reaction Wheels Control Torque'),
hold on, grid on, box on
plot(outDetumbling.Mc_RW, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([outDetumbling.tout(first_detumbling_time/stepsize), outDetumbling.tout(end)])
% saveFigAsPdf('detum_RWtorque',0.5,2)

figure('Name','Detumbling - total Control Torque'),
hold on, grid on, box on
plot(outDetumbling.Mc_tot, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([0, outDetumbling.tout(end)])
% saveFigAsPdf('detum_totaltorquw',0.5,2)

figure('Name','Detumbling - delta omega'),
hold on, grid on, box on
plot(outDetumbling.delta_omega, 'linewidth',1.5);
xlabel('t $[s]$'), ylabel('$[rad/s]$')
legend('$(\omega_x - \bar{\omega}_x)$',...
    '$(\omega_y - \bar{\omega}_y)$',...
    '$(\omega_z - \bar{\omega}_z)$')
xlim([0, outDetumbling.tout(end)])
% saveFigAsPdf('detum_deltaomega',0.5,2)

figure('Name','Detumbling - delta Md'),
hold on, grid on, box on
plot(outDetumbling.delta_Md, 'linewidth',1.5);
legend('$(M_{d,x} - \bar{M}_{d,x})$',...
    '$(M_{d,y} - \bar{M}_{d,y})$',...
    '$(M_{d,z} - \bar{M}_{d,z})$')
xlabel('t $[s]$'), ylabel('$[N m]$')
xlim([0, outDetumbling.tout(end)])
ylim([-8e-4, 7e-4])
% saveFigAsPdf('detum_deltaMd',0.5,2)



%% Tracking 4 Reaction Wheels:

A_BN0 = finalDetumbling_A_BN; 
w0 = finalDetumbling_w;
theta0 = finalDetumbling_theta;
hr0 = finalDetumbling_hr;
startTime = 70000;
stopTime = 110000;

outSlew = sim('tracking');

% Final Conditions:
finalSlew_A_BN = outSlew.A_BN(:,:,end); 
finalSlew_w = outSlew.w(:,1,end);
finalSlew_theta = outSlew.theta(end);

% plots tracking with 4 reaction wheels

% Setup for the plots:
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

figure('Name','Tracking - Angular Velocity'),
hold on, grid on, box on
plot(outSlew.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
% saveFigAsPdf('tracking_omega',0.5,2)

figure('Name','Tracking - Angular Velocity Detail'),
hold on, grid on, box on
plot(outSlew.w, 'linewidth',1.5);
plot(outSlew.tout, n*outSlew.tout./outSlew.tout,'k--', 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$','n')
xlim([outSlew.tout(30000/stepsize), outSlew.tout(end)])
% saveFigAsPdf('tracking_omegadetail',0.5,2)

figure('Name','Slew - Pointing Error'),
hold on, grid on, box on
plot(outSlew.point_error, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')
xlim([outDetumbling.tout(end), outSlew.tout(20000/stepsize)])
% saveFigAsPdf('tracking_slewpointing',0.5,2)

figure('Name','Earth Pointing - Pointing Error'),
hold on, grid on, box on
plot(outSlew.point_error, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')
xlim([outSlew.tout(20000/stepsize), outSlew.tout(end)])
% saveFigAsPdf('tracking_earthpointerror',0.5,2)

figure('Name','Tracking - total Control Torque'),
hold on, grid on, box on
plot(outSlew.Mc_tot, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$M_c [N  m]$')
legend('$Mc_x$','$Mc_y$','$Mc_z$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
% saveFigAsPdf('tracking_totaltorque',0.5,2)

figure('Name','Tracking - delta omega'),
hold on, grid on, box on
plot(outSlew.delta_omega, 'linewidth',1.5);
xlabel('t $[s]$'), ylabel('$[rad/s]$')
legend('$(\omega_x - \bar{\omega}_x)$',...
    '$(\omega_y - \bar{\omega}_y)$',...
    '$(\omega_z - \bar{\omega}_z)$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
% saveFigAsPdf('tracking_deltaomega',0.5,2)

figure('Name','Tracking - delta Md'),
hold on, grid on, box on
plot(outSlew.delta_Md, 'linewidth',1.5);
legend('$(M_{d,x} - \bar{M}_{d,x})$',...
    '$(M_{d,y} - \bar{M}_{d,y})$',...
    '$(M_{d,z} - \bar{M}_{d,z})$')
xlabel('t $[s]$'), ylabel('$[N m]$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
ylim([-8e-4, 7e-4])
% saveFigAsPdf('tracking_deltaMd',0.5,2)

figure('Name','Tracking - delta Mc'),
hold on, grid on, box on
plot(outSlew.delta_Mc, 'linewidth',1.5);
legend('$(M_{c,x} - \bar{M}_{c,x})$',...
    '$(M_{c,y} - \bar{M}_{c,y})$',...
    '$(M_{c,z} - \bar{M}_{c,z})$')
xlabel('t $[s]$'), ylabel('$[N m]$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
% saveFigAsPdf('tracking_deltaMc',0.5,2)

figure('Name','Tracking - A error'),
hold on, grid on, box on
plot(outSlew.A_error, 'linewidth',1.5);
xlabel('t $[s]$'), ylabel('$[deg]$')
xlim([outDetumbling.tout(end), outSlew.tout(end)])
ylim([-1.5, 1.5])
% saveFigAsPdf('tracking_Aerror',0.5,2)











