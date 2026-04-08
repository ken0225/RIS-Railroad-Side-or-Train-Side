% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/03/01), v0.7(2022/01/02)


%% Clean All & Timer
close all;
clear;
tic;

%% System Parameters Initialization

global c0 fc speed M N; % Global parameters

% -------- Changeable parameters begin -------- 
% You CAN change the parameters below to obtain different figures
% In this simulation we IGNORE the phase drifts


M = 64; N = M; % IRS has M x N elements



%  -------- Changeable parameters end -------- 

c0 = 299792458; % Light speed

fc= 2.4e9; % Carrier frequency is 2.4 GHz

disp(['Step 1: Parameter Initialization begin. In this setup, the IRS has ' num2str(M) ' x ' num2str(N) ' elements.']);

dx = (c0/fc)/(2*sqrt(pi)); dy = dx; % dxdy=\lambda^2/(4\pi), this means the element gain is 1. 

h_IRS = 3; % Height of the IRS

h_User = 1.5; % Height of the user

h_BS = 20; % Height of the BS

%rotate_angle = -pi/4;

speed = 55; % Vehicle speed, note that 300km/h = 83.3333m/s

total_time = 200; % Total time for one vehicle pass

Pt = 0.1; % Transmit power is 1W = 30dBm

N0 = 10e-12; % N0 = -80dBm

% The coordinate of the vehicle, i.e., the receiver
x_Vehicle_start = -5000;
y_Vehicle = 0+h_User;
z_Vehicle = 50; 

% The coordinate of the Base Station, i.e., the transmitter
x_BS = 0;
y_BS = 0+h_BS;
z_BS = 20; 

%% Compute the trajectory for one vehicle pass

disp(['Step 2: Compute the trajectory for one vehicle pass. The speed is ' num2str(speed) ' m/s.']);

p_BS = [x_BS, y_BS, z_BS];

p_vehicle_start = [x_Vehicle_start, y_Vehicle, z_Vehicle];

% Total trajectory = starting point + moving trajectory per second (we start from 1s rather than 0s)
p_vehicle_trajectory = [];

% The vehicle starts at p_vehicle_start, and travels at the speed during
% total_time period
for t = 1 : total_time
    
    temp_p_vehicle_trajectory =  p_vehicle_start+function_vehicle_moving_xdir(speed, t);
    
    p_vehicle_trajectory = [p_vehicle_trajectory; temp_p_vehicle_trajectory];
    
end

% 'A_0' is the gain for the direct path. It is a theoretical result.
[A_0] = function_A0(p_BS, p_vehicle_trajectory);

%% Compute the gain when the IRS is non-isotropic, after phase optimization, and w/o HWI
disp('Step 3: Compute the gain when the IRS is isotropic, after phase optimization, and w/o HWI.');

% Locate the positions for M and N
[centers_IRS] = function_centers_IRS(M, N, dx, dy);

% 'h_IRS' is the height of IRS
centers_IRS(:, 2) = centers_IRS(:, 2) + h_IRS;

% 'A_mn_theoretical' is the gain for the IRS path, and it is a theoretical result.
[A_mn, A_mn_matrix] = function_Amn_Bside(p_BS, centers_IRS, p_vehicle_trajectory, total_time);

% Compute the delays
[tau_0, tau_mn] = function_time_delay_Bside(centers_IRS, p_BS, p_vehicle_trajectory); 

% Optimal phase shift
phi_optimal = 2*pi*(fc*(tau_0-tau_mn) + ceil(-(fc * (tau_0-tau_mn))));

% Final phase part
exp_phi_optimal = exp(-1i*2*pi*fc.*tau_mn-1i.*phi_optimal);

% Total gain of theoretical result
A_mn_plus_A_0 = A_mn.'+A_0;

Td_1RIS_Bside = max(tau_mn+(phi_optimal./(2*pi*fc))) - tau_0;

%% Compute the  Gain/SE Approximation, i.e., Theorem 1
disp('Step 5: Compute the gain/SE Approximation, i.e., Theorem 1');

A_star = [];
Q = [];
a_HWI = 0; kappa_t = 0; kappa_r = 0; 

% Compute the A*
for t = 1 : total_time
    
    temp_A_star=function_Astar(A_mn_matrix(t,:));
    A_star = [A_star; temp_A_star];
    
end

% Compute the Q
for t = 1 : total_time
    
    temp_Q = function_Q(A_0(t), A_star(t), A_mn_matrix(t,:), a_HWI);
    Q(t) = temp_Q;
    
end

Q = Pt*Q;

% Compute the SNR approximation
SNR_approximation = (Q ./ ((kappa_t+kappa_r)*Q + N0))'; 

% Compute the SE approximation
SE_approximation = log2(1+SNR_approximation);


%% Plot the simulation results
close all;

disp('Step 7: Plot the simulation results.');

inst_received_power_no_IRS = Pt * (abs(A_0)) .^ 2;
inst_SNR_no_IRS = inst_received_power_no_IRS ./ (N0);
SE_no_IRS = log2(1 + inst_SNR_no_IRS);

inst_received_power_IRS_theoretical = Pt * (abs(A_mn_plus_A_0)) .^ 2;
inst_SNR_IRS_theoretical = inst_received_power_IRS_theoretical ./ (N0);
SE_IRS_no_HWI = log2(1 + inst_SNR_IRS_theoretical);



% Figure 1
figure(1); hold on; box on; grid on;

p1 = plot(1 : total_time, SE_no_IRS, 'r-', 'LineWidth', 1.25);
p2 = plot(1 : total_time, SE_IRS_no_HWI, 'g-', 'LineWidth', 1.25);
p5 = plot(1 : total_time, SE_approximation, 'b+', 'LineWidth', 1.25);

xlabel('Vehicle Moving Time (s)','Interpreter','LaTex');
ylabel('Spectral Efficiency (bit/s/Hz)','Interpreter','LaTex');
title('SE with RIS and Transceiver HWIs for One Vehicle Pass','Interpreter','LaTex');

toc;

save 'SE_IRS_no_HWI_Bside' SE_IRS_no_HWI
save 'Td_1RIS_Bside' Td_1RIS_Bside
%save 'EE_IRS_no_HWI_Bside' EE_IRS_no_HWI_Bside
