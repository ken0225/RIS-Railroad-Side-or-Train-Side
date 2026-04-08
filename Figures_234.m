% Author: Ke(Ken) WANG, Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Version: v0.1 (Jan 07, 2022)

%% Clean & Timer
close all;
clear;
tic;

%% System Parameters

global c0 fc speed M N total_time;

% Adjustable parameters
M = 64; N = M; % IRS size (M x N)

c0 = 299792458; % Speed of light
fc = 2.4e9;     % Carrier frequency

disp(['Step 1: Initialization. IRS size = ' num2str(M) ' x ' num2str(N)]);

dx = (c0/fc)/(2*sqrt(pi));
dy = dx; % Element size (unit gain)

h_IRS = 3;
h_User = 1.5;
h_BS = 20;

rotate_angle = pi/4;
vehicle_length = 3;

speed = 55;        % m/s
total_time = 200;  % s

Pt = 0.1;          % Transmit power
N0 = 10e-12;       % Noise power

% User (receiver)
x_User_start = -5000;
y_User = 0 + h_User;
z_User = 50;

% Base station (transmitter)
x_BS = 0;
y_BS = 0 + h_BS;
z_BS = 20;

%% Vehicle Trajectory

disp(['Step 2: Trajectory computation. Speed = ' num2str(speed) ' m/s']);

p_BS = [x_BS, y_BS, z_BS];
p_User_start = [x_User_start, y_User, z_User];

p_User_trajectory = [];

for t = 1:total_time
    temp = p_User_start + function_vehicle_moving_xdir(speed, t);
    p_User_trajectory = [p_User_trajectory; temp];
end

[A_0] = function_A0(p_BS, p_User_trajectory);

%% IRS Gain (no HWI, phase optimized)

disp('Step 3: IRS gain (vehicle-side, ideal)');

[centers_IRS] = function_centers_IRS_zy(M, N, dx, dy);
centers_IRS(:,2) = centers_IRS(:,2) + h_IRS;
centers_IRS = centers_IRS * function_rotate_IRS(rotate_angle);

p_IRS_start = [x_User_start + vehicle_length, y_User, z_User];

p_IRS_trajectory = zeros(M*N, 3, total_time);

for t = 1:total_time
    temp = p_IRS_start + centers_IRS + function_vehicle_moving_xdir(speed, t);
    p_IRS_trajectory(:,:,t) = temp;
end

[A_mn, A_mn_matrix] = function_Amn_Vside(p_BS, p_User_trajectory, p_IRS_trajectory, total_time);

[tau_0, tau_mn] = function_time_delay_Vside(p_IRS_trajectory, p_BS, p_User_trajectory);

k_imin = ceil(-(fc*(tau_0 - tau_mn)));
phi_optimal = 2*pi*(fc*(tau_0 - tau_mn) + k_imin);

exp_phi_optimal = exp(-1i*2*pi*fc.*tau_mn - 1i.*phi_optimal);

A_mn_plus_A_0 = A_mn + A_0;

%% Delay Spread

Td_1RIS_Vside = max(tau_mn + (phi_optimal./(2*pi*fc))) - tau_0;

%% SE Approximation (Theorem 1)

disp('Step 5: SE approximation');

a_HWI = 0;
kappa_t = 0;
kappa_r = 0;
A_star = [];
Q = [];

for t = 1:total_time
    A_star = [A_star; function_Astar(A_mn_matrix(t,:))];
end

for t = 1:total_time
    Q(t) = function_Q(A_0(t), A_star(t), A_mn_matrix(t,:), a_HWI);
end

Q = Pt * Q;

SNR_approximation = (Q ./ ((kappa_t + kappa_r)*Q + N0))';
SE_approximation = log2(1 + SNR_approximation);

%% Plot

close all;
disp('Step 7: Plot results');

moving_distance = linspace(1, total_time, total_time) .* speed;

inst_received_power_no_IRS = Pt * (abs(A_0)).^2;
SE_no_IRS = log2(1 + inst_received_power_no_IRS./N0);

inst_received_power_IRS = Pt * (abs(A_mn_plus_A_0)).^2;
SE_IRS_no_HWI_Vside = log2(1 + inst_received_power_IRS./N0);

% Figure 1: SE vs distance
figure(1); hold on; box on; grid on;

a = load('SE_IRS_no_HWI_Bside'); SE_IRS_no_HWI_Bside = a.SE_IRS_no_HWI;

plot(moving_distance, SE_no_IRS,           'k:',  'LineWidth', 3);
plot(moving_distance, SE_IRS_no_HWI_Vside, 'b-',  'LineWidth', 3);
plot(moving_distance, SE_IRS_no_HWI_Bside, 'g-.', 'LineWidth', 3);
plot(moving_distance, SE_approximation,    'r+',  'LineWidth', 2.5, 'MarkerSize', 14);

xlabel('Vehicle Moving Distance (m)', 'Interpreter', 'LaTex');
ylabel('Spectral Efficiency $\mathrm{SE}(t)$ (bit/s/Hz)', 'Interpreter', 'LaTex');
legend({'w/o RIS', ...
        'HST Side w/ $\phi_{m,n}^\mathrm{opt}(t)$', ...
        'Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$', ...
        'Analytical Result in Eq. (20)'}, ...
        'Interpreter', 'LaTex');
axis([4600, 5400, 5, 20]);

% Figure 2: Delay spread vs distance
figure(2); hold on; box on; grid on;

b = load('Td_1RIS_Bside');                       Td_1RIS_Bside = b.Td_1RIS_Bside;
c = load('Td_1RIS_Vside_k_imin_plus');           Td_1RIS_Vside_k_imin_plus = c.Td_1RIS_Vside;
d = load('Td_1RIS_Vside_vehicle_length_minus3'); Td_1RIS_Vside_vehicle_length_minus3 = d.Td_1RIS_Vside;

plot(moving_distance, Td_1RIS_Vside,                      'b-',  'LineWidth', 3);
plot(moving_distance, Td_1RIS_Vside_vehicle_length_minus3, 'g-.', 'LineWidth', 3);
plot(moving_distance, Td_1RIS_Vside_k_imin_plus,           'r--', 'LineWidth', 3);
plot(moving_distance, Td_1RIS_Bside,                       'k:',  'LineWidth', 3);

xlabel('Vehicle Moving Distance (m)', 'Interpreter', 'LaTex');
ylabel('Delay Spread $T(t)$ (s)', 'Interpreter', 'LaTex');
legend({'HST Side w/ $\phi_{m,n}^\mathrm{opt}(t), \Delta=3 \ \mathrm{m}$', ...
        'HST Side w/ $\phi_{m,n}^\mathrm{opt}(t), \Delta=-3 \ \mathrm{m}$', ...
        'HST Side w/ $\phi_{m,n}^*(t)$', ...
        'Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$'}, ...
        'Interpreter', 'LaTex', 'Location', 'NorthEast');
axis([4e3, 6e3, -0.2e-7, 1.6e-7]);

% Figure 3: CDF of SE
figure(3); hold on; box on; grid on;

plot(sort(SE_no_IRS(:),          'ascend'), linspace(0,1,total_time), 'b--', 'LineWidth', 3);
plot(sort(SE_IRS_no_HWI_Bside(:),'ascend'), linspace(0,1,total_time), 'k',   'LineWidth', 3);
plot(sort(SE_IRS_no_HWI_Vside(:),'ascend'), linspace(0,1,total_time), 'r-.', 'LineWidth', 3);
plot(linspace(0,14,15), 0.05*ones(15,1), 'k--', 'LineWidth', 1);

xlabel('Spectral Efficiency $\mathrm{SE}(t)$ (bit/s/Hz)', 'Interpreter', 'LaTex');
ylabel('CDF', 'Interpreter', 'latex');
legend({'w/o RIS', ...
        'Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$', ...
        'HST Side w/ $\phi_{m,n}^\mathrm{opt}(t)$'}, ...
        'Interpreter', 'latex', 'Location', 'SouthEast');

toc;