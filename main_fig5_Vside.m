% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1 (Jan 07, 2022)


%% Clean All & Timer
close all;
clear;
tic;

%% System Parameters Initialization

global c0 fc speed M N total_time; % Global parameters

% -------- Changeable parameters begin -------- 
% You CAN change the parameters below to obtain different figures
% In this simulation we IGNORE the phase drifts


M = 64; N = M; % IRS has M x N elements

%  -------- Changeable parameters end -------- 

c0 = 299792458; % Light speed

fc= 2.4e9; % Carrier frequency is 2.4 GHz

disp(['Step 1: Parameter Initialization begin. In this setup, the IRS has ' num2str(M) ' x ' num2str(N) ' elements.']);

dx = (c0/fc)/(2*sqrt(pi)); dy = dx; % dxdy=\lambda^2/(4\pi), this means the element gain is 1. 
% RIS的总面积为dx*dy*32^2

h_IRS = 3; % Height of the IRS

h_User = 1.5; % Height of the vehicle

h_BS = 20; % Height of the BS

rotate_angle = pi/4;

vehicle_length = 3; % 高铁：8节车厢一共是220米，我们设每个车厢一个IRS，那么就是27.5米一个车厢
%小轿车：1.5 from RIS_HighMobility_CE_ICC_RZ20

speed = 55; % Vehicle speed, note that 400 km/h \approx 110 m/s

total_time = 200; % Total time for one vehicle pass

Pt = 0.1; % Transmit power is 1W = 20dBm

N0 = 10e-12; % N0 = -80dBm

% The coordinate of the vehicle, i.e., the receiver
x_User_start = -5000;
y_User = 0+h_User;
z_User = 50; 

% The coordinate of the Base Station, i.e., the transmitter
x_BS = 0;
y_BS = 0+h_BS;
z_BS = 20; 

%% Compute the trajectory for one vehicle pass

disp(['Step 2: Compute the trajectory for one vehicle pass. The speed is ' num2str(speed) ' m/s.']);

p_BS = [x_BS, y_BS, z_BS]; % size(p_BS) = 1 x 3

p_User_start = [x_User_start, y_User, z_User];

% Total trajectory = starting point + moving trajectory per second (we start from 1s rather than 0s)
p_User_trajectory = [];

% The vehicle starts at p_vehicle_start, and travels at the speed during
% total_time period
for t = 1 : total_time
    
    temp_p_IRS_trajectory =  p_User_start+function_vehicle_moving_xdir(speed, t);
    
    p_User_trajectory = [p_User_trajectory; temp_p_IRS_trajectory]; % size(p_vehicle_trajectory) = total_time x 3

end

% 'A_0' is the gain for the direct path. It is a theoretical result.
[A_0] = function_A0(p_BS, p_User_trajectory);

%% Compute the gain when the IRS is non-isotropic, after phase optimization, and w/o HWI
disp('Step 3: Compute the gain when the IRS is non-isotropic, after phase optimization, and w/o HWI.');

% 处理IRS矩阵，注意现在要计算的是每个时间t中的IRS坐标
% 那么最终的结果应该是size(centers_IRS_trajector) = MN x 3 

% isotropic = 0; % 1/0 stands for isotropic IRS/practical IRS
% 
% % Compute some parameters for IRS
% [M_new, N_new, dx_new, dy_new] = function_parameters_IRS(M, N, dx, dy, isotropic);

% Locate the positions for M and N
[centers_IRS] = function_centers_IRS_zy(M, N, dx, dy);

% 'h_IRS' is the height of IRS
centers_IRS(:, 2) = centers_IRS(:, 2) + h_IRS;

centers_IRS = centers_IRS * function_rotate_IRS(rotate_angle);

p_IRS_start = [x_User_start + vehicle_length, y_User, z_User]; %x_User_start = -5000;

p_IRS_trajectory = zeros(M*N, 3, total_time); 
%计算在total_time中，每秒的IRS移动的坐标；有MN个元素，每个元素都有一个三维坐标
% size(p_IRS_trajectory) = MN x 3 x total_time

for t = 1 : total_time
    
    % 注意在移动中只有x轴坐标在变，y和z是不变的
    temp_p_IRS_trajectory =  p_IRS_start + centers_IRS + function_vehicle_moving_xdir(speed, t);
    
    p_IRS_trajectory(:,:,t) = temp_p_IRS_trajectory;
    
end

%%220106
% 'A_mn_theoretical' is the gain for the IRS path, and it is a theoretical result.
[A_mn, A_mn_matrix] = function_Amn_Vside(p_BS, p_User_trajectory, p_IRS_trajectory, total_time);

%%
% Compute the delays
[tau_0, tau_mn] = function_time_delay_Vside(p_IRS_trajectory, p_BS, p_User_trajectory); 


% % Optimal phase shift
% phi_optimal = 2 * pi * (fc * (tau_0 - tau_mn) + ceil(-(fc * (tau_0 - tau_mn))));

k_imin = ceil(-(fc * (tau_0-tau_mn))); 

% delay_k_imin = 6;
% 
% k_imin = max(k_imin(:))+delay_k_imin;

phi_optimal = 2*pi*(fc*(tau_0-tau_mn)+k_imin);

% Final phase part
exp_phi_optimal = exp(-1i*2*pi*fc.*tau_mn-1i.*phi_optimal);

A_mn_plus_A_0 = A_mn+A_0;

%% Td

Td_1RIS_Vside = max(tau_mn + (phi_optimal ./ (2*pi*fc)))-tau_0;
%save 'Td_1RIS_Vside_k_imin_plus' Td_1RIS_Vside
%save 'Td_1RIS_Vside_vehicle_length_minus3' Td_1RIS_Vside

%% Compute the  Gain/SE Approximation, i.e., Theorem 1
disp('Step 5: Compute the gain/SE Approximation, i.e., Theorem 1');

a_HWI = 0;
kappa_t = 0; 
kappa_r = 0; 
A_star = [];
Q = [];

% Compute the A*
for t = 1 : total_time
    
    temp_A_star=function_Astar(A_mn_matrix(t,:));
    A_star = [A_star; temp_A_star];
    
end

% Compute the Q
for t = 1 : total_time
    
    %temp_Q = function_Q(A_0(:,t), A_star(t,:), A_mn_matrix(t,:), a_HWI);
    temp_Q = function_Q(A_0(t), A_star(t), A_mn_matrix(t,:), a_HWI);
    %Q = [Q; temp_Q];
    Q(t) = temp_Q;
    
end

Q = Pt*Q;

% Compute the SNR approximation
SNR_approximation = (Q ./ ((kappa_t+kappa_r)*Q + N0))'; 

% Compute the SE approximation
SE_approximation = log2(1+SNR_approximation);

%% Plot the simulation results
close all;

moving_distance = linspace(1,total_time, total_time) .* speed;

disp('Step 7: Plot the simulation results.');

inst_received_power_no_IRS = Pt * (abs(A_0)) .^ 2;
inst_SNR_no_IRS = inst_received_power_no_IRS ./ (N0);
SE_no_IRS = log2(1 + inst_SNR_no_IRS);

inst_received_power_IRS_theoretical = Pt * (abs(A_mn_plus_A_0)) .^ 2;
inst_SNR_IRS_theoretical = inst_received_power_IRS_theoretical ./ (N0);
SE_IRS_no_HWI_Vside = log2(1 + inst_SNR_IRS_theoretical);

% inst_received_power_isotropic_IRS = Pt * (abs(A_mn_total)) .^ 2;
% inst_SNR_isotropic_IRS = inst_received_power_isotropic_IRS ./ ((N0));
% SE_isotropic_IRS = log2(1 + inst_SNR_isotropic_IRS);

%% Figure 1
figure(1); hold on; box on; grid on;

a = load('SE_IRS_no_HWI_Bside'); SE_IRS_no_HWI_Bside = a.SE_IRS_no_HWI;

plot(moving_distance, SE_no_IRS,'k:','LineWidth',3);
plot(moving_distance, SE_IRS_no_HWI_Vside,'b-','LineWidth',3);
plot(moving_distance, SE_IRS_no_HWI_Bside,'g-.','LineWidth',3);
plot(moving_distance, SE_approximation,'r+','LineWidth',2.5, 'MarkerSize',14);


xlabel('Vehicle Moving Distance (m)','Interpreter','LaTex');
ylabel('Spectral Efficiency $\mathrm{SE}(t)$ (bit/s/Hz)','Interpreter','LaTex');
%title('Spectral Efficiency for One Vehicle Pass','Interpreter','LaTex');
legend({'w/o RIS','HST Side w/ $\phi_{m,n}^\mathrm{opt}(t)$','Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$','Analytical Result in Eq. (20)'},'Interpreter','LaTex');
axis([4600, 5400, 5, 20]);

%% Figure 2
figure(2); hold on; box on; grid on;

b = load('Td_1RIS_Bside'); Td_1RIS_Bside = b.Td_1RIS_Bside;
c = load('Td_1RIS_Vside_k_imin_plus'); Td_1RIS_Vside_k_imin_plus = c.Td_1RIS_Vside;
d = load('Td_1RIS_Vside_vehicle_length_minus3'); Td_1RIS_Vside_vehicle_length_minus3 = d.Td_1RIS_Vside;

plot(moving_distance, Td_1RIS_Vside,'b-','LineWidth',3);
plot(moving_distance, Td_1RIS_Vside_vehicle_length_minus3,'g-.','LineWidth',3);
plot(moving_distance, Td_1RIS_Vside_k_imin_plus, 'r--','LineWidth',3);
plot(moving_distance, Td_1RIS_Bside,'k:','LineWidth',3);


xlabel('Vehicle Moving Distance (m)','Interpreter','LaTex');
ylabel('Delay Spread $T(t)$ (s)','Interpreter','LaTex');
%title('Delay Spread for One Vehicle Pass','Interpreter','LaTex');
legend({'HST Side w/ $\phi_{m,n}^\mathrm{opt}(t), \Delta=3 \ \mathrm{m}$','HST Side w/ $\phi_{m,n}^\mathrm{opt}(t), \Delta=-3 \ \mathrm{m}$'...
    'HST Side w/ $\phi_{m,n}^*(t)$','Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$'},'Interpreter','LaTex','Location','NorthEast');
axis([4e3, 6e3, -0.2e-7, 1.6e-7]);

%% Figure 3
figure(3); hold on; box on; grid on;

plot(sort(SE_no_IRS(:),'ascend'),linspace(0,1,total_time),'b--','LineWidth',3);
plot(sort(SE_IRS_no_HWI_Bside(:),'ascend'),linspace(0,1,total_time),'k','LineWidth',3);
plot(sort(SE_IRS_no_HWI_Vside(:),'ascend'),linspace(0,1,total_time),'r-.','LineWidth',3);
plot(linspace(0,14,15),0.05*ones(15,1),'k--','LineWidth',1);

xlabel('Spectral Efficiency $\mathrm{SE}(t)$ (bit/s/Hz)','Interpreter','LaTex');
ylabel('CDF','Interpreter','latex');
%title('CDF of Spectral Efficiency for One Vehicle Pass','Interpreter','LaTex');
legend({'w/o RIS','Railroad Side w/ $\phi_{m,n}^\mathrm{opt}(t)$','HST Side w/ $\phi_{m,n}^\mathrm{opt}(t)$'},'Interpreter','latex','Location','SouthEast');

% % Figure 4
% figure(4); hold on; box on; grid on;
% 
% plot(sort(Td_1RIS_Vside(:),'ascend'),linspace(0,1,total_time),'b--','LineWidth',3);
% plot(sort(Td_1RIS_Vside_K100(:),'ascend'),linspace(0,1,total_time),'k','LineWidth',3);
% plot(sort(Td_1RIS_Bside(:),'ascend'),linspace(0,1,total_time),'r-.','LineWidth',3);
% 
% xlabel('Delay Spread $T_d(t)$ (s)','Interpreter','LaTex');
% ylabel('CDF','Interpreter','latex');
% legend({'Vehicle-Side w/ $(19)$','Vehicle-Side w/ $k_{i_m,i_n}=100$','Roadside'},'Interpreter','LaTex','Location','NorthEast');

% % Figure 5
% 
% e = load('EE_IRS_no_HWI_Bside'); EE_IRS_no_HWI_Bside = e.EE_IRS_no_HWI_Bside;
% 
% B = 1e6; % Bandwidth is 1MHz
% 
% Pe = 1e-3; % The power consumption of each phase shifter (element), the typical value is 0.8 dBm
% 
% nu = 0.5; % The efficiency of the transmit power amplifier, the typical value is 0.5
% 
% PU =  3e-3; % The hardware static power in the receiver, the typical value is 5 dBm
% 
% PBS = 3e-3; % The total hardware static power cunsumption at BS, the typical value is 5 dBm
% 
% EE_IRS_no_HWI_Vside = (B .* SE_IRS_no_HWI_Vside) ./ (Pt/nu + PU + PBS + M*N*Pe);
% EE_no_IRS = (B .* SE_no_IRS) ./ (Pt/nu + PU + PBS + M*N*Pe);
% 
% figure(5); hold on; box on; grid on;
% 
% plot(1 : total_time, EE_IRS_no_HWI_Bside./1e6,'b--','LineWidth',3);
% plot(1 : total_time, EE_IRS_no_HWI_Vside./1e6,'k','LineWidth',3);
% 
% 
% xlabel('Vehicle Moving Time (s)','Interpreter','LaTex');
% ylabel('Energy Efficiency (Mbit/Joule)','Interpreter','LaTex');
% title('Energy Efficiency for One Vehicle Pass','Interpreter','LaTex');

% % Figure 6
% figure(6); hold on; box on; grid on;
% 
% plot(sort(EE_no_IRS(:)./1e6,'ascend'),linspace(0,1,total_time),'b--','LineWidth',3);
% plot(sort(EE_IRS_no_HWI_Bside(:)./1e6,'ascend'),linspace(0,1,total_time),'k','LineWidth',3);
% plot(sort(EE_IRS_no_HWI_Vside(:)./1e6,'ascend'),linspace(0,1,total_time),'r-.','LineWidth',3);
% plot(linspace(0,5,6),0.05*ones(6,1),'k--','LineWidth',1);
% 
% xlabel('Energy Efficiency (Mbit/Joule)','Interpreter','LaTex');
% ylabel('CDF','Interpreter','latex');
% %title('CDF of Energy Efficiency for One Vehicle Pass','Interpreter','LaTex');
% legend({'w/o RIS','Roadside w/ $\phi_{m,n}^\mathrm{opt}(t)$ [6]','Vehicle Side w/ $\phi_{m,n}^\mathrm{opt}(t)$'},'Interpreter','latex','Location','SouthEast');

toc;





