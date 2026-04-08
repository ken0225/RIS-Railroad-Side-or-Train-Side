function [A_mn, A_mn_matrix]=function_Amn_Vside(p_BS, p_User_trajectory, p_IRS_trajectory, total_time)
% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update information: v0.1 (Jan 07, 2022)
%
% This function aims to calculate A_mn if the IRS is on the vehicle side
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example:
%
% TBD

global c0 fc;

lambda_c = c0/fc;

%============== BS to IRS ==============
A_BS_IRS = []; % size(A_BS_IRS) = 1 x MN x total_time

for tt = 1 : total_time
    
    temp_A_BS_IRS = lambda_c/(4*pi) * sqrt(1*1)./vecnorm((p_IRS_trajectory(:,:,tt) - p_BS)');
    % size(p_IRS_trajectory(:,:,t)) = MN x 3; size(p_BS) = 1 x 3;
    % size(p_IRS_trajectory(:,:,t)-p_BS) = MN x 3
    % size(vecnorm((p_IRS_trajectory(:,:,1)-p_BS)')) = 1 x MN
    % size(lambda_c/(4*pi) * sqrt(1*1)./vecnorm((p_IRS_trajectory(:,:,1)-p_BS)')) = 1 x MN
    
    A_BS_IRS(:,:,tt) = temp_A_BS_IRS;
    
end
%============== BS to IRS ==============


%============== IRS to USER ==============
A_IRS_USER = []; % size(A_IRS_USER) = 1 x MN x total_time

for tt = 1 : total_time
    
    temp_A_IRS_USER = lambda_c/(4*pi) * sqrt(1*1)./vecnorm((p_User_trajectory(tt,:) - p_IRS_trajectory(:,:,tt))');
    % size(p_vehicle_trajectory(t,:)) = 1 x 3; size(p_IRS_trajectory(:,:,t)) = MN x 3
    
    A_IRS_USER(:,:,tt) = temp_A_IRS_USER;
    
end
%============== IRS to USER ==============

% Total gain, A_BS_IRS .* A_IRS_USER

A_mn = sum(A_BS_IRS .* A_IRS_USER);

A_mn = reshape(A_mn, [1, total_time]);

temp_A_mn_matrix = A_BS_IRS .* A_IRS_USER;


for ii = 1 : total_time
    
    A_mn_matrix(ii,:) = temp_A_mn_matrix(:,:,ii);
    
end

end