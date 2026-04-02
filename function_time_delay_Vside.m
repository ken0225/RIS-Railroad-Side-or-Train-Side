function [tau_0, tau_mn]=function_time_delay_Vside(p_IRS_trajectory, p_BS, p_User_trajectory)
% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/11/20), v0.2(2021/08/26), v0.3(2021/09/04)
%
% This function aims to calculate tau_0 and tau_mn, i.e., eq(5) and eq(7), in the GC paper 2021
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example:
%
% TBD

global c0 total_time;

% if ismatrix(p_IRS_trajectory) && ismatrix(p_BS) && ismatrix(p_User_trajectory)
%
%     %Check the center_IRS is suitable or not
%     centers_IRS = function_check_dim(centers_IRS);
%
%     %Check the coordinate of BS is suitable or not
%     p_BS = function_check_dim(p_BS);
%
%     %Check the coordinate of the trajectory of vehicle is suitable or not
%     p_User_trajectory = function_check_dim(p_User_trajectory);
%
% else
%
%     error('Only matrix supported! ')
%
% end

%tau_0
tau_0 = vecnorm((p_User_trajectory-p_BS).') ./ c0;


% 计算t_mn的话应该需要两个for循环嵌套，首先是以时间tt = 1：total_time为变量的外部for循环
% 之后是以ii = 1 : MN为变量的内部for循环
% 注意size(p_IRS_trajectory) = MN x 3 x total_time， 那么我们首先要用
% p_IRS_trajectory(:,:,tt)

tau_mn = [];

for tt = 1 : total_time
    
    tau_mn_the_ii_ele_every_tt = [];
    
    for ii = 1 : size(p_IRS_trajectory,1) % size(p_IRS_trajectory,1) = MN
        
        %Calculate the propagation time between the i-th element of the IRS
        %and the BS.
        temp1 = vecnorm((p_IRS_trajectory(ii,:,tt)-p_BS).') ./ c0;
        
        %Calculate the propagation time between the user and the i-th element of the IRS
        %every second
        temp2 = vecnorm((p_User_trajectory(tt,:)-p_IRS_trajectory(ii,:,tt)).') ./ c0; % size(p_User_trajectory(1,:)) = 1 x 3, size(p_IRS_trajectory(ii,:,tt)) = 1 x 3
        
        tau_mn_the_ii_ele_every_tt(ii) = temp1+temp2; % 第tt秒中的第ii个元素所造成的时延
        
    end
    
    tau_mn = [tau_mn; tau_mn_the_ii_ele_every_tt];
    
    clear tau_mn_the_ii_ele_every_tt
    
end

tau_mn = tau_mn';

end
