function [output_centers_IRS]=function_centers_IRS_zy(M, N, dx, dy)
% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update information: v0.1(2020/11/19), v0.2(2021/08/26)
%
% This function aims to calculate an IRS with M x N elements, and each of
% which with the size of dx x dy. Note that before using this function, you
% have to know if the IRS is isotropic or not, and ensure the inputs are correct.
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% original article.
%
% Example: 
%
% [output_function_centers_IRS]=function_centers_IRS(2, 2, 1, 2)
%
%    -0.5000   -1.0000         0
%    -0.5000    1.0000         0
%     0.5000   -1.0000         0
%     0.5000    1.0000         0

M_interval = function_interval(M);

N_interval = function_interval(N);

M_range = (M_interval(1) : M_interval(2));

N_range = (N_interval(1) : N_interval(2));

for m = 1 : length(M_range)
    
    gm(m) = M_range(m)*dx - 0.5*dx*mod(M+1, 2);
    
end

for n = 1 : length(N_range)
    
    gn(n) = N_range(n)*dy - 0.5*dy*mod(N+1, 2);
    
end

output_centers_IRS = function_cartprod(gm, gn);

output_centers_IRS(:, end+1) = 0;

output_centers_IRS = sortrows(output_centers_IRS, 1);

output_centers_IRS(:, 3) = output_centers_IRS(:, 1); % Convert X-axis coordinates to Z-axis

output_centers_IRS(:, 1) = 0; % Then set X-axis coordinates to 0

end