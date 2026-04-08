function [A_star]=function_Astar(A_mn_vector)

% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/11/19), v0.2(2021/08/26), v0.3(2021/09/04)

temp_A_star_2 = [];

for k = 1 : size(A_mn_vector,2)

temp_A_star_1 = A_mn_vector(k) * (sum(A_mn_vector)-A_mn_vector(k));

temp_A_star_2 = [temp_A_star_2; temp_A_star_1 ];

end

A_star = sum(temp_A_star_2);

end