function [Q]=function_Q(A_0, A_star, A_mn_vector, a)

% Author: Ke(Ken)WANG from Macao Polytechnic Institute
% Email: ke.wang@ipm.edu.mo, kewang0225@gmail.com
% Update infomation: v0.1(2020/11/19), v0.2(2021/08/26), v0.3(2021/09/04)

Q = (A_0)^2 + (sinc(a/pi))^2 * A_star + sum(A_mn_vector.^2) + 2*A_0*sinc(a/pi)*sum(A_mn_vector); 

end