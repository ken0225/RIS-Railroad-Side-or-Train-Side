function [Q]=function_Q(A_0, A_star, A_mn_vector, a)

Q = (A_0)^2 + (sinc(a/pi))^2 * A_star + sum(A_mn_vector.^2) + 2*A_0*sinc(a/pi)*sum(A_mn_vector); 

end