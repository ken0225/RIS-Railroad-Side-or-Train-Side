function [output_function_rotate_IRS]=function_rotate_IRS(rotate_angle)
% author: Ke(Ken)WANG from Macao Polytechnic Institute
% email: p1909883@ipm.edu.mo
% update information: v0.1(2020/11/19)
% 
% This function rotates the 3D IRS coordinate counterclockwise around
% the Y-axis by rotate_angle radians.
%
% The INPUT is the angle of rotation and the OUTPUT is a 3x3 rotation matrix.
%

    
    % Counterclockwise rotation around the Y-axis (when rotate_angle > 0)
    output_function_rotate_IRS = [cos(rotate_angle), 0, sin(rotate_angle); 0, 1, 0; -sin(rotate_angle), 0, cos(rotate_angle)]';

end


% Example Usage:
% IRS = [10, 10, 0]; IRS * function_rotate_IRS(pi/4)
% 
% ans =
% 
%     7.0711   10.0000   -7.0711