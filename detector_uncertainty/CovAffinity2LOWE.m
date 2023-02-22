function [ C_A ] = CovAffinity2LOWE( cov_in, r1, r2, s1, s2, d1, d2, f1, f2)
%CovAffinity2LOWE propagates the uncertainty of LOWE keypoint into Affinity matrix
% Input: 
% cov_in ... the covariance of input parameters, 
%            if no given shear/scale difference: 
%                   diag([sigma_r1, sigma_r2, sigma_s1, sigma_s2]).^2;
%            otherwise
%                   diag([sigma_r1, sigma_r2, sigma_s1, sigma_s2,sigma_d1, sigma_d2, sigma_f1, sigma_f2]).^2;
% [r1,r2] ... rotation of the first and second Lowe keypoint
% [s1,s2] ... scale of the first and second Lowe keypoint
% [d1,d2] ... shear of the first and second Lowe keypoint
% [f1,f2] ... scale difference of the first and second Lowe keypoint
%
% Output:
% C_A ... covariance of Affinity matrix

    switch nargin
        case 9
            J = srdf_par_derivative(r1, r2, s1, s2, d1, d2, f1, f2);
        case 5
            J = sr_par_derivative(r1, r2, s1, s2);
        case 6
            % set scale differences and shears
            d1=0; d2=0; f1=0; f2=0;            
            J = srdf_par_derivative(r1, r2, s1, s2, d1, d2, f1, f2);
            % extend cov_in
            var_df = trace(cov_in)/4;
            cov = [cov_in, zeros(4); zeros(4) var_df*eye(4)];
            cov_in=cov;
            
    end

    C_A = J' * cov_in * J;
end

function J = sr_par_derivative(r1, r2, s1, s2)
    t1 = -r2 + r1;
    t2 = sin(t1);
    t1 = cos(t1);
    t3 = 1 / s1;
    t4 = t3 * t2;
    t5 = t4 * s2;
    t6 = t3 * t1;
    t7 = t6 * s2;
    t3 = s2 * t3 ^ 2;
    t1 = t3 * t1;
    t2 = t3 * t2;
    J = [-t5 -t7 t7 -t5; t5 t7 -t7 t5; -t1 t2 -t2 -t1; t6 -t4 t4 t6];
end

function J = srdf_par_derivative(r1, r2, s1, s2, d1, d2, f1, f2)
    t1 = d1 - d2;
    t2 = exp(-t1);
    t3 = -r2 + r1;
    t4 = sin(t3);
    t5 = f1 - f2;
    t3 = cos(t3);
    t6 = 0.1e1 / s1;
    t7 = t5 * t3;
    t8 = t7 - t4;
    t9 = s2 * t6;
    t10 = t9 * t2;
    t11 = t10 * t8;
    t1 = exp(t1);
    t5 = t5 * t4;
    t12 = t5 - t3;
    t9 = t9 * t1;
    t13 = t9 * t12;
    t5 = t5 + t3;
    t14 = t10 * t5;
    t7 = t7 + t4;
    t15 = t9 * t7;
    t16 = s2 * t6 ^ 2;
    t17 = t16 * t2;
    t16 = t16 * t1;
    t2 = t6 * t2;
    t1 = t6 * t1;
    t6 = t10 * t4;
    t18 = t9 * t3;
    t3 = t10 * t3;
    t4 = t9 * t4;
    J = [t11 t13 t14 -t15; -t11 -t13 -t14 t15; -t17 * t5 t16 * t7 t17 * t8 t16 * t12; t2 * t5 -t1 * t7 -t2 * t8 -t1 * t12; -t14 -t15 t11 -t13; t14 t15 -t11 t13; t6 -t18 -t3 -t4; -t6 t18 t3 t4];
end
