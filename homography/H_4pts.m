function [ C_H ] = H_4pts( p1, p2, H, C )
%H_4pts - propagate the uncertainty from 4 pts into Homography matrix paramers
% Input: 
%   p1  ...  2x4pts [u11 u12 u13 u14; v11 v12 v13 v14] (first image)
%   p2  ...  2x4pts [u21 u22 u23 u24; v21 v22 v23 v24] (second image)
%   H   ... the homography matrix
%   C   ... 16x16 covariance matrix of input points [u11 v11 u12 v12 u13 v13 u14 v14   u21 v21 u22 v22 u23 v23 u24 v24]
% Output:
%   C_H ... 9x9 covariance matrix of homography parameters [h11 h21 h31 h12 h22 h32 h13 h23 h33]


    % assume h33 = 1
    H = normalizeMatrix(H);
    
    % compute the derivatives
    A = deriv_measurements( p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),...
                            p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p2(7),p2(8),...
                            H(1),H(2),H(3),H(4),H(5),H(6),H(7),H(8),H(9));
    B = deriv_params(   p1(1),p1(2),p1(3),p1(4),p1(5),p1(6),p1(7),p1(8),...
                            p2(1),p2(2),p2(3),p2(4),p2(5),p2(6),p2(7),p2(8),...
                            H(1),H(2),H(3),H(4),H(5),H(6),H(7),H(8),H(9));
    % the propagation
    iBA = B \ A;
    C_H = iBA * C * iBA';
end

function A = deriv_measurements(u11, v11, u12, v12, u13, v13, u14, v14, ...
                                u21, v21, u22, v22, u23, v23, u24, v24,...
                                h11, h21, h31, h12, h22, h32, h13, h23, h33)
    A = [-h31 * u21 + h11 -h32 * u21 + h12 0 0 0 0 0 0 -h31 * u11 - h32 * v11 - h33 0 0 0 0 0 0 0; -h31 * v21 + h21 -h32 * v21 + h22 0 0 0 0 0 0 0 -h31 * u11 - h32 * v11 - h33 0 0 0 0 0 0; 0 0 -h31 * u22 + h11 -h32 * u22 + h12 0 0 0 0 0 0 -h31 * u12 - h32 * v12 - h33 0 0 0 0 0; 0 0 -h31 * v22 + h21 -h32 * v22 + h22 0 0 0 0 0 0 0 -h31 * u12 - h32 * v12 - h33 0 0 0 0; 0 0 0 0 -h31 * u23 + h11 -h32 * u23 + h12 0 0 0 0 0 0 -h31 * u13 - h32 * v13 - h33 0 0 0; 0 0 0 0 -h31 * v23 + h21 -h32 * v23 + h22 0 0 0 0 0 0 0 -h31 * u13 - h32 * v13 - h33 0 0; 0 0 0 0 0 0 -h31 * u24 + h11 -h32 * u24 + h12 0 0 0 0 0 0 -h31 * u14 - h32 * v14 - h33 0; 0 0 0 0 0 0 -h31 * v24 + h21 -h32 * v24 + h22 0 0 0 0 0 0 0 -h31 * u14 - h32 * v14 - h33; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

function B = deriv_params(  u11, v11, u12, v12, u13, v13, u14, v14, ...
                            u21, v21, u22, v22, u23, v23, u24, v24,...
                            h11, h21, h31, h12, h22, h32, h13, h23, h33)
    B = [u11 0 -u11 * u21 v11 0 -u21 * v11 1 0 -u21; 0 u11 -u11 * v21 0 v11 -v11 * v21 0 1 -v21; u12 0 -u12 * u22 v12 0 -u22 * v12 1 0 -u22; 0 u12 -u12 * v22 0 v12 -v12 * v22 0 1 -v22; u13 0 -u13 * u23 v13 0 -u23 * v13 1 0 -u23; 0 u13 -u13 * v23 0 v13 -v13 * v23 0 1 -v23; u14 0 -u14 * u24 v14 0 -u24 * v14 1 0 -u24; 0 u14 -u14 * v24 0 v14 -v14 * v24 0 1 -v24; (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h11 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h21 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h31 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h12 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h22 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h32 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h13 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h23 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h33];
end