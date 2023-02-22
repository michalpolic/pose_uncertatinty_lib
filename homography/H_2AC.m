function [ C_H ] = H_2AC( p1, p2, A1, A2, H, C, sel_eqns )
%H_2AC - propagate the uncertainty from 2 affine correspondences into the
% uncertatinty of homography matrix
% Input: 
%   p1  ...  2x2pts [u11 u12; v11 v12] (first image)
%   p2  ...  2x2pts [u21 u22; v21 v22] (second image)
%   A1  ...  2x2 affinity matrix [a11_1 a12_1 a21_1 a22_1]
%   A2  ...  2x2 affinity matrix [a11_2 a12_2 a21_2 a22_2]
%   H   ...  the homography matrix
%   C   ...  16x16 covariance matrix of input points [u11 v11 u12 v12 u21 v21 u22 v22  a11_1 a12_1 a21_1 a22_1  a11_2 a12_2 a21_2 a22_2]
%   sel_eqns ... we need >=8 equations out of 12 for the uncertatinty propagation, 
%               this vector represents ids of equations that will be used for 
%               uncertatinty propagation. Note that equations [1,2,3,4] are
%               nesscessary.
% Output:
%   C_H ... 8x8 covariance matrix of homography parameters [h11 h12 h13 h21 h22 h23 h31 h32]


    % normalize H
    H = normalizeMatrix(H);
    
    % compute the derivatives
    A = deriv_measurements(p1(1),p1(2),p1(3),p1(4),...
                           p2(1),p2(2),p2(3),p2(4),...
                           A1(1),A1(2),A1(3),A1(4),...
                           A2(1),A2(2),A2(3),A2(4),...
                           H(1),H(2),H(3),H(4),H(5),H(6),H(7),H(8),H(9));
    
	B = deriv_params(p1(1),p1(2),p1(3),p1(4),...
                           p2(1),p2(2),p2(3),p2(4),...
                           A1(1),A1(2),A1(3),A1(4),...
                           A2(1),A2(2),A2(3),A2(4),...
                           H(1),H(2),H(3),H(4),H(5),H(6),H(7),H(8),H(9));

    % select the equtations used for the propagation of the ucertainty
    if nargin > 6 && ~isempty(sel_eqns) 
        if length(sel_eqns) ~= length(unique([sel_eqns 1 2 3 4])) || length(sel_eqns) < 9    
            error('H_2AC: sel_eqns require at least 8 equation ids and the ids [1 2 3 4] are mandatory.');
        else
            A = A(sel_eqns,:);
            B = B(sel_eqns,:);
        end
    end
                         
    % the propagation
    iBA = B \ A;
    C_H = iBA * C * iBA';
end

function A = deriv_measurements(u11, v11, u12, v12, ...
                                u21, v21, u22, v22,...
                                a11_1, a12_1, a21_1, a22_1,...
                                a11_2, a12_2, a21_2, a22_2,...
                                h11, h21, h31, h12, h22, h32, h13, h23, h33)
    A = [-h31 * u21 + h11 -h32 * u21 + h12 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0 0 0 0 0 0 0 0; -h31 * v21 + h21 -h32 * v21 + h22 0 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0 0 0 0 0 0 0; 0 0 -h31 * u22 + h11 -h32 * u22 + h12 0 0 -u12 * h31 - v12 * h32 - h33 0 0 0 0 0 0 0 0 0; 0 0 -h31 * v22 + h21 -h32 * v22 + h22 0 0 0 -u12 * h31 - v12 * h32 - h33 0 0 0 0 0 0 0 0; -a11_1 * h31 -a11_1 * h32 0 0 -h31 0 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0 0 0 0; -a12_1 * h31 -a12_1 * h32 0 0 -h32 0 0 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0 0 0; -a21_1 * h31 -a21_1 * h32 0 0 0 -h31 0 0 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0 0; -a22_1 * h31 -a22_1 * h32 0 0 0 -h32 0 0 0 0 0 -u11 * h31 - v11 * h32 - h33 0 0 0 0; 0 0 -a11_2 * h31 -a11_2 * h32 0 0 -h31 0 0 0 0 0 -u12 * h31 - v12 * h32 - h33 0 0 0; 0 0 -a12_2 * h31 -a12_2 * h32 0 0 -h32 0 0 0 0 0 0 -u12 * h31 - v12 * h32 - h33 0 0; 0 0 -a21_2 * h31 -a21_2 * h32 0 0 0 -h31 0 0 0 0 0 0 -u12 * h31 - v12 * h32 - h33 0; 0 0 -a22_2 * h31 -a22_2 * h32 0 0 0 -h32 0 0 0 0 0 0 0 -u12 * h31 - v12 * h32 - h33; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

function B = deriv_params(  u11, v11, u12, v12, ...
                                u21, v21, u22, v22,...
                                a11_1, a12_1, a21_1, a22_1,...
                                a11_2, a12_2, a21_2, a22_2,...
                                h11, h21, h31, h12, h22, h32, h13, h23, h33)
    B = [u11 0 -u11 * u21 v11 0 -u21 * v11 1 0 -u21; 0 u11 -u11 * v21 0 v11 -v11 * v21 0 1 -v21; u12 0 -u12 * u22 v12 0 -u22 * v12 1 0 -u22; 0 u12 -u12 * v22 0 v12 -v12 * v22 0 1 -v22; 1 0 -a11_1 * u11 - u21 0 0 -a11_1 * v11 0 0 -a11_1; 0 0 -a12_1 * u11 1 0 -a12_1 * v11 - u21 0 0 -a12_1; 0 1 -a21_1 * u11 - v21 0 0 -a21_1 * v11 0 0 -a21_1; 0 0 -a22_1 * u11 0 1 -a22_1 * v11 - v21 0 0 -a22_1; 1 0 -a11_2 * u12 - u22 0 0 -a11_2 * v12 0 0 -a11_2; 0 0 -a12_2 * u12 1 0 -a12_2 * v12 - u22 0 0 -a12_2; 0 1 -a21_2 * u12 - v22 0 0 -a21_2 * v12 0 0 -a21_2; 0 0 -a22_2 * u12 0 1 -a22_2 * v12 - v22 0 0 -a22_2; (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h11 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h21 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h31 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h12 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h22 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h32 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h13 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h23 (h11 ^ 2 + h12 ^ 2 + h13 ^ 2 + h21 ^ 2 + h22 ^ 2 + h23 ^ 2 + h31 ^ 2 + h32 ^ 2 + h33 ^ 2) ^ (-0.1e1 / 0.2e1) * h33];
end