%% b and R from E using set of pairs of directions
%
% Alg. 20
%
% [b, R, code] = bR_from_E_uv(E,u,v)
%
% E     = 3x3 E-matrix
% u,v   = Nx3 matrices of directions
%
% b     = 3x1 normalized base vector
% R     = 3x3 rotation matrix
% code  = code for signes of b and R
%           0 S(b)= UZU', R = UWV'
%           1 S(b)= UZ'U', R = UWV
%           2 S(b)= UZU', R = UW'V'
%           3 S(b)= UZ'U', R = UW'V'
%
% Wolfgang Förstner 8/2013
% Modified by Michal Polic
% wfoerstn@uni-bonn.de
%
% See also sugr_E_Matrix

function [b, R, aa, code] = bR_from_E_uv(E, u, v)
 
    % set basic matrices
    W = [0 1 0; - 1, 0, 0; 0, 0, 1];
    Z = [0 1 0; - 1, 0, 0; 0, 0, 0];
    I = size(u, 1);

    % svd of E
    [U, S, V] = svd(E);
    % enforce U and V to be proper rotations
    U = U * det(U);
    V = V * det(V);
    % check sign of b and R
    code = - 1;
    count = 0;
    for c = 0:3 % for all four cases set b and R

         count = count + 1;
         switch c
             case 0
                 S = U*Z'*U';   R = U*W*V';
             case 1
                 S = U*Z*U';    R = U*W*V';
             case 2
                 S = U*Z'*U';   R = U*W'*V';
             case 3
                 S = U*Z*U';    R = U*W'*V';
         end
         b = [S(3, 2); S(1, 3); S(2, 1)];
         b = b / sqrt(b'*b);
         % check whether all 3D points are in direction of u and v
         sign_s = zeros(I, 1);
         sign_r = zeros(I, 1);
         for i = 1:I
             ui = u(i, :)'; % * (1/norm(u(i, :)))
             vi = v(i, :)'; %* (1/norm(v(i, :)))
             wi = R*vi;
             m = cross(cross(b, ui), b);
             sign_s(i) = sign(det([b, m, cross(ui, wi)]));
             sign_r(i) = sign_s(i) * sign(m'*wi);
         end
         % compute the angle * axis
         axis = X2v(R - R');
         if sqrt(axis' * axis) < 100 * eps
             axis = null(R - eye(3));
         end
         axis = axis(:,1) * (1/sqrt(axis(:,1)'*axis(:,1)));
         aa = acos(0.5*(trace(R) - 1)) * axis;
           
         % check: the majority of points need to be in direction of u and v
         %     signs = [sign_s,sign_r];
         %     correct_sign = [mean(sign_s),mean(sign_r)];
         if mean(sign_s) > 0 && mean(sign_r) > 0
             code = c;
             return
         end
    end
end
