% test CovM
%
% [lambda,T,F] = check_CovM(C_empir,C_theor,m,S)
%
% C_empir = empirical covariance matrix
% C_theor = theoretical covariance matrix
% m       = sample size for determining empirical CovM
% S       = significance level (close to 1)
%
% lambda  = test statistic
% T       = chi^2 threshold
% F       = refined F threshold
%
% Wolfgang Förstner 2/2012
% wfoerstn@uni-bonn.de 

function [lambda, T, F] = check_CovM(C_empir, C_theor, m, S)

p = size(C_empir, 1);

lambda = m * (log(det(C_theor) / det(C_empir)) - p + trace(C_empir * inv(C_theor))); %#ok<MINV>

Tup = chi2inv(1 - (1 - S) / 2, 0.5 * p * (p + 1));
Tlow = chi2inv((1 - S) / 2, 0.5 * p * (p + 1));
T = [Tlow, Tup];

D1 = (2 * p + 1 - 2 / (p + 1)) / (6 * m);
D2 = (p - 1) * (p + 2) / (6 * m ^ 2);
q1 = p * (p + 1) / 2;
q2 = (q1 + 2) / (D2 - D1 ^ 2);
b = q1 / (1 - D1 - q1 / q2);
q2 = min(10000, q2);
Fup = b * finv(1 - (1 - S) / 2, q1, q2);
Flow = b * finv((1 - S) / 2, q1, q2);
F = [Flow, Fup];



