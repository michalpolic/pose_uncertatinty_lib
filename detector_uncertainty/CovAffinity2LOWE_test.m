% Unity test of: CovAffinity2LOWE
% Input
s1 =  2.997189760208130;    s2 = 3.096155166625977;
a1 = -0.448401885418503;    a2 = -0.510097079744290;
d1 = 1;     d2 = 1;
f1 = 1;     f2 = 1;
sigma_s1 = 0.1*s1;         
sigma_s2 = 0.1*s2;
sigma_r1 = 10 / 180 * pi;
sigma_r2 = 10 / 180 * pi;
sigma_d1 = 0.1;
sigma_d2 = 0.1;
sigma_f1 = 0.1;
sigma_f2 = 0.1;
in_cov_sr = diag([sigma_r1, sigma_r2, sigma_s1, sigma_s2]).^2;
in_cov_srdf = diag([sigma_r1, sigma_r2, sigma_s1, sigma_s2, ...
                    sigma_d1, sigma_d2, sigma_f1, sigma_f2]).^2;
%A = [1.014441919086481, 0.062665719626614;...
%    -0.062665719626614, 1.014441919086481];

%% Test parameters 
k     = 1e4;
sigma = 1e-7;

%% Test for rotation and scale
% the assumed covariance
C_A1 = CovAffinity2LOWE( in_cov_sr, a1, a2, s1, s2 ) * sigma^2;

% compute A with deviations 
change1 = zeros(4,k);
for i = 1:k
    a = (a2 + randn(1) * sigma_r2 * sigma) - (a1 + randn(1) * sigma_r1 * sigma);
    s = (s2 + randn(1) * sigma_s2 * sigma) / (s1 + randn(1) * sigma_s1 * sigma);
    A_i = s * [cos(a) -sin(a); sin(a) cos(a)];
    change1(:,i) = A_i(:);
end 

% the result
C_A1 
C_A1_MC=cov(change1')

H = null(C_A1);
J = null(H');
C_A1_r    = J' * C_A1    * J;
C_A1_MC_r = J' * C_A1_MC * J;

[L,T,~]=check_CovM(C_A1_MC_r,C_A1_r,k,0.999);
if L > T(1) && L < T(2)
    display(['No reason to doubt covariance matrices are the same: ',...
        num2str(L),' in [',num2str(T(1)),',',num2str(T(2)),']']);
else
    display(['Reason to doubt covariance matrices are the same: ',...     
        num2str(L),' not in [',num2str(T(1)),',',num2str(T(2)),'] ********']); 
end

%% Test for rotation, scale, shear, scale difference
% the assumed covariance
C_A2 = CovAffinity2LOWE( in_cov_srdf, a1, a2, s1, s2, d1, d2, f1, f2 ) * sigma^2;

% compute A with deviations 
change2 = zeros(4,k);
for i = 1:k
    a = (a2 + randn(1) * sigma_r2 * sigma) - (a1 + randn(1) * sigma_r1 * sigma);
    s = (s2 + randn(1) * sigma_s2 * sigma) / (s1 + randn(1) * sigma_s1 * sigma);
    d = (d2 + randn(1) * sigma_d2 * sigma) - (d1 + randn(1) * sigma_d1 * sigma);
    f = (f2 + randn(1) * sigma_f2 * sigma) - (f1 + randn(1) * sigma_f1 * sigma);
    A_i = s * [exp(d) 0 ; 0 exp(-d)] * [1 f; f 1] * [cos(a) -sin(a); sin(a) cos(a)];
    change2(:,i) = A_i(:);
end 

% the result
C_A2_MC=cov(change2')
C_A2


H = null(C_A2);
J = null(H');
C_A2_r    = J' * C_A2    * J;
C_A2_MC_r = J' * C_A2_MC * J;

[L,T,~]=check_CovM(C_A2_MC_r,C_A2_r,k,0.999);
if L > T(1) && L < T(2)
    display(['No reason to doubt covariance matrices are the same: ',...
        num2str(L),' in [',num2str(T(1)),',',num2str(T(2)),']']);
else
    display(['Reason to doubt covariance matrices are the same: ',...     
        num2str(L),' not in [',num2str(T(1)),',',num2str(T(2)),'] ********']); 
end

%% Test for rotation, scale zero shear, scale difference
% the assumed covariance
C_A3 = CovAffinity2LOWE_modified( in_cov_sr, a1, a2, s1, s2, 1 ) * sigma^2;


% compute A with deviations 
change2 = zeros(4,k);
std_df = sqrt(trace(in_cov_sr)/4);
for i = 1:k
    a = (a2 + randn(1) * sigma_r2 * sigma) - (a1 + randn(1) * sigma_r1 * sigma);
    s = (s2 + randn(1) * sigma_s2 * sigma) / (s1 + randn(1) * sigma_s1 * sigma);
    d = (randn(1) * std_df - randn(1) * std_df) * sigma;
    f = (randn(1) * std_df - randn(1) * std_df) * sigma;
    A_i = s * [exp(d) 0 ; 0 exp(-d)] * [1 f; f 1] * [cos(a) -sin(a); sin(a) cos(a)];
    change3(:,i) = A_i(:);
end 

% the result
C_A3_MC=cov(change3')
C_A3


H = null(C_A3);
J = null(H');
C_A3_r    = J' * C_A3    * J;
C_A3_MC_r = J' * C_A3_MC * J;

[L,T,~]=check_CovM(C_A3_MC_r,C_A3_r,k,0.999);
if L > T(1) && L < T(2)
    display(['No reason to doubt covariance matrices are the same: ',...
        num2str(L),' in [',num2str(T(1)),',',num2str(T(2)),']']);
else
    display(['Reason to doubt covariance matrices are the same: ',...     
        num2str(L),' not in [',num2str(T(1)),',',num2str(T(2)),'] ********']); 
end








