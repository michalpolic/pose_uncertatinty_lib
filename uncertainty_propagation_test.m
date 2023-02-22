% The summary of the uncertainty propagation unit tests
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% init
clear; close; clc;
rng(1)
num_simulations = 500;
num_iner_iterations = 100;
sigma2 = 10^(-13);

% add that folder plus all subfolders to the path.
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

test_results = struct();

% %% HOMOGRAPHY
H_4pts_test; test_results.H_4PC = res;
H_2AC_test; test_results.H_2AC = res;
% 
% %% ESSENTIAL MATRIX
E_5pts_test_br; test_results.E_5PC = res;
E_2AC_test_br; test_results.E_2PC = res;
% 
% %% ESSENTIAL MATRIX + focal
Ef_6pts_test; test_results.Ef_6PC = res;
Ef_2AC_test; test_results.Ef_2AC = res;
% Ef_2AC_br_test;

%% FUNDAMENTAL MATRIX
F_7pts_test; test_results.F_7PC = res;
F_3AC_test; test_results.F_3AC = res;



%% GRAPH
figure(); hold on;
methods = fields(test_results);
K = size(fields(test_results),1);
y = zeros(K,1); err = zeros(K,1);
legend_txt = {};
for i = 1:K
    data = test_results.(methods{i});
    inliers = data(1,:) >= data(2,:) & data(1,:) <= data(3,:);
    plot([i i], [data(2,1) data(3,1)], '+'); 
    y(i) = mean(data(1,inliers));
    err(i) = std(data(1,inliers));
    methods{i} = [strrep(methods{i},'_',' (') ')'];
    legend_txt{i} = ['chi2 thresholds for ' methods{i} ' [' sprintf('%.0f%%',round(100*sum(inliers)/num_simulations)) ']'];
end
errorbar(1:K,y,err)
legend_txt{end+1} = 'test statistic (mean, std)';
xticks(1:K); xticklabels(methods); axis([0.5 8.5 0 80]);
title('Statistical test of the uncertatinty propagation (significance level = 0.999)'); 
ylabel('statistics (chi2 threshold, mean, std)');
legend(legend_txt);

