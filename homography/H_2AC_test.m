% Unity test for uncertainty propagation from 2 AC to Homography matrix 
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% settings
num_model_params = 9;
num_sample_pts = 2;
% num_simulations = 5;
% num_iner_iterations = 100;
% sigma2 = 10^(-13);
in_cov = sigma2*eye(16);

% propagation type: [] - overconstrained, [1,2,...] - select 8 equations
selected_rows = [1:8 13];

% test data
load('test_data/homography_correspondences_AC.mat');


%% Monte Carlo simulation
res = zeros(3, num_simulations);
fprintf('Run %d unit tests: 2AC -> H\n',num_simulations);
for i = 1:num_simulations
    
    % select correspondences
    sample_pts = [];
    while length(unique(sample_pts)) ~= num_sample_pts
        sample_pts = randi(size(correspondences,1),num_sample_pts,1);
    end
    sample_correspondences = correspondences(sample_pts,:);
   
    % propagate the covariance using H_4pts
    H = normalizeMatrix(homographyFromAffine(sample_correspondences));
    covH = H_2AC( sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', ...
        sample_correspondences(1,5:8), sample_correspondences(2,5:8), H, in_cov , ...
        selected_rows);
    
    % simulate the covariance propagation by MC   
    distorted_models = zeros(num_iner_iterations,9);
    for j = 1:num_iner_iterations
        % distort points by in_cov
        distorted_correspondences = sample_correspondences + sqrt(sigma2)*randn(size(sample_correspondences));
        
        % compute distorted model
        distorted_H = normalizeMatrix(homographyFromAffine(distorted_correspondences));
        distorted_models(j,:) = distorted_H(:)';
    end
    MC_covH = cov(distorted_models - repmat(H(:)',num_iner_iterations,1));
    
    % measure the error of uncertatinty propagation
    H = null(covH);
    J = null(H');
    [L,T,~]=check_CovM( J'*MC_covH*J, J'*covH*J, num_iner_iterations, 0.999);
    if L > T(1) && L < T(2)
        fprintf('> No reason to doubt: %.2f in [%.2f,%.2f]\n',L,T(1),T(2));
    else
        fprintf('> Reason to doubt: %.2f not in [%.2f,%.2f] **********\n',L,T(1),T(2)); 
    end
    res(:,i) = [L; T(1); T(2)];
end





