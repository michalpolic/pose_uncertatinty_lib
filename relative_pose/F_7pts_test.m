% Unity test for uncertainty propagation from 7 poitns to Fundamental matrix 
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% settings
num_model_params = 9;
num_sample_pts = 7;
% num_simulations = 5;
% num_iner_iterations = 100;
% sigma2 = 10^(-13);
in_cov = sigma2*eye(28);
inliers_threshold = 3 / 1000;

% test data
load('test_data/fundamental_matrix_correspondences.mat')


%% Monte Carlo simulation
res = zeros(3, num_simulations);
fprintf('Run %d unit tests: 7pts -> F\n',num_simulations);
for i = 1:num_simulations
    
    % select correspondences
    sample_pts = [];
    while length(unique(sample_pts)) ~= num_sample_pts
        sample_pts = randi(size(correspondences,1),num_sample_pts,1);
    end
    sample_correspondences = correspondences(sample_pts,:);
   
    % estimate the fundamental matrix
    F_models = fundamentalMatrixFromPoints(sample_correspondences);
    [F, num_of_inliers] = verifyCorrectModel( F_models, correspondences, inliers_threshold );
    F = normalizeMatrix(F);
    
    % propagate the covariance
    covF = F_7pts( sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', F, in_cov );

    % simulate the covariance propagation by MC   
    distorted_models = zeros(num_iner_iterations,num_model_params);
    for j = 1:num_iner_iterations
        % distort points by in_cov
        distorted_correspondences = sample_correspondences + sqrt(sigma2)*randn(size(sample_correspondences));

        % compute distorted model
        distorted_F_models = fundamentalMatrixFromPoints(distorted_correspondences);
        distorted_F = verifyCorrectModel( distorted_F_models, correspondences, inliers_threshold );
        dF = normalizeMatrix(distorted_F);
        distorted_models(j,:) = dF(:)';
    end
    MC_covF = cov(distorted_models - repmat(F(:)',num_iner_iterations,1));


    % measure the error of uncertatinty propagation
    H = null(covF);
    J = null(H');
    [L,T,~]=check_CovM( J'*MC_covF*J, J'*covF*J, num_iner_iterations, 0.999);
    if L > T(1) && L < T(2)
        fprintf('> No reason to doubt: %.2f in [%.2f,%.2f]\n',L,T(1),T(2));
    else
        fprintf('> Reason to doubt: %.2f not in [%.2f,%.2f] **********\n',L,T(1),T(2)); 
    end
    res(:,i) = [L; T(1); T(2)];
end


