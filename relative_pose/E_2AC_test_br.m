% Unity test for uncertainty propagation from 2AC to Essential matrix 
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% settings
model_params = 6;
num_sample_pts = 2;
% num_simulations = 5;
% num_iner_iterations = 100;
% sigma2 = 10^(-13);
in_cov = diag(sigma2*ones(1,16));
inliers_threshold = 3 / 1500;

% test data
load('test_data/relative_pose_correspondences.mat');


%% Monte Carlo simulation
res = zeros(3, num_simulations);
fprintf('Run %d unit tests: 2AC -> E\n',num_simulations);
for i = 1:num_simulations
    
    % select correspondences
    sample_pts = [24 44];  
    while length(unique(sample_pts)) ~= num_sample_pts
        sample_pts = randi(size(correspondences,1),num_sample_pts,1);
    end
    sample_correspondences = correspondences(sample_pts,:);
   
    % estimate the essential matrix
    E_models = essentialMatrixFromAffine(sample_correspondences);
    [E,num_of_inliers] = verifyCorrectModel(E_models, correspondences, inliers_threshold );
    [b, ~, aa] = bR_from_E_uv(E, a2h(correspondences(:,1:2)')',a2h(correspondences(:,3:4)')');
                                    
    % uncertatinty propagation
    covE = E_2AC_br( sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', ...
        sample_correspondences(1,5:8), sample_correspondences(2,5:8), b, aa, in_cov);

    % simulate the covariance propagation by MC   
    distorted_models = zeros(num_iner_iterations,model_params);
    for j = 1:num_iner_iterations
        % distort points by in_cov
        distorted_correspondences = sample_correspondences + sqrt(sigma2)*randn(2,8);

        % compute distorted model
        distorted_E_models = essentialMatrixFromAffine(distorted_correspondences);
        dE = verifyCorrectModel(distorted_E_models, correspondences, inliers_threshold );
        [db, ~, daa] = bR_from_E_uv(dE, a2h(correspondences(:,1:2)')',a2h(correspondences(:,3:4)')');
        distorted_models(j,:) = [db; daa]';
    end
    MC_covE = cov(distorted_models - repmat([b;aa]',num_iner_iterations,1));

    
    % measure the error of uncertatinty propagation
    H = null(covE);
    J = null(H');
    [L,T,~]=check_CovM( J'*MC_covE*J, J'*covE*J, num_iner_iterations, 0.999);
    if L > T(1) && L < T(2)
        fprintf('> No reason to doubt: %.2f in [%.2f,%.2f]\n',L,T(1),T(2));
    else
        fprintf('> Reason to doubt: %.2f not in [%.2f,%.2f] **********\n',L,T(1),T(2)); 
    end
    res(:,i) = [L; T(1); T(2)];
end




