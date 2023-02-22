% Unity test for uncertainty propagation from 5 poitns to Essential matrix 
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% settings
num_model_params = 6;
num_sample_pts = 5;
% num_simulations = 5;
% num_iner_iterations = 100;
% sigma2 = 10^(-13);
in_cov = sigma2 * eye(20);
inliers_threshold = 3 / 1500;

% test data
load('test_data/scene_dataset1.mat');
% ShowCameras(P, diag([1 1 1]), m, M, true, false, true, 1:size(m{1},2), m);

    
%% Monte Carlo simulation
res = zeros(3, num_simulations);
% select tested camera pair / correspondences  [] = random
camera_pair = [];   %[2 4];     
sample_pts = [];    %[33 57 69 38 14];  

fprintf('Run %d unit tests: 5pts -> E\n',num_simulations);
for i = 1:num_simulations

    % select correspondences
    while length(unique(camera_pair)) ~= 2
        camera_pair = randi(size(P,2),1,2);
    end
    [gtE, gtb, gtR, gtaa] = essentialMatrix4TwoCameras(P{camera_pair(1)}, P{camera_pair(2)});
    correspondences = [h2a(P{camera_pair(1)}*a2h(M))' h2a(P{camera_pair(2)}*a2h(M))'];
    corresp_in_image = correspondences < 1 & correspondences > -1;
    while length(unique(sample_pts)) ~= num_sample_pts
        sample_pts = randi(size(correspondences,1),num_sample_pts,1);
        if sum(sum(corresp_in_image(sample_pts,:))) ~= 4*num_sample_pts   % check if all the sample points are in images 
            sample_pts = [];
        end
    end
    sample_correspondences = correspondences(sample_pts,:);  


    % % show selected correspondences
    % fprintf('Cameras: %d - %d\n',camera_pair(1),camera_pair(2));
    % mp = cell(1,2); mp{1} = sample_correspondences(:,1:2)'; mp{2} = sample_correspondences(:,3:4)';
    % ShowCameras(P(camera_pair), diag([1 1 1]), mp, M, true, false, true, 1:size(mp{1},2), mp);


    % estimate the essential matrix
    E_models = essentialMatrixFromPoints(sample_correspondences);
    E = verifyCorrectModel(E_models, correspondences, inliers_threshold);

    % estimated model
    [b, ~, aa] = bR_from_E_uv(E, a2h(sample_correspondences(:,1:2)')', a2h(sample_correspondences(:,3:4)')');

    % uncertatinty propagation
    covE = E_5pts_br(sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', E, in_cov);

    % simulate the covariance propagation by MC   
    distorted_models = zeros(num_iner_iterations,num_model_params);
    min_model_distances = zeros(num_iner_iterations,1);
    for j = 1:num_iner_iterations
        % distort points by in_cov
        distorted_correspondences = sample_correspondences + sqrt(sigma2)*randn(size(sample_correspondences));

        % compute distorted model
        distorted_E_models = essentialMatrixFromPoints(distorted_correspondences);
        distorted_E = verifyCorrectModel(distorted_E_models, correspondences, inliers_threshold);
        [d_b, ~, d_aa] = bR_from_E_uv(distorted_E,  a2h(distorted_correspondences(:,1:2)')', ...
                        a2h(distorted_correspondences(:,3:4)')');         

        % store the change
        distorted_models(j,:) = [d_b; d_aa]'; 
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
