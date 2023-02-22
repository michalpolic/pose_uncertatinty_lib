% Unity test for uncertainty propagation from 2AC to Fundamental matrix 
% with assumption of one focal length
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

num_model_params = 9;
num_sample_pts = 2;
% num_simulations = 5;
% num_iner_iterations = 100;
% sigma2 = 10^(-13);
in_cov = sigma2*eye(16);
units_scale = 1;
inliers_threshold = 1 * units_scale;
min_num_inl = 50;

% test data
load('test_data/Ef_correspondences_AC.mat')
correspondences = units_scale * correspondences;
% correspondences = (correspondences - repmat(0.001*[960 540 960 540 0 0 0 0],size(correspondences,1),1));

%% Monte Carlo simulation

fprintf('Run %d unit tests: 2AC -> F (E + f)\n',num_simulations);
i = 1;
while (i <= num_simulations)
    
    % select correspondences
    sample_pts = []; %[33 45];
    while length(unique(sample_pts)) ~= num_sample_pts
        sample_pts = randi(size(correspondences,1),num_sample_pts,1);
    end
    sample_correspondences = correspondences(sample_pts,:);
   
    % estimate the fundamental matrix
    [Fs,fs] = EfMatrixFromAffine(sample_correspondences);
    if isempty(Fs)
        num_simulations = num_simulations + 1;
        continue;
    end
    [F, num_of_inliers, model_id] = verifyCorrectModel( Fs, correspondences, inliers_threshold );
    if num_of_inliers < min_num_inl
        num_simulations = num_simulations + 1;
        continue;
    end
    f = fs(model_id);
    F = normalizeMatrix(F);
%     F = normalizeMatrix(diag([f,f,1]) * F * diag([f,f,1]));   
%     [b, ~, aa] = bR_from_E_uv(diag([f,f,1]) * F * diag([f,f,1]), a2h(correspondences(:,1:2)')',a2h(correspondences(:,3:4)')');
    
    
    % propagate the covariance
    covF = Ef_2AC( sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', ... 
        sample_correspondences(1,5:end), sample_correspondences(2,5:end), F, in_cov );
%     covF = Ef_2AC_br( sample_correspondences(:,1:2)', sample_correspondences(:,3:4)', ... 
%         sample_correspondences(1,5:end), sample_correspondences(2,5:end), b, aa, f, in_cov );

    
    % simulate the covariance propagation by MC   
    distorted_models = zeros(num_iner_iterations,num_model_params);
    models_filter = zeros(num_iner_iterations,1);
    %for j = 1:num_iner_iterations
    j = 1;
    while sum(models_filter) < num_iner_iterations
        % distort points by in_cov
        distorted_correspondences = sample_correspondences + sqrt(sigma2)*randn(size(sample_correspondences));

        % compute distorted model
        [dF_models, dfs] = EfMatrixFromAffine(distorted_correspondences);
        if isempty(dF_models)
           continue; 
        end
%         [dF, d_inl, dmodel_id] = verifyCorrectModel( dF_models, correspondences, inliers_threshold );
%         mod_dist = arrayfun(@(k)norm( ...
%             abs(normalizeMatrix(diag([dfs(k),dfs(k),1]) * dF_models(:,:,k) * diag([dfs(k),dfs(k),1]))) - abs(F)),1:size(dF_models,3));
        mod_dist = arrayfun(@(k)norm(abs(normalizeMatrix(dF_models(:,:,k))) - abs(F)),1:size(dF_models,3));
        [min_dist, dmodel_id] = min(mod_dist);
        if min_dist < 0.01
            models_filter(j) = 1;
        end
        dF = dF_models(:,:,dmodel_id);
        df = dfs(dmodel_id);
        dF = normalizeMatrix(dF);     %diag([df,df,1]) * dF * diag([df,df,1])
%         dF = normalizeMatrix(diag([df,df,1]) * dF * diag([df,df,1]));  
        if sign(dF(1)) ~= sign(F(1))
            dF = -dF;
        end
%         df = dfs(dmodel_id);
%         [db, ~, daa] = bR_from_E_uv(diag([df,df,1]) * dF * diag([df,df,1]), a2h(correspondences(:,1:2)')',a2h(correspondences(:,3:4)')');
        distorted_models(j,:) = dF(:); %[db; daa; df];
        j = j + 1;
    end
%     distorted_models = diag(sign(distorted_models(:,1))) * distorted_models;
%     if sign(distorted_models(1,1)) ~= sign(F(1))
%         distorted_models = - distorted_models;
%     end
    MC_covF = cov(distorted_models(logical(models_filter),:));   % - repmat(F(:)',num_iner_iterations,1) [b;aa;f]'

    % measure the error of uncertatinty propagation
    H = null(covF);
    J = null(H');
    scaleCovF = (J'*MC_covF*J)  ./ (J'*covF*J);
    [L,T,~]=check_CovM( J'*MC_covF*J, (J'*covF*J) * mean(scaleCovF(:)), sum(models_filter), 0.999);
    if L > T(1) && L < T(2)
        fprintf('> No reason to doubt: %.2f in [%.2f,%.2f]\n',L,T(1),T(2));
    else
        fprintf('> Reason to doubt: %.2f not in [%.2f,%.2f] **********\n',L,T(1),T(2)); 
%         display(J'*MC_covF*J)
%         display(J'*covF*J)
    end
    i = i + 1;
end





