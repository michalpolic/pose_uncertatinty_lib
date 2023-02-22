function [ E, num_of_inliers, best_model_id ] = verifyCorrectModel( Es, correspondences, theshold, fs )
    num_of_inliers = 0;
    best_model_id = 1;
    p1 = a2h(correspondences(:,1:2)');
    p2 = a2h(correspondences(:,3:4)');
    for i = 1:size(Es,3)
        if nargin > 3
           p1(3,:) = fs(i); 
           p2(3,:) = fs(i); 
        end
        err = sum(p2 .* (Es(:,:,i) * p1)).^2;
        inliers = err < theshold;
        %sum(inliers)
        if sum(inliers) > num_of_inliers
            num_of_inliers = sum(inliers);
            E = Es(:,:,i);
            best_model_id = i;
        end
    end
end

