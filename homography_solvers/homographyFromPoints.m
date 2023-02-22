function H = homographyFromPoints(correspondences)
    n = size(correspondences,1);
    A = [];
    
    % Normalize the points
    [T1,offset1,s1,s1_1,s1_2] = getNormalizingTransformation(correspondences(:,1:2));
    [T2,offset2,s2,s2_1,s2_2] = getNormalizingTransformation(correspondences(:,3:4));
    
    pts1Tr = (T1*[correspondences(:,1:2)';ones(1,n)])';
    pts2Tr = (T2*[correspondences(:,3:4)';ones(1,n)])';
    
    normalized_correspondences = [pts1Tr(:,1:2), pts2Tr(:,1:2)];
    
    for i = 1 : n
        x1 = normalized_correspondences(i,1);
        y1 = normalized_correspondences(i,2);

        x2 = normalized_correspondences(i,3);
        y2 = normalized_correspondences(i,4);

        A(2*i-1,:) = [-x1,-y1,-1.0,0,0,0,x2*x1,x2*y1,x2];
        A(2*i,:) = [0,0,0,-x1,-y1,-1.0,y2*x1,y2*y1,y2];
    end

    [tmp1,tmp2,V] = svd(A);
    h = V(:,9);
    H = reshape(h, 3, 3)';
    
    H = inv(T2) * H * T1;
end
