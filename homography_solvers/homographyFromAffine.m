%-----------------------------------------------------------------%
% Copyright 2014-2016 , Daniel Barath barath.daniel@sztaki.mta.hu %
%-----------------------------------------------------------------%
% Affines
% 	Type: Matrix
%	Size: N * 4
%	Structure: Each row contains the parameteters of an affine transformation.
%	a11, a12, a21, a22
% pts1, pts2
%	Type: Vector
%	Size: N * 2
%	The points on the first, and second images
function H = homographyFromAffine(correspondences)

    N   = size(correspondences, 1); 
    A   = [];
    
    % Normalize the points
    [T1,offset1,s1,s1_1,s1_2] = getNormalizingTransformation(correspondences(:,1:2));
    [T2,offset2,s2,s2_1,s2_2] = getNormalizingTransformation(correspondences(:,3:4));
    
    pts1Tr = (T1*[correspondences(:,1:2)';ones(1,N)])';
    pts2Tr = (T2*[correspondences(:,3:4)';ones(1,N)])';
    
    % Normalize the affine transformations
    affineTr = normalizeAffineTransformations( correspondences(:,5:8), T1, T2 );
    
    normalized_correspondences = [pts1Tr(:,1:2), pts2Tr(:,1:2), affineTr];

    for i = 1 : 6 : N*6
        currRow		= (i-1) / 6 + 1;
        pt1			= normalized_correspondences(currRow, 1:2);
        pt2			= normalized_correspondences(currRow, 3:4);

        a11         = normalized_correspondences(currRow, 5);
        a12         = normalized_correspondences(currRow, 6);
        a21         = normalized_correspondences(currRow, 7);
        a22         = normalized_correspondences(currRow, 8);

        A = [A; [0, 0, 0, pt1(1), pt1(2), 1 , -pt1(1) * pt2(2), -pt1(2) * pt2(2), -pt2(2)]];
        A = [A; [pt1(1), pt1(2), 1, 0, 0, 0, -pt1(1) * pt2(1), -pt1(2) * pt2(1), -pt2(1)]];
        A = [A; [1, 0, 0, 0, 0, 0, -pt2(1) - a11 * pt1(1), -a11 * pt1(2), -a11]];
        A = [A; [0, 1, 0, 0, 0, 0, -a12 * pt1(1), -pt2(1) - a12 * pt1(2), -a12]];
        A = [A; [0, 0, 0, 1, 0, 0, -pt2(2) - a21 * pt1(1), -a21 * pt1(2), -a21]];
        A = [A; [0, 0, 0, 0, 1, 0, -a22 * pt1(1), -pt2(2) - a22 * pt1(2), -a22]];
    end;
        
    [tmp1, tmp2, V]     = svd(A(1:8,:));
    h                   = V(:, 9);
    H                   = reshape(h, 3, 3)';    
    H                   = inv(T2) * H * T1;
end

