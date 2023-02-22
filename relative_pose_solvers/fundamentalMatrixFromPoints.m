function Fmodels = fundamentalMatrixFromPoints(correspondences)
    n = size(correspondences,1);
    A = [];
    Fs = {};
    
    % Normalize the points
    [T1,offset1,s1,s1_1,s1_2] = getNormalizingTransformation(correspondences(:,1:2));
    [T2,offset2,s2,s2_1,s2_2] = getNormalizingTransformation(correspondences(:,3:4));
    
    pts1Tr = (T1*[correspondences(:,1:2)';ones(1,n)])';
    pts2Tr = (T2*[correspondences(:,3:4)';ones(1,n)])';
    
    normalized_correspondences = [pts1Tr(:,1:2), pts2Tr(:,1:2)];
    
    for i = 1 : n
        x0 = normalized_correspondences(i,1);
        y0 = normalized_correspondences(i,2);

        x1 = normalized_correspondences(i,3);
        y1 = normalized_correspondences(i,4);

        A(i, 1) = x1 * x0;
        A(i, 2) = x1 * y0;
        A(i, 3) = x1;
        A(i, 4) = y1 * x0;
        A(i, 5) = y1 * y0;
        A(i, 6) = y1;
        A(i, 7) = x0;
        A(i, 8) = y0;
        A(i, 9) = 1;
    end

    %% Apply the 7-point algorithm when a minimal sample is given
    if n == 7
        [tmp1,tmp2,V] = svd(A);
        f1 = V(:,9);
        f2 = V(:,8);

        % f1, f2 is a basis => lambda*f1 + mu*f2 is an arbitrary f. matrix.
        % as it is determined up to a scale, normalize lambda & mu (lambda + mu = 1),
        % so f ~ lambda*f1 + (1 - lambda)*f2.
        % use the additional constraint det(f) = det(lambda*f1 + (1-lambda)*f2) to find lambda.
        % it will be a cubic equation.
        % find c - polynomial coefficients.
        for i = 1 : 9
           f1(i) = f1(i) - f2(i);
        end

        c = [];
        t0 = f2(5) * f2(9) - f2(6) * f2(8);
        t1 = f2(4) * f2(9) - f2(6) * f2(7);
        t2 = f2(4) * f2(8) - f2(5) * f2(7);

        c(4) = f2(1) * t0 - f2(2) * t1 + f2(3) * t2;

        c(3) = f1(1) * t0 - f1(2) * t1 + f1(3) * t2 - ...
            f1(4) * (f2(2) * f2(9) - f2(3) * f2(8)) + ...
            f1(5) * (f2(1) * f2(9) - f2(3) * f2(7)) - ...
            f1(6) * (f2(1) * f2(8) - f2(2) * f2(7)) + ...
            f1(7) * (f2(2) * f2(6) - f2(3) * f2(5)) - ...
            f1(8) * (f2(1) * f2(6) - f2(3) * f2(4)) + ...
            f1(9) * (f2(1) * f2(5) - f2(2) * f2(4));

        t0 = f1(5) * f1(9) - f1(6) * f1(8);
        t1 = f1(4) * f1(9) - f1(6) * f1(7);
        t2 = f1(4) * f1(8) - f1(5) * f1(7);

        c(2) = f2(1) * t0 - f2(2) * t1 + f2(3) * t2 - ...
            f2(4) * (f1(2) * f1(9) - f1(3) * f1(8)) + ...
            f2(5) * (f1(1) * f1(9) - f1(3) * f1(7)) - ...
            f2(6) * (f1(1) * f1(8) - f1(2) * f1(7)) + ...
            f2(7) * (f1(2) * f1(6) - f1(3) * f1(5)) - ...
            f2(8) * (f1(1) * f1(6) - f1(3) * f1(4)) + ...
            f2(9) * (f1(1) * f1(5) - f1(2) * f1(4));

        c(1) = f1(1) * t0 - f1(2) * t1 + f1(3) * t2;

        r = roots(c);

        for ri = 1 : length(r)

            if abs(imag(r(ri))) > 0
                continue;
            end

            % for each root form the fundamental matrix
            root = real(r(ri));
            lambda = root;
            mu = 1.0;
            s = f1(9) * root + f2(9);

            % normalize each matrix, so that F(3,3) (~fmatrix[8]) == 1
            if abs(s) > 1e-5

                mu = 1.0 / s;
                lambda = lambda * mu;

                for i = 1 : 8
                    f(i) = f1(i) * lambda + f2(i) * mu;
                end
                f(9) = 1.0;

                F = reshape(f, 3, 3)';
                Fs{end + 1} = T2' * F * T1;
           end
        end 
    else 
        %% Apply the 8-point algorithm otherwise
        [tmp1,tmp2,V] = svd(A);
        f = V(:,9);
        
        F = reshape(f, 3, 3)';
        Fs{1} = T2' * F * T1;
    end
    
    % rewrite Fs into unified format
    Fmodels = zeros(3,3,length(Fs));
    for i = 1:length(Fs)
        Fmodels(:,:,i) = Fs{i};
    end
end
