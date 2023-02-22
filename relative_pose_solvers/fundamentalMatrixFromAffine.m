function Fmodels = fundamentalMatrixFromAffine(correspondences)
    n = size(correspondences,1);
    A = [];
    Fs = {};
    row_idx = 1;
    eqNum = 3;
        
    for i = 1 : n
        u1          = correspondences(i,1);
        v1          = correspondences(i,2);

        u2          = correspondences(i,3);
        v2          = correspondences(i,4);
        
        a1         = correspondences(i, 5);
        a2         = correspondences(i, 6);
        a3         = correspondences(i, 7);
        a4         = correspondences(i, 8);
        
        from       = eqNum * (i - 1) + 1;
        to         = eqNum * i;
        
        A(from : to, :) = [u1 * u2, v1 * u2, u2, u1 * v2, v1 * v2, v2, u1, v1, 1; ...
            u2 + a1 * u1, a1 * v1, a1, v2 + a3 * u1, a3 * v1, a3, 1, 0, 0; ...
            a2 * u1, u2 + a2 * v1, a2, a4 * u1, v2 + a4 * v1, a4, 0, 1, 0];
    end

    % select the minimal num. of constrains
    A = A(1:end-2,:);
    
    %% Apply the 3-affine algorithm when a minimal sample is given
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
            Fs{end + 1} = F;
       end
    end 
    
    % rewrite Fs into unified format
    Fmodels = zeros(3,3,length(Fs));
    for i = 1:length(Fs)
        Fmodels(:,:,i) = Fs{i};
    end
end
