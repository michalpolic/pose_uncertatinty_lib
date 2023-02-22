function [ ResAffines ] = normalizeAffineTransformations( Affines, T1, T2 )
    
    N = size(Affines, 1);    
    ResAffines = zeros(N,4);
        
    for i = 1 : N        
        A = [Affines(i, 1), Affines(i, 2), 0;
            Affines(i, 3), Affines(i, 4), 0;
            0, 0, 1];
        
        A = T2 * A * inv(T1);

        ResAffines(i,:) = [A(1,1), A(1,2), A(2,1), A(2,2)];
    end
end

