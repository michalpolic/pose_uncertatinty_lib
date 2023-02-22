function [ Es, fs ] = EfMatrixFromAffine( corresp )
    [Es, fs] = solveEf_affine(a2h(corresp(:,1:2)'), a2h(corresp(:,3:4)'), corresp(:,5:8)');
    for i = 1:size(Es,3)
%         Ki = [1/fs(i) 0 0; 0 1/fs(i) 0; 0 0 1];
%         Es(:,:,i) = Ki' * Es(:,:,i)' * Ki;
        Es(:,:,i) = Es(:,:,i)';
    end
end

