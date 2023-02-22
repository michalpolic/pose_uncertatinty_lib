function E = essentialMatrixFromAffine( corresp )
%ESSENTIALFROMAFFINE - wrapper that run the modified Nister cpp solver
    E = solveE_nister_LAF(corresp(:,3:4)', corresp(:,1:2)', ...
        reshape(corresp(2,5:8),2,2), reshape(corresp(1,5:8),2,2));
end

