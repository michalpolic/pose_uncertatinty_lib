function [ E ] = essentialMatrixFromPoints( corresp )
%ESSENTIALFROMPOINTS - wrapper that run the modified Nister cpp solver
    E = solveE_nister(a2h(corresp(:,3:4)'), a2h(corresp(:,1:2)'));
end