% Unity test for uncertainty propagation from 2 AC to Homography matrix 
% Implemented by: Michal Polic, michal.polic(at)cvut.cz

% settings
path_to_eigen = 'd:/School/Research/Reconstruction3D/Software/vcpkg/installed/x64-windows/include/eigen3';


%% process
if exist(path_to_eigen,'dir')

    mex(['-I' '"' path_to_eigen '"'],'-g','-DBUILD_MEX', 'solveE_nister_LAF.cc','Ematrix_5pt.cc',...
        'polyquotient.cc', 'polydet.cc', 'sturm.cc', 'Ematrix_6pt.cc')

    mex(['-I' '"' path_to_eigen '"'],'-g','-DBUILD_MEX', 'solveE_nister.cc','Ematrix_5pt.cc',...
        'polyquotient.cc', 'polydet.cc', 'sturm.cc', 'Ematrix_6pt.cc')

    mex(['-I' '"' path_to_eigen '"'],'-g','-DBUILD_MEX', 'solveEf_affine.cc','Efmatrix_5pt.cc',...
        'polyquotient.cc', 'polydet.cc', 'sturm.cc', 'Efmatrix_6pt.cc')

    mex(['-I' '"' path_to_eigen '"'],'-g','-DBUILD_MEX', 'partsolveEf6affine.cc','Efmatrix_6pt.cc',...
        'polyquotient.cc', 'polydet.cc', 'sturm.cc')
    
else
   error('Please select the path to eigen library in ./relative_pose_solvers/compile_mex.m');
end