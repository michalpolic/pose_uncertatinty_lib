function [Fmodels, fs] = EfMatrixFromPoints(correspondences)
    
    [Fs, fs] = solve_f_E_f_elimIdeal(a2h(correspondences(:,1:2)')',...
                                     a2h(correspondences(:,3:4)')',eye(3),eye(3),2);
    
    % rewrite Fs into unified format
    Fmodels = zeros(3,3,length(Fs));
    for i = 1:length(Fs)
        Fmodels(:,:,i) = Fs{i};
    end
end
