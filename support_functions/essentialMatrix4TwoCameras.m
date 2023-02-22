function [ E, b, R, aa ] = essentialMatrix4TwoCameras( P1, P2 )
    R1 = P1(:,1:3);
    R2 = P2(:,1:3);
    R = R2 * R1';
    b = - P2(:,4) + R * P1(:,4);
    E = v2X(b) * R;
    axis = X2v(R - R');
    if sqrt(axis' * axis) < 100 * eps
        axis = null(R - eye(3));
    end
    axis = axis(:,1) * (1/sqrt(axis(:,1)'*axis(:,1)));
    aa = acos(0.5*(trace(R) - 1)) * axis;
end

