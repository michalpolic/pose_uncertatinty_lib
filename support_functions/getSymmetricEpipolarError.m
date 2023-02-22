function residual = getSymmetricEpipolarError(F, correspondence)

    pt1 = [correspondence(1:2)'; 1];
    pt2 = [correspondence(3:4)'; 1];
    
    l1 = F' * pt2;
    l2 = F * pt1;

    d1 = abs(pt1' * l1 / norm(l1(1:2)));
    d2 = abs(pt2' * l2 / norm(l2(1:2)));
    residual = 0.5 * (d1 + d2);
end

