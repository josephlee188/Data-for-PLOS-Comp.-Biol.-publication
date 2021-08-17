function idx = findMemberIdx(xy, unit, n, BoxSize, scanR, threshAngle)
    dr = xy - xy(:, n);
    dr = sqrt(dr(1, :).^2 + dr(2, :).^2);
    dr(1, n) = nan;
    dR = xy - (xy(:, n) - BoxSize*(xy(:, n) > (BoxSize-scanR)) + ...
        BoxSize*(xy(:, n) < (scanR)));
    dR = sqrt(dR(1,:).^2 + dR(2,:).^2);
    dR(1, n) = nan;
    
    A = squeeze(dr <= scanR); 
    B = squeeze(dR <= scanR);
    C = A | B;
    
    c = unit(1, :) .* unit(1, n) + unit(2, :) .* unit(2, n);
    c(1, n) = nan;
    c = c >= threshAngle;
    
    D = c & C;
    
    idx = find(D);

end