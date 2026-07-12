function W = assembled_field(P, Ap, Am, delta)
    z = P(3,:);
    q = [P(1,:); P(2,:); zeros(1, size(P,2))];
    qn = vecnorm(q); qn(qn < 1e-12) = 1; q = q ./ qn;
    lam = 0.5 + z/(2*delta);
    lam = max(min(lam, 1), 0);
    Wtube = lam.*(Ap*q) + (1-lam).*(Am*q);
    Wp = Ap*P;  Wm = Am*P;
    W = Wtube;
    up = z >  delta;  lo = z < -delta;
    W(:, up) = Wp(:, up);
    W(:, lo) = Wm(:, lo);
end
