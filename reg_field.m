function W = reg_field(P, Ap, Am, eps)
    z = P(3,:);
    ph = sin(pi * max(min(z/eps, 1), -1) / 2);
    lam = 0.5*(1 + ph);
    W = lam.*(Ap*P) + (1-lam).*(Am*P);
end
