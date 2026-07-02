function W = reg_field(P, Ap, Am, eps)
%REG_FIELD  Sotomayor-Teixeira regularization of Z = (F^+, F^-)_f on S^2,
%   eq. (reg), with f(x,y,z) = z and transition function phi(s)=sin(pi*s/2)
%   clipped to [-1,1].  P : 3xK points; Ap, Am : 3x3 with F^+- = Ap*x, Am*x.
    z = P(3,:);
    ph = sin(pi * max(min(z/eps, 1), -1) / 2);
    lam = 0.5*(1 + ph);
    W = lam.*(Ap*P) + (1-lam).*(Am*P);
end
