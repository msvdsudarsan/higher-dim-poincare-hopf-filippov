function W = assembled_field(P, Ap, Am, delta)
%ASSEMBLED_FIELD  The assembled boundary map used to define the candidate
%   local index in dimension n = 3 (paper Sec. 'Construction of the candidate
%   local index').  On the two polar caps of S^2 (|z| > delta) it is the Gauss
%   map of F^+ / F^-; on the equatorial tube (|z| <= delta) it is the convex
%   interpolation field h of eq. (interp), with F^+ , F^- evaluated at the
%   PROJECTED point pi(x) on the equator E (this is exactly the paper's h,
%   which freezes F^+ , F^- along E).  f(x,y,z) = z, so the fiber coordinate
%   is t = z.
%
%   P : 3xK points on S^2.   Ap, Am : 3x3 matrices with F^+- = Ap*x , Am*x.
    z = P(3,:);
    % projected equator point pi(x) = normalize (x, y, 0)
    q = [P(1,:); P(2,:); zeros(1, size(P,2))];
    qn = vecnorm(q); qn(qn < 1e-12) = 1; q = q ./ qn;
    lam = 0.5 + z/(2*delta);          % weight of F^+ , eq. (interp)
    lam = max(min(lam, 1), 0);
    Wtube = lam.*(Ap*q) + (1-lam).*(Am*q);
    Wp = Ap*P;  Wm = Am*P;            % Gauss maps on the caps (at x itself)
    W = Wtube;
    up = z >  delta;  lo = z < -delta;
    W(:, up) = Wp(:, up);
    W(:, lo) = Wm(:, lo);
end
