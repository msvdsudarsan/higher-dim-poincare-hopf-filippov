function d = deg_sphere(gfun, Nu, Nv)
    if nargin < 2, Nu = 400; end
    if nargin < 3, Nv = 200; end
    uu = linspace(0, 2*pi, Nu+1); uu = uu(1:end-1);
    vv = linspace(1e-3, pi-1e-3, Nv);
    du = uu(2)-uu(1); dv = vv(2)-vv(1);
    [U, V] = meshgrid(uu, vv);
    h = 1e-6;
    G  = local_ghat(gfun, U,     V);
    Gu = (local_ghat(gfun, U+h,  V) - G) / h;
    Gv = (local_ghat(gfun, U,   V+h) - G) / h;
    cx = Gv(:,:,2).*Gu(:,:,3) - Gv(:,:,3).*Gu(:,:,2);
    cy = Gv(:,:,3).*Gu(:,:,1) - Gv(:,:,1).*Gu(:,:,3);
    cz = Gv(:,:,1).*Gu(:,:,2) - Gv(:,:,2).*Gu(:,:,1);
    jac = G(:,:,1).*cx + G(:,:,2).*cy + G(:,:,3).*cz;
    d = sum(jac(:)) * du * dv / (4*pi);
end

function Gn = local_ghat(gfun, U, V)
    [m, n] = size(U);
    X = sin(V).*cos(U);  Y = sin(V).*sin(U);  Zc = cos(V);
    pts = [X(:).'; Y(:).'; Zc(:).'];
    W = gfun(pts);
    W = W ./ vecnorm(W);
    Gn = zeros(m, n, 3);
    Gn(:,:,1) = reshape(W(1,:), m, n);
    Gn(:,:,2) = reshape(W(2,:), m, n);
    Gn(:,:,3) = reshape(W(3,:), m, n);
end
