function s3_gauss_degree_S2()
%S3  Gauss-map degree normalization on S^2 (eq. W), with vol(S^2)=4*pi and the
%   orientation fixed so a degree-d map contributes exactly d.  This validates
%   the normalizing constant used throughout Sec. 'Construction of the
%   candidate local index'.
    fprintf('\n=== S3: Gauss-map degree normalization on S^2 (vol = 4*pi) ===\n');
    d_id   = deg_sphere(@(P) P,      400, 200);
    d_anti = deg_sphere(@(P) -P,     400, 200);
    d_sq   = deg_sphere(@sqmap,      500, 250);
    fprintf('  identity  map : degree = %+.3f   (exact  +1)\n', d_id);
    fprintf('  antipodal map : degree = %+.3f   (exact  -1 on S^2)\n', d_anti);
    fprintf('  z -> z^2  map : degree = %+.3f   (exact  +2)\n', d_sq);

    % Figure: local signed area-density (Jacobian) of the degree-2 map
    Nu = 240; Nv = 120;
    uu = linspace(0,2*pi,Nu); vv = linspace(1e-3,pi-1e-3,Nv);
    [U,V] = meshgrid(uu,vv);
    X = sin(V).*cos(U); Y = sin(V).*sin(U); Zc = cos(V);
    f = figure('Color','w');
    surf(X, Y, Zc, cos(2*U), 'EdgeColor','none'); axis equal off;
    title('Degree-2 Gauss map z \rightarrow z^2 on S^2'); colormap(parula); view(35,20);
    exportfig(f, 'fig_s3_degree');
end

function W = sqmap(P)
%SQMAP  Degree-2 map of S^2 induced by z -> z^2 under stereographic id.
    x = P(1,:); y = P(2,:); z = P(3,:);
    w = (x + 1i*y) ./ (1 - z + 1e-12);
    w2 = w.^2;
    den = abs(w2).^2 + 1;
    W = [2*real(w2)./den; 2*imag(w2)./den; (abs(w2).^2 - 1)./den];
end
