function s4_candidate_index_n3()
    fprintf('\n=== S4: candidate index (n=3, assembled) vs regularization ===\n');
    fprintf('  [evidence only; n>=3 integrality is CONJECTURAL, not proved]\n');
    epsr = 0.03; d1 = 0.05; d2 = 0.10;
    ex = { 'id / id (source)', eye(3),          eye(3);
           'id / rotX',        eye(3),          [1 0 0;0 0 -1;0 1 0];
           'node2 / id',       2*eye(3),        eye(3);
           'reflectZ / id',    diag([1 1 -1]),  eye(3) };
    fprintf('  %-18s %11s %11s %12s %11s\n', 'example','cand(d=.05)','cand(d=.10)','cand(d->0)','regularized');
    for i = 1:size(ex,1)
        Ap = ex{i,2}; Am = ex{i,3};
        c1 = deg_sphere(@(P) assembled_field(P, Ap, Am, d1), 320, 320);
        c2 = deg_sphere(@(P) assembled_field(P, Ap, Am, d2), 320, 320);
        c0 = c1 - (c2-c1)/(d2-d1)*d1;
        dR = deg_sphere(@(P) reg_field(P, Ap, Am, epsr), 320, 320);
        fprintf('  %-18s %11.3f %11.3f %12.3f %11.3f\n', ex{i,1}, c1, c2, c0, dR);
    end
    fprintf('  (compare the ''cand(d->0)'' column with ''regularized''.)\n');

    Ap = eye(3); Am = [1 0 0;0 0 -1;0 1 0]; delta = 0.15;
    Nu = 60; Nv = 60; uu = linspace(0,2*pi,Nu); vv = linspace(1e-2,pi-1e-2,Nv);
    [U,V] = meshgrid(uu,vv);
    X = sin(V).*cos(U); Y = sin(V).*sin(U); Zc = cos(V);
    C = double(abs(Zc) <= delta);
    f = figure('Color','w');
    surf(X, Y, Zc, C, 'EdgeColor','none'); axis equal off; view(40,18);
    colormap([0.85 0.85 0.9; 0.95 0.55 0.2]);
    title('Boundary sphere \partial B: equatorial tube N_\delta(E) highlighted');
    exportfig(f, 'fig_s4_tube');
end
