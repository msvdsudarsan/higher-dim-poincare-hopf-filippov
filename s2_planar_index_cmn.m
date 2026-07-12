function s2_planar_index_cmn()
    fprintf('\n=== S2: planar CMN index vs regularization (Thm cmn) ===\n');
    ex = { 'rot+ / id',       [0 -1;1 0],  eye(2);
           'source id/id',    eye(2),      eye(2);
           'node2 / rot',     [2 0;0 2],   [0 -1;1 0];
           'sink / sink',    -eye(2),     -eye(2);
           'saddle / saddle', [1 0;0 -1],  [1 0;0 -1] };
    fprintf('  %-16s %7s %7s %7s %8s %7s\n', 'example','W+','W-','corr','I_index','I_reg');
    for i = 1:size(ex,1)
        [Wp,Wm,corr,Iidx,Ireg] = planar_one(ex{i,2}, ex{i,3});
        fprintf('  %-16s %7.3f %7.3f %7.3f %8.3f %7.3f\n', ex{i,1}, Wp, Wm, corr, Iidx, Ireg);
    end
    fprintf('  (generic examples only; antiparallel F^+ = -lambda F^- on E is excluded)\n');

    plot_portrait([2 0;0 2], [0 -1;1 0], 'fig_s2_planar');
end

function [Wp,Wm,corr,Iidx,Ireg] = planar_one(Ap, Am)
    r = 1; N = 6000;
    th_up = linspace(0, pi, N);
    th_lo = linspace(pi, 2*pi, N);
    Pup = [r*cos(th_up); r*sin(th_up)];  Plo = [r*cos(th_lo); r*sin(th_lo)];
    Fp = Ap*Pup;  Fm = Am*Plo;
    aUp = unwrap(atan2(Fp(2,:), Fp(1,:)));
    aLo = unwrap(atan2(Fm(2,:), Fm(1,:)));
    Wp = (aUp(end)-aUp(1)) / (2*pi);
    Wm = (aLo(end)-aLo(1)) / (2*pi);
    brg = @(u,v,M) (1-linspace(0,1,M)).*u + linspace(0,1,M).*v;
    W = [Fp, brg(Fp(:,end),Fm(:,1),600), Fm, brg(Fm(:,end),Fp(:,1),600)];
    ang = unwrap(atan2(W(2,:), W(1,:)));
    Iidx = (ang(end)-ang(1)) / (2*pi);
    corr = Iidx - Wp - Wm;
    Iidx = round(Iidx*1e6)/1e6;
    eps = 0.02; M = 8*N; th = linspace(0, 2*pi, M); P = [r*cos(th); r*sin(th)];
    ph = sin(pi*max(min(P(2,:)/eps,1),-1)/2);
    Z = (1+ph)/2.*(Ap*P) + (1-ph)/2.*(Am*P);
    az = unwrap(atan2(Z(2,:), Z(1,:)));
    Ireg = (az(end)-az(1)) / (2*pi);
end

function plot_portrait(Ap, Am, name)
    eps = 0.05;
    [X,Y] = meshgrid(linspace(-1.3,1.3,26), linspace(-1.3,1.3,26));
    P = [X(:).'; Y(:).'];
    ph = sin(pi*max(min(P(2,:)/eps,1),-1)/2);
    Z = (1+ph)/2.*(Ap*P) + (1-ph)/2.*(Am*P);
    U = reshape(Z(1,:), size(X)); Vv = reshape(Z(2,:), size(Y));
    f = figure('Color','w');
    quiver(X, Y, U, Vv, 1.2, 'Color', [0.2 0.3 0.8]); hold on;
    th = linspace(0,2*pi,400); plot(cos(th), sin(th), 'k-', 'LineWidth', 1.4);
    plot([-1.3 1.3], [0 0], 'r--', 'LineWidth', 1.2);
    plot([1 -1], [0 0], 'ro', 'MarkerFaceColor', 'r');
    text(1.02, 0.12, 'p_+'); text(-1.15, 0.12, 'p_-');
    text(0, 1.15, '\Gamma^+ (F^+)', 'HorizontalAlignment','center');
    text(0, -1.2, '\Gamma^- (F^-)', 'HorizontalAlignment','center');
    axis equal; axis([-1.35 1.35 -1.35 1.35]); grid on;
    title('Regularized Filippov field, switching line \Sigma, and \partial B');
    exportfig(f, name);
end
