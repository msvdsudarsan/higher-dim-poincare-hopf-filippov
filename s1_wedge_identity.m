function s1_wedge_identity()
    fprintf('\n=== S1: wedge-product magnitude identity (eq. wedge) ===\n');
    rng(1);
    dims = [2 3 4 5];  N = 20000;  maxerr = 0;
    for n = dims
        e = 0;
        for k = 1:N
            u = randn(n,1); v = randn(n,1);
            lhs = wedge_norm_sq(u,v) + (u.'*v)^2;
            rhs = (u.'*u)*(v.'*v);
            e = max(e, abs(lhs-rhs));
        end
        fprintf('  n=%d : max|LHS-RHS| = %.3e\n', n, e);
        maxerr = max(maxerr, e);
    end
    fprintf('  OVERALL max error = %.3e  (expected ~1e-12)\n', maxerr);

    n = 3; M = 500; wv = zeros(M,1); pred = zeros(M,1);
    for k = 1:M
        u = randn(n,1); v = randn(n,1);
        wv(k) = sqrt(wedge_norm_sq(u,v));
        th = acos( (u.'*v)/(norm(u)*norm(v)) );
        pred(k) = norm(u)*norm(v)*sin(th);
    end
    f = figure('Color','w');
    plot(pred, wv, '.', 'MarkerSize', 8); hold on;
    mx = max(pred);
    plot([0 mx], [0 mx], 'r-', 'LineWidth', 1.3);
    xlabel('norm(u) norm(v) sin(theta)'); ylabel('norm(u wedge v)');
    title('Wedge magnitude identity in R^3'); axis equal tight; grid on;
    legend('samples','y = x','Location','SouthEast');
    exportfig(f, 'fig_s1_wedge');
end

function s = wedge_norm_sq(u, v)
    n = numel(u); s = 0;
    for i = 1:n-1
        for j = i+1:n
            s = s + (u(i)*v(j) - u(j)*v(i))^2;
        end
    end
end
