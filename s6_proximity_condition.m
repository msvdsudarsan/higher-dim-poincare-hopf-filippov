function s6_proximity_condition()
    fprintf('\n=== S6: proximity condition eq.(proximity) as sufficient nonvanishing ===\n');
    a = [1;0.2;0]; b = [0.2;1;0]; delta = 0.2; epsr = 0.05;
    ts = linspace(-delta, delta, 400);
    hfun = @(t) (0.5 + t/(2*delta)).*a + (0.5 - t/(2*delta)).*b;
    mu = inf; for t = ts, mu = min(mu, norm(hfun(t))); end
    ap0 = [0;-1;0]; bp0 = [1;0;0];
    scales = linspace(0, 4, 17);
    LHS = zeros(size(scales)); Hmin = zeros(size(scales));
    phi = @(t) sin(pi*max(min(t/epsr,1),-1)/2);
    for k = 1:numel(scales)
        ap = scales(k)*ap0; bp = scales(k)*bp0;
        Zf = @(t) (0.5*(1+phi(t))).*(a + t*ap) + (1 - 0.5*(1+phi(t))).*(b + t*bp);
        LHS(k) = norm(a-b) + delta*(norm(ap)+norm(bp));
        hm = inf;
        for s = linspace(0,1,40)
            for t = ts
                hm = min(hm, norm((1-s)*hfun(t) + s*Zf(t)));
            end
        end
        Hmin(k) = hm;
    end
    fprintf('  mu = min||h|| = %.3f\n', mu);
    for k = 1:numel(scales)
        tag = 'fail'; if LHS(k) < mu, tag = 'HOLD'; end
        fprintf('  scale=%.2f  LHS=%.3f  min||H||=%.3f  [%s]\n', scales(k), LHS(k), Hmin(k), tag);
    end
    f = figure('Color','w');
    plot(LHS, Hmin, 'o-', 'LineWidth', 1.4); hold on;
    yl = ylim; plot([mu mu], [0 max(Hmin)*1.15], 'r--', 'LineWidth', 1.3);
    xlabel('||a-b|| + \delta(||a''|| + ||b''||)   (LHS of proximity)');
    ylabel('min_{s,t} ||H(s,t)||'); grid on;
    title('Proximity condition: LHS < \mu guarantees nonvanishing homotopy');
    legend('min||H||', '\mu = min||h||', 'Location','best');
    exportfig(f, 'fig_s6_proximity');
end
