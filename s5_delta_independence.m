function s5_delta_independence()
%S5  delta-dependence of the candidate index and its delta -> 0 limit
%   (evidence for the conjecture that the candidate index is well defined via
%   the delta -> 0 limit; paper Conjectures 'limit'/'deltaindep').  Evidence only.
%
%   At finite delta the assembled degree is NOT exactly an integer: it carries
%   an O(delta) bias. What the paper conjectures is that the delta -> 0 LIMIT is
%   the correct integer. The linear fit below estimates that limit.
    fprintf('\n=== S5: candidate index vs delta and its delta->0 limit ===\n');
    fprintf('  [evidence only; NOT a proof]\n');
    Ap = eye(3); Am = [1 0 0;0 0 -1;0 1 0];
    deltas = [0.05 0.08 0.10 0.15 0.20 0.30];
    vals = zeros(size(deltas));
    for i = 1:numel(deltas)
        G = @(P) assembled_field(P, Ap, Am, deltas(i));
        vals(i) = deg_sphere(G, 320, 320);
        fprintf('  delta = %.2f   candidate index = %+.4f\n', deltas(i), vals(i));
    end
    p = polyfit(deltas, vals, 1);   % candidate ~ p(2) + p(1)*delta
    fprintf('  linear fit: candidate ~ %+.3f  %+.3f*delta\n', p(2), p(1));
    fprintf('  ==> delta->0 estimate = %+.3f   (compare with integer degree +1)\n', p(2));

    f = figure('Color','w');
    plot(deltas, vals, 'o', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 7); hold on;
    xf = [0 max(deltas)]; plot(xf, polyval(p, xf), 'r--', 'LineWidth', 1.3);
    plot(0, p(2), 'rp', 'MarkerSize', 13, 'MarkerFaceColor', 'r');
    grid on; xlabel('\delta  (tube radius)'); ylabel('candidate index');
    ylim([min([vals p(2)])-0.3, max([vals p(2)])+0.3]);
    legend('computed', 'linear fit', '\delta\rightarrow0 estimate', 'Location','best');
    title('Candidate index vs \delta and its \delta\rightarrow0 limit');
    exportfig(f, 'fig_s5_delta');
end
