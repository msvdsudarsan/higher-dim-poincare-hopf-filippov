function exportfig(f, name)
%EXPORTFIG  Save figure f as <name>.pdf and .png to the first writable folder.
%   Never throws: if no folder is writable, it just prints a note and returns,
%   so the calling script always finishes.
    thisdir = fileparts(mfilename('fullpath'));
    cands = { fullfile(thisdir, '..', 'figures'), ...
              fullfile(pwd, 'figures'), ...
              pwd, ...
              tempdir };
    outdir = '';
    for i = 1:numel(cands)
        d = cands{i};
        try
            if ~exist(d, 'dir'); mkdir(d); end
            % probe writability
            probe = fullfile(d, '.probe_write');
            fid = fopen(probe, 'w');
            if fid ~= -1
                fclose(fid); delete(probe);
                outdir = d; break;
            end
        catch
            % try next candidate
        end
    end
    if isempty(outdir); outdir = pwd; end
    saved = false;
    try
        exportgraphics(f, fullfile(outdir, [name '.pdf']), 'ContentType', 'vector');
        exportgraphics(f, fullfile(outdir, [name '.png']), 'Resolution', 200);
        saved = true;
    catch
        try
            saveas(f, fullfile(outdir, [name '.png']));
            saved = true;
        catch
        end
    end
    if saved
        fprintf('  [saved figure] %s.png\n', fullfile(outdir, name));
    else
        fprintf('  [figure NOT saved (read-only fs) - values above are still valid] %s\n', name);
    end
end
