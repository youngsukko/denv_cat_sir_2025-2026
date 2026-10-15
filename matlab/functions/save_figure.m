function save_figure(cfg, name, fig)
%SAVE_FIGURE Export a figure to matlab/output/figures/<name>.png at 600 dpi.

if nargin < 3; fig = gcf; end
f = fullfile(cfg.fig_dir, [name '.png']);
exportgraphics(fig, f, 'Resolution', 600)
fprintf('Saved: %s\n', f);
end
