function tl = panel_map_layout(cfg, width_cm, height_cm)
%PANEL_MAP_LAYOUT New figure with one tile per state arranged like a map.
%   Place state si with nexttile(cfg.state_tile(si)); tile cfg.legend_tile
%   is free for a legend.

figure('Units', 'centimeters', 'Position', [0 0 width_cm height_cm], 'PaperPositionMode', 'auto');
tl = tiledlayout(cfg.tile_rows, cfg.tile_cols, 'TileSpacing', 'tight', 'Padding', 'tight');
end
