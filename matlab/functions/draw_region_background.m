function draw_region_background(cfg, ymax, alpha)
%DRAW_REGION_BACKGROUND Shade the x-range of each region on a state-level axis
%   whose x positions follow cfg.reg_order.

if nargin < 3; alpha = 0.15; end
b = cfg.reg_bound;
for r = 1:5
    fill([b(r) b(r+1) b(r+1) b(r)] + 0.5, [0 0 ymax ymax], cfg.region_colors(r, :), ...
         'EdgeColor', 'none', 'FaceAlpha', alpha)
    hold on
end
end
