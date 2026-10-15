function plot_region_map_inset(cfg)
%PLOT_REGION_MAP_INSET Small map of Brazil coloured by region (Fig. 3F inset).
%   Drawn at its final printed size so that font sizes are preserved when
%   placed at 100%. Small north-eastern coastal states are labelled with a
%   fan of leader lines; DF, ES, RJ and SC have individual leader lines.

map_w_cm = 3.3;        % printed map width (cm)
lab_fs   = 6;          % label font size (pt)
xl1      = -74.3;      % west boundary (longitude)
yl       = [-34 5.8];  % south/north boundaries (latitude)

fan_codes = {'RN', 'PB', 'PE', 'AL', 'SE'};
free_lab  = {'DF', -2.6, -2.6, 'center';     % {code, dx, dy, alignment}
             'ES',  4.5, -0.5, 'left';
             'RJ',  0.0, -3.2, 'center';
             'SC',  5.4, -0.9, 'left'};
in_shift  = {'GO', -1.2, 1.8};                % label shift inside the state

raw  = shaperead(cfg.file_shape);
n    = numel(raw);
code = {raw.SIGLA_UF};

% region colour, label anchor (centroid of the largest ring) and area for drawing order
reg = zeros(n, 1); cx = zeros(n, 1); cy = zeros(n, 1); xmax = zeros(n, 1); area = zeros(n, 1);
for i = 1:n
    reg(i) = cfg.region_of(strcmp(cfg.state_codes, code{i}));
    X = raw(i).X(:); Y = raw(i).Y(:);
    area(i) = (max(X) - min(X)) * (max(Y) - min(Y));
    brk = [0; find(isnan(X)); numel(X) + 1];
    best = 0;
    for k = 1:numel(brk) - 1
        xs = X(brk(k)+1 : brk(k+1)-1); ys = Y(brk(k)+1 : brk(k+1)-1);
        if numel(xs) < 3, continue, end
        xa = [xs; xs(1)]; ya = [ys; ys(1)];
        cr = xa(1:end-1) .* ya(2:end) - xa(2:end) .* ya(1:end-1);
        a  = sum(cr) / 2;
        if abs(a) > best
            best    = abs(a);
            cx(i)   = sum((xa(1:end-1) + xa(2:end)) .* cr) / (6*a);
            cy(i)   = sum((ya(1:end-1) + ya(2:end)) .* cr) / (6*a);
            xmax(i) = max(xs);
        end
    end
end

% x position of the coastal label column and matching east boundary
[~, fan_i] = ismember(fan_codes, code);
[~, ord]   = sort(cy(fan_i), 'descend');
fan_i      = fan_i(ord);
x_lab      = max(xmax(fan_i)) + 1.3;
Lfrac      = (lab_fs / 72 * 2.54) / map_w_cm;      % label height / map width
xr         = (x_lab + 0.3 - 1.8 * Lfrac * xl1) / (1 - 1.8 * Lfrac);
lab_h      = Lfrac * (xr - xl1);                    % label height (degrees)
map_h_cm   = map_w_cm * diff(yl) / (xr - xl1);

figure('Units', 'centimeters', 'Position', [0 0 map_w_cm map_h_cm], 'PaperPositionMode', 'auto');
axes('Position', [0 0 1 1]); hold on

% large states first so that small ones (e.g. DF) stay on top
[~, draw_ord] = sort(area, 'descend');
for ii = draw_ord(:)'
    X = raw(ii).X(:); Y = raw(ii).Y(:);
    brk = [0; find(isnan(X)); numel(X) + 1];
    for k = 1:numel(brk) - 1
        xs = X(brk(k)+1 : brk(k+1)-1); ys = Y(brk(k)+1 : brk(k+1)-1);
        if numel(xs) < 3, continue, end
        patch(xs, ys, cfg.region_colors(reg(ii), :), 'EdgeColor', 'w', 'LineWidth', 0.4);
    end
end
axis equal off
xlim([xl1 xr]); ylim(yl)

txt = @(x, y, c, ha) text(x, y, c, 'FontName', 'SansSerif', 'FontSize', lab_fs, ...
                          'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
dot = @(x, y) plot(x, y, 'ko', 'MarkerSize', 1.8, 'MarkerFaceColor', 'k', 'LineWidth', 0.1);

% (1) labels inside the states
for i = 1:n
    c = code{i};
    if any(strcmp(c, fan_codes)) || any(strcmp(c, free_lab(:, 1))), continue, end
    dxy = [0 0];
    j = find(strcmp(in_shift(:, 1), c));
    if ~isempty(j), dxy = [in_shift{j, 2}, in_shift{j, 3}]; end
    txt(cx(i) + dxy(1), cy(i) + dxy(2), c, 'center');
end

% (2) coastal north-eastern states: fan of labels around their mean latitude
nf    = numel(fan_i);
y_mid = mean(cy(fan_i));
gap   = 1.15 * lab_h;
for k = 1:nf
    i   = fan_i(k);
    y_l = y_mid + ((nf - 1)/2 - (k - 1)) * gap;
    plot([cx(i) x_lab - 0.3], [cy(i) y_l], 'k-', 'LineWidth', 0.4)
    dot(cx(i), cy(i))
    txt(x_lab, y_l, code{i}, 'left');
end

% (3) individual leader lines
for j = 1:size(free_lab, 1)
    i  = find(strcmp(code, free_lab{j, 1}));
    lx = cx(i) + free_lab{j, 2}; ly = cy(i) + free_lab{j, 3};
    if strcmp(free_lab{j, 4}, 'left')
        ex = lx - 0.3; ey = ly;
    else   % stop the line at the edge of the label box
        vx = cx(i) - lx; vy = cy(i) - ly; vn = hypot(vx, vy);
        ux = vx / vn;    uy = vy / vn;
        tb = 1.15 * min(0.75 * lab_h / max(abs(ux), 1e-6), 0.45 * lab_h / max(abs(uy), 1e-6));
        ex = lx + ux * tb; ey = ly + uy * tb;
    end
    plot([cx(i) ex], [cy(i) ey], 'k-', 'LineWidth', 0.4)
    dot(cx(i), cy(i))
    txt(lx, ly, free_lab{j, 1}, free_lab{j, 4});
end
end
