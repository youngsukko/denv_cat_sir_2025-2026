%% ═══════════════════════════════════════════════════════════════
%  Yearly proportion of 0/1/2/3/4 prior infections, with 95% CrI
%  Population-weighted across age bands, for one state
%  Input:
%    results/catalytic_ppc_draws.mat   (op_seroposall_0 ... op_seroposall_4)
%    results/state_order.csv           (idx, uf)
%    results/pop_age_state_year.csv    (rows = (s-1)*T+t, cols = age bands)
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

target_state = 'MS';   % change this to run for a different state

S = 27; T = 26; A5 = 17;
years       = 2000:2025;
agerangeidx = 0:5:80;

agelabels = arrayfun(@(i) sprintf('%d-',agerangeidx(i)), 1:16, 'UniformOutput',false);
agelabels{17} = '80+';

cmap = [150 150 150; 255 221 102; 255 153 51; 255 102 51; 204 0 0] / 255;

pct025 = @(x) quantile(x, 0.025);
pct975 = @(x) quantile(x, 0.975);

%% ── 1. Load draws and locate the state ───────────────────────────
load('results/catalytic_ppc_draws.mat', ...
     'op_seroposall_0','op_seroposall_1','op_seroposall_2', ...
     'op_seroposall_3','op_seroposall_4');
% each op_seroposall_k: [n_draw x S x T x A5]

seropos = {op_seroposall_0, op_seroposall_1, op_seroposall_2, ...
           op_seroposall_3, op_seroposall_4};
n_draw  = size(op_seroposall_0, 1);
n_order = 5;

state_order = readtable('results/state_order.csv');   % columns: idx, uf
s_idx = state_order.idx(strcmp(state_order.uf, target_state));
if isempty(s_idx)
    error('State "%s" not found in state_order.csv', target_state);
end

%% ── 2. Population, with the age-column order corrected ───────────
% The CSV columns come out of R's sort(unique(age_group)), which is
% alphabetical, so "5-9" lands between "45-49" and "50-54". Reorder to
% true age order before using the population as age weights.
pop_mat = table2array(readtable('results/pop_age_state_year.csv'));   % [S*T x A5]

wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
wrong_order_short = cellfun(@(s) strsplit(s,'-'), wrong_order, 'UniformOutput', false);
wrong_order_short = cellfun(@(c) [c{1} '-'], wrong_order_short, 'UniformOutput', false);
wrong_order_short{strcmp(wrong_order,'80+')} = '80+';

[~, reorder_idx] = ismember(agelabels, wrong_order_short);
assert(all(reorder_idx > 0), 'CSV column mapping failed — check wrong_order');
pop_mat = pop_mat(:, reorder_idx);

pop_state = pop_mat((s_idx-1)*T + (1:T), :);   % [T x A5]

%% ── 3. Population-weighted average across age bands, per draw ─────
weighted_draws = zeros(n_draw, T, n_order);

for t = 1:T
    w = pop_state(t, :);
    w = w / sum(w);
    for k = 1:n_order
        draws_ord = squeeze(seropos{k}(:, s_idx, t, :));   % [n_draw x A5]
        weighted_draws(:, t, k) = draws_ord * w';
    end
end

% Normalize per draw so the five orders sum to 1
weighted_draws = weighted_draws ./ sum(weighted_draws, 3);

%% ── 4. Median and 95% CrI, wide table (rows = years, cols = 0-4) ──
med_pct = squeeze(median(weighted_draws, 1)) * 100;   % [T x order]
lo_pct  = squeeze(pct025(weighted_draws))     * 100;
hi_pct  = squeeze(pct975(weighted_draws))     * 100;

cell_str = strings(T, n_order);
for t = 1:T
    for k = 1:n_order
        cell_str(t, k) = sprintf('%.2f {%.2f-%.2f}', ...
            med_pct(t, k), lo_pct(t, k), hi_pct(t, k));
    end
end

result = table(years', cell_str(:,1), cell_str(:,2), cell_str(:,3), ...
               cell_str(:,4), cell_str(:,5), ...
               'VariableNames', {'Year','0','1','2','3','4'});
disp(result);

out_file = sprintf('results/seropos_by_year_%s.csv', target_state);
writetable(result, out_file);
fprintf('Saved: %s\n', out_file);

%% ── 5. Stacked bar chart with error bars ─────────────────────────
% CrI of the CUMULATIVE stack boundary, computed per draw
cum_draws = cumsum(weighted_draws, 3);                % [n_draw x T x order]
cum_med   = squeeze(median(cum_draws, 1)) * 100;      % [T x order]
cum_lo    = squeeze(pct025(cum_draws))    * 100;
cum_hi    = squeeze(pct975(cum_draws))    * 100;

% Segment heights that reproduce the median cumulative boundaries exactly
seg_height = [cum_med(:,1), diff(cum_med, 1, 2)];

figure('Units','centimeters','Position',[0 0 18 8],'PaperPositionMode','auto','Color','white');
hold on
bb = bar(years, seg_height, 0.95, 'stacked', 'EdgeColor', 'none');
for k = 1:n_order; bb(k).FaceColor = cmap(k,:); end

% Error bars sit ON the stack boundaries (last boundary is always 100%)
for k = 1:(n_order-1)
    errorbar(years, cum_med(:,k), ...
             cum_med(:,k) - cum_lo(:,k), cum_hi(:,k) - cum_med(:,k), ...
             'k', 'LineStyle','none', 'CapSize', 2, 'LineWidth', 0.6);
end

xlim([years(1)-0.5 years(end)+0.5]); ylim([0 100])
xticks(2000:4:2025)
xlabel('Year'); ylabel('Proportion (%)')
% title(sprintf('Infection count distribution — %s', target_state), ...
%       'FontSize', 8, 'FontWeight', 'normal')
lg = legend(bb, {'0','1','2','3','4'}, 'Location','northoutside', 'NumColumns', 5);
legend('boxoff'); title(lg, 'Prior exposure', 'FontWeight', 'normal')
set(gca, 'FontName','SansSerif', 'FontSize', 11)

mkdir('figures')
exportgraphics(gcf, sprintf('figures/fig_seropos_by_year_%s.jpg', target_state), ...
               'Resolution', 600);
fprintf('Saved: figures/fig_seropos_by_year_%s.jpg\n', target_state);

%% ═══════════════════════════════════════════════════════════════
%  Proportion of 0/1/2/3/4 prior infections by age group, most recent year
%  Stacked bars with 95% CrI error bars on the stack boundaries
%  Input:
%    results/catalytic_ppc_draws.mat   (op_seroposall_0 ... op_seroposall_4)
%    results/state_order.csv           (idx, uf)
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

target_state = 'MS';   % change this to run for a different state

S = 27; T = 26; A5 = 17;
years       = 2000:2025;
agerangeidx = 0:5:80;

agelabels = arrayfun(@(i) sprintf('%d-',agerangeidx(i)), 1:16, 'UniformOutput',false);
agelabels{17} = '80+';

cmap = [150 150 150; 255 221 102; 255 153 51; 255 102 51; 204 0 0] / 255;

pct025 = @(x) quantile(x, 0.025);
pct975 = @(x) quantile(x, 0.975);

ti = T;   % most recent year

%% ── 1. Load draws and locate the state ───────────────────────────
load('results/catalytic_ppc_draws.mat', ...
     'op_seroposall_0','op_seroposall_1','op_seroposall_2', ...
     'op_seroposall_3','op_seroposall_4');
% each op_seroposall_k: [n_draw x S x T x A5]

seropos = {op_seroposall_0, op_seroposall_1, op_seroposall_2, ...
           op_seroposall_3, op_seroposall_4};
n_draw  = size(op_seroposall_0, 1);
n_order = 5;

state_order = readtable('results/state_order.csv');   % columns: idx, uf
s_idx = state_order.idx(strcmp(state_order.uf, target_state));
if isempty(s_idx)
    error('State "%s" not found in state_order.csv', target_state);
end

%% ── 2. Pull the age profile for the most recent year, per draw ────
prof_draws = zeros(n_draw, A5, n_order);
for k = 1:n_order
    prof_draws(:, :, k) = squeeze(seropos{k}(:, s_idx, ti, :));   % [n_draw x A5]
end

% Normalize per draw so the five orders sum to 1 at each age band
prof_draws = prof_draws ./ sum(prof_draws, 3);

%% ── 3. Median and 95% CrI across draws ───────────────────────────
med_pct = squeeze(median(prof_draws, 1)) * 100;   % [A5 x order]
lo_pct  = squeeze(pct025(prof_draws))    * 100;
hi_pct  = squeeze(pct975(prof_draws))    * 100;

%% ── 4. CSV (rows = age group, cols = order 0-4) ──────────────────
cell_str = strings(A5, n_order);
for b = 1:A5
    for k = 1:n_order
        cell_str(b, k) = sprintf('%.2f {%.2f-%.2f}', ...
            med_pct(b, k), lo_pct(b, k), hi_pct(b, k));
    end
end

result = table(agelabels', cell_str(:,1), cell_str(:,2), cell_str(:,3), ...
               cell_str(:,4), cell_str(:,5), ...
               'VariableNames', {'AgeGroup','0','1','2','3','4'});
disp(result);

out_file = sprintf('results/seropos_by_age_%s_%d.csv', target_state, years(ti));
writetable(result, out_file);
fprintf('Saved: %s\n', out_file);

%% ── 5. Stacked bar chart with error bars ─────────────────────────
% Stack order is reversed: 4 infections at the bottom, 0 at the top.
% Colours stay tied to their infection count, so cmap is flipped to match.
prof_rev = flip(prof_draws, 3);       % order dim now runs 4,3,2,1,0
cmap_rev = flipud(cmap);

% CrI of the CUMULATIVE stack boundary, computed per draw
cum_draws = cumsum(prof_rev, 3);                      % [n_draw x A5 x order]
cum_med   = squeeze(median(cum_draws, 1)) * 100;      % [A5 x order]
cum_lo    = squeeze(pct025(cum_draws))    * 100;
cum_hi    = squeeze(pct975(cum_draws))    * 100;

% Segment heights that reproduce the median cumulative boundaries exactly
seg_height = [cum_med(:,1), diff(cum_med, 1, 2)];

figure('Units','centimeters','Position',[0 0 9 7],'PaperPositionMode','auto','Color','white');
hold on
bb = bar(agerangeidx, seg_height, 0.95, 'stacked', 'EdgeColor', 'none');
for k = 1:n_order; bb(k).FaceColor = cmap_rev(k,:); end

% Error bars sit ON the stack boundaries (last boundary is always 100%),
% dodged slightly along x so the four boundaries don't overlap
dodge = 0.8;   % x offset step, in age-axis units (bands are 5 apart)
for k = 1:(n_order-1)
    xk = agerangeidx + (k - (n_order-1)/2 - 0.5) * dodge;
    errorbar(xk, cum_med(:,k), ...
             cum_med(:,k) - cum_lo(:,k), cum_hi(:,k) - cum_med(:,k), ...
             'k', 'LineStyle','none', 'CapSize', 2, 'LineWidth', 0.6);
end

xlim([-2.5 82.5]); ylim([0 100])
xticks(agerangeidx); xticklabels(agelabels)
xlabel('Age group'); ylabel('Proportion (%)')
% title(sprintf('%s — %d', target_state, years(ti)), ...
%       'FontSize', 8, 'FontWeight', 'normal')
lg = legend(flip(bb), {'0','1','2','3','4'}, 'Location','northoutside', 'NumColumns', 5);
legend('boxoff'); title(lg, 'Prior exposure', 'FontWeight', 'normal')
set(gca, 'FontName','SansSerif', 'FontSize', 11, 'XTickLabelRotation', 90)

mkdir('figures')
exportgraphics(gcf, sprintf('figures/fig_seropos_by_age_%s_%d.jpg', target_state, years(ti)), ...
               'Resolution', 600);
fprintf('Saved: figures/fig_seropos_by_age_%s_%d.jpg\n', target_state, years(ti));