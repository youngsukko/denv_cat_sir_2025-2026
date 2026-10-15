%% Catalytic model results: Fig. 2A-F and Figs. S1, S2, S4-S8, plus the
%  national proportion with >= 1 prior infection reported in the Results.
%
%  Input : data/processed (R pipeline), results_big/catalytic_*_draws.mat
%  Output: output/figures/fig2*.jpg, figS1-S8*.jpg

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
obs = load_observed_data(cfg);
d   = load_catalytic_draws(cfg, 'all');

S = cfg.S; T = cfg.T; A5 = cfg.A5;
years = cfg.years;
ages  = cfg.age_lower;
n_draw = size(d.kappa_draws, 1);

col_foi    = [0.839 0.153 0.157];
col_report = [0.580 0.404 0.741];
col_model  = [0.15 0.15 0.15];
col_data   = [0.85 0.10 0.10];
lighten    = @(c, frac) c + (1 - c) * frac;
cmap_inf   = [150 150 150; 255 221 102; 255 153 51; 255 102 51; 204 0 0] / 255;  % 0..4 prior infections

% FOI causing primary infection = 4 x per-serotype FOI
foi     = 4 * d.lambda_draws;       % [n x S x T]
foi_pre = 4 * d.lambda_pre_draws;   % [n x S]

% Observed incidence per 100 by state and year, and totals over 2000-2025
inc_obs  = obs.cases_st ./ obs.pop_st * 100;   % [T x S]
inc_max  = max(inc_obs(:)) * 1.05;
cases_age_obs = zeros(S, A5);
for s = 1:S
    cases_age_obs(s, :) = sum(obs.cases((s-1)*T + (1:T), :), 1);
end
pop_total = sum(obs.pop_st, 1)';                % [S x 1] person-years, 2000-2025

%% Fig. 2A: FOI over time in three example states (P00 = pre-2000 constant FOI)
sel   = cellfun(@(c) find(strcmp(cfg.state_codes, c)), {'BA', 'SP', 'SC'});
cols  = {'r', 'b', 'k'};
lsty  = {'-', '-', '--'};
xlabs = ['P00', arrayfun(@(y) sprintf('%02d', mod(y, 100)), years, 'UniformOutput', false)];

figure('Units', 'centimeters', 'Position', [0 0 7 6], 'PaperPositionMode', 'auto');
hold on
h = gobjects(3, 1);
for k = 1:3
    si = sel(k);
    med = [median(foi_pre(:, si)), median(squeeze(foi(:, si, :)), 1)];
    h(k) = plot(1:T+1, med, lsty{k}, 'Color', cols{k}, 'LineWidth', 1.5);
end
xline(1.5, 'k:', 'LineWidth', 0.8)
xticks(1:2:T+1); xticklabels(xlabs(1:2:end)); xlim([0.5, T+1.5]); ylim([0 inf])
ylabel('Force of infection (FOI)')
legend(h, cfg.state_names(sel), 'Location', 'northeast', 'FontSize', 9); legend('boxoff')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
save_figure(cfg, 'fig2A_foi_selected_states')

%% Fig. 2B / 2C: age-specific FOI and reporting multipliers (median, 95% CrI)
multipliers = {d.kappa_draws, d.omega_draws};
mult_cols   = {col_foi, col_report};
mult_ylabel = {'FOI multiplier by age', 'Report rate multiplier'};
mult_names  = {'fig2B_kappa_age', 'fig2C_omega_age'};
for m = 1:2
    figure('Units', 'centimeters', 'Position', [0 0 7 6], 'PaperPositionMode', 'auto');
    hold on
    for a = 1:A5
        v = multipliers{m}(:, a);
        plot([a a], quantile(v, [0.025 0.975]), 'Color', mult_cols{m}, 'LineWidth', 1.5);
        plot(a, median(v), 'ko', 'LineWidth', 1.5, 'MarkerFaceColor', mult_cols{m}, 'MarkerSize', 3);
    end
    h95  = plot(nan, nan, '-', 'Color', mult_cols{m}, 'LineWidth', 1.5);
    hmed = plot(nan, nan, 'ko', 'MarkerFaceColor', mult_cols{m}, 'MarkerSize', 3);
    legend([hmed h95], {'Median', '95% CrI'}, 'Location', 'south'); legend('boxoff')
    xlim([0 18]); if m == 2; ylim([0 1]); end
    xlabel('Age group'); ylabel(mult_ylabel{m})
    set(gca, 'XTick', 1:A5, 'XTickLabel', cfg.age_labels, 'FontName', 'SansSerif', 'FontSize', 11)
    save_figure(cfg, mult_names{m})
end

%% Fig. 2D: observed vs model-predicted average annual incidence per 100,000, 2000-2025
inc_obs_state = sum(cases_age_obs, 2) ./ pop_total * 1e5;           % [S x 1]
inc_model     = sum(d.cum_cases_age, 3) ./ pop_total' * 1e5;        % [n x S]
o  = inc_obs_state(cfg.reg_order);
mm = mean(inc_model(:, cfg.reg_order), 1);
lo = quantile(inc_model(:, cfg.reg_order), 0.025, 1);
hi = quantile(inc_model(:, cfg.reg_order), 0.975, 1);

figure('Units', 'centimeters', 'Position', [0 0 17 5], 'PaperPositionMode', 'auto');
ymax = max([o(:); hi(:)]) * 1.15;
draw_region_background(cfg, ymax, 0.15)
bar(1:S, o, 0.6, 'FaceColor', col_data, 'EdgeColor', 'none', 'FaceAlpha', 0.75)
for k = 1:S
    plot([k k], [lo(k) hi(k)], 'k-', 'LineWidth', 1.2)
    plot(k, mm(k), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 0.5)
end
h1 = fill(nan, nan, col_data, 'FaceAlpha', 0.75, 'EdgeColor', 'none');
h2 = plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5, 'LineWidth', 0.5);
h3 = plot(nan, nan, 'k-', 'LineWidth', 1.2);
legend([h1 h2 h3], {'Observed (data)', 'Model mean', 'Model 95% CrI'}, 'Location', 'north', 'FontSize', 9)
legend('boxoff')
xlim([0.25 S+0.75]); ylim([0 ymax])
xticks(1:S); xticklabels(cfg.state_codes(cfg.reg_order))
ylabel({'Annual incidence', ['(2000' char(8211) '2025)']})
set(gca, 'XTickLabelRotation', 90, 'FontName', 'SansSerif', 'FontSize', 11); box on
save_figure(cfg, 'fig2D_model_calibration')

%% Fig. 2E / 2F: immune profiles at the start of 2000 and of 2026
profiles = {0, T};
profile_names = {'fig2E_immune_profile_2000', 'fig2F_immune_profile_2026'};
for f = 1:2
    panel_map_layout(cfg, 18, 22);
    for si = 1:S
        nexttile(cfg.state_tile(si))
        plot_immune_panel(seropos_draws(d, si, profiles{f}), ages, cmap_inf, cfg.state_codes{si})
    end
    nexttile(cfg.legend_tile); axis off; hold on
    bb = bar(nan(1, A5), nan(5, A5), 0.95, 'stacked', 'EdgeColor', 'none');
    for k = 1:5; bb(k).FaceColor = cmap_inf(k, :); end
    bb(1).BaseLine.Visible = 'off';
    lg = legend({'0', '1', '2', '3', '4'}, 'Location', 'northwest', 'FontSize', 11, 'NumColumns', 3);
    legend('boxoff'); title(lg, 'Prior exposure', 'FontWeight', 'normal')
    save_figure(cfg, profile_names{f})
end

%% Fig. S1: FOI by state and year (median, IQR) with reported incidence per 100
panel_map_layout(cfg, 18, 22);
for si = 1:S
    nexttile(cfg.state_tile(si))
    v = squeeze(foi(:, si, :));
    yyaxis right
    plot(years, inc_obs(:, si), 'k:', 'LineWidth', 1.2)
    ylim([0 inc_max]); ax = gca; ax.YAxis(2).Color = 'k'; ax.YAxis(1).Color = col_foi;
    yyaxis left
    fill([years, fliplr(years)], [quantile(v, 0.75, 1), fliplr(quantile(v, 0.25, 1))], ...
         col_foi, 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    plot(years, median(v, 1), 'Color', col_foi, 'LineWidth', 1.2)
    ylim([0 0.3]); xlim([years(1)-0.5 years(end)+0.5])
    xticks([2000 2010 2020]); xticklabels({'00', '10', '20'})
    title(cfg.state_names{si}, 'FontSize', 8, 'FontWeight', 'normal')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 8)
end
nexttile(cfg.legend_tile); axis off; hold on
h2 = fill(nan, nan, col_foi, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h3 = plot(nan, nan, '-', 'Color', col_foi, 'LineWidth', 1.5);
h4 = plot(nan, nan, 'k:', 'LineWidth', 1.5);
legend([h3 h2 h4], {'Median FOI', 'IQR (25-75%)', 'Incidence (data, per 100)'}, ...
       'Location', 'north', 'FontName', 'SansSerif', 'FontSize', 11); legend('boxoff')
save_figure(cfg, 'figS1_foi_panel_map')

%% Fig. S2: pre-2000 FOI by state (box chart) and region map
group_s = repmat(1:S, n_draw, 1);
group_r = repmat(cfg.region_of(cfg.reg_order)', n_draw, 1);
vals    = foi_pre(:, cfg.reg_order);

figure('Units', 'centimeters', 'Position', [0 0 9 12], 'PaperPositionMode', 'auto');
hb = boxchart(group_s(:), vals(:), 'GroupByColor', group_r(:), ...
              'MarkerStyle', 'none', 'BoxWidth', 3, 'BoxFaceAlpha', 1);
for r = 1:5; set(hb(r), 'BoxFaceColor', cfg.region_colors(r, :)); end
xticks(linspace(0.55, S+0.5, S)); xticklabels(cfg.state_codes(cfg.reg_order)); xlim([0 S+1])
ylabel('Pre-2000 FOI'); ylim([0 0.05])
view([90 90])
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
save_figure(cfg, 'figS2_pre2000_foi')

shp = load_brazil_shapes(cfg);
label_offset = containers.Map({'DF', 'PE', 'RJ', 'GO'}, {[0 0.8], [0 -0.8], [0 -0.8], [-1 0]});
figure('Units', 'centimeters', 'Position', [0 0 10 12], 'PaperPositionMode', 'auto');
hold on
for i = 1:numel(shp)
    for k = 1:numel(shp(i).parts)
        xy = shp(i).parts{k};
        patch(xy(:,1), xy(:,2), cfg.region_colors(cfg.region_of(shp(i).si), :), ...
              'EdgeColor', 'w', 'LineWidth', 0.8, 'FaceAlpha', 1)
    end
end
for i = 1:numel(shp)
    off = [0 0];
    if isKey(label_offset, shp(i).code); off = label_offset(shp(i).code); end
    text(shp(i).cx + off(1), shp(i).cy + off(2), shp(i).code, ...
         'HorizontalAlignment', 'center', 'FontName', 'SansSerif')
end
axis equal off
set(gca, 'FontName', 'SansSerif', 'FontSize', 14)
save_figure(cfg, 'figS2_region_map')

%% Fig. S4: infection-to-apparent illness ratio by infection order
figure('Units', 'centimeters', 'Position', [0 0 7 8], 'PaperPositionMode', 'auto');
hold on
for k = 1:4
    boxchart(k * ones(n_draw, 1), d.delta_draws(:, k), 'BoxFaceColor', col_report, 'MarkerStyle', 'none');
end
xlim([0.5 4.5]); xticks(1:4); xticklabels({'Primary', 'Secondary', 'Tertiary', 'Quaternary'})
ylabel('Infection to apparent illness ratio')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
save_figure(cfg, 'figS4_delta')

%% Fig. S5: apparent illness-to-case report ratio by state and year
panel_map_layout(cfg, 18, 22);
for si = 1:S
    nexttile(cfg.state_tile(si))
    v = squeeze(d.rho_draws(:, si, :));
    yyaxis right
    plot(years, inc_obs(:, si), 'k:', 'LineWidth', 1.2)
    ylim([0 inc_max]); ax = gca; ax.YAxis(2).Color = 'k'; ax.YAxis(1).Color = col_report;
    yyaxis left
    fill([years, fliplr(years)], [quantile(v, 0.975, 1), fliplr(quantile(v, 0.025, 1))], ...
         lighten(col_report, 0.65), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    fill([years, fliplr(years)], [quantile(v, 0.75, 1), fliplr(quantile(v, 0.25, 1))], ...
         col_report, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(years, median(v, 1), 'Color', col_report, 'LineWidth', 1.2)
    ylim([0 1]); xlim([years(1)-0.5 years(end)+0.5])
    xticks([2000 2010 2020]); xticklabels({'00', '10', '20'})
    title(cfg.state_names{si}, 'FontSize', 8, 'FontWeight', 'normal')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 8)
end
nexttile(cfg.legend_tile); axis off; hold on
h1 = fill(nan, nan, lighten(col_report, 0.65), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h2 = fill(nan, nan, col_report, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h3 = plot(nan, nan, '-', 'Color', col_report, 'LineWidth', 1.5);
h4 = plot(nan, nan, 'k:', 'LineWidth', 1.5);
legend([h3 h2 h1 h4], {'Median report rate', 'IQR (25-75%)', '95% CrI', 'Incidence (data, per 100)'}, ...
       'Location', 'north', 'FontName', 'SansSerif', 'FontSize', 11); legend('boxoff')
save_figure(cfg, 'figS5_report_ratio_panel_map')

%% Fig. S6: cumulative cases by age per 100,000 (state person-years), model vs data
cum_inc_model = d.cum_cases_age ./ reshape(pop_total, 1, S) * 1e5;   % [n x S x A5]
cum_inc_obs   = cases_age_obs ./ pop_total * 1e5;                    % [S x A5]

panel_map_layout(cfg, 18, 22);
for si = 1:S
    nexttile(cfg.state_tile(si))
    v = squeeze(cum_inc_model(:, si, :));
    fill([ages, fliplr(ages)], [quantile(v, 0.975, 1), fliplr(quantile(v, 0.025, 1))], ...
         lighten(col_model, 0.65), 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    fill([ages, fliplr(ages)], [quantile(v, 0.75, 1), fliplr(quantile(v, 0.25, 1))], ...
         col_model, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    plot(ages, median(v, 1), 'o-', 'Color', col_model, 'MarkerFaceColor', col_model, ...
         'MarkerSize', 2.5, 'LineWidth', 1)
    bar(ages, cum_inc_obs(si, :), 0.6, 'FaceColor', col_data, 'EdgeColor', 'none', 'FaceAlpha', 0.55)
    xlim([-2.5 82.5]); xticks([0 20 40 60 80])
    title(cfg.state_names{si}, 'FontSize', 8, 'FontWeight', 'normal')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 8)
end
nexttile(cfg.legend_tile); axis off; hold on
h1 = plot(nan, nan, 'o-', 'Color', col_model, 'MarkerFaceColor', col_model, 'LineWidth', 1);
h2 = fill(nan, nan, col_model, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h3 = fill(nan, nan, lighten(col_model, 0.65), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
h4 = bar(nan, nan, 'FaceColor', col_data, 'EdgeColor', 'none', 'FaceAlpha', 0.55);
legend([h1 h2 h3 h4], {'Model median', 'IQR (25-75%)', '95% CrI', 'Observed (data)'}, ...
       'Location', 'northwest', 'FontSize', 11); legend('boxoff')
save_figure(cfg, 'figS6_case_age_panel_map')

%% Fig. S7: reported cases per infection by state, 2000-2025
ratio = case_infection_ratio(d);                                     % [n x S]
r_med = median(ratio(:, cfg.reg_order), 1);
r_lo  = quantile(ratio(:, cfg.reg_order), 0.025, 1);
r_hi  = quantile(ratio(:, cfg.reg_order), 0.975, 1);
reg_c = cfg.region_of(cfg.reg_order);

figure('Units', 'centimeters', 'Position', [0 0 20 7], 'PaperPositionMode', 'auto');
draw_region_background(cfg, max(r_hi) * 1.3, 0.15)
for k = 1:S
    c = cfg.region_colors(reg_c(k), :);
    plot([k k], [r_lo(k) r_hi(k)], '-', 'Color', c, 'LineWidth', 1.5)
    plot(k, r_med(k), 's', 'MarkerFaceColor', c, 'MarkerEdgeColor', c * 0.7, 'MarkerSize', 7, 'LineWidth', 0.5)
end
xlim([0.25 S+0.75]); ylim([0 0.15])
xticks(1:S); xticklabels(cfg.state_codes(cfg.reg_order))
ylabel('Case / Infection ratio')
set(gca, 'XTickLabelRotation', 90, 'FontName', 'SansSerif', 'FontSize', 12); box on
save_figure(cfg, 'figS7_case_infection_ratio')

%% Fig. S8: proportion with <= 1 prior infection by age, 2000 (cyan) to 2025 (magenta)
figure('Units', 'centimeters', 'Position', [0 0 20 23], 'PaperPositionMode', 'auto', 'Color', 'white');
tiledlayout(cfg.tile_rows, cfg.tile_cols, 'TileSpacing', 'tight', 'Padding', 'compact')
for si = 1:S
    nexttile(cfg.state_tile(si)); hold on
    for ti = 1:T
        Xk = seropos_draws(d, si, ti);
        denom = Xk{1} + Xk{2} + Xk{3} + Xk{4} + Xk{5};
        c = (ti-1) / (T-1);
        plot(ages, mean(Xk{1} ./ denom, 1) + mean(Xk{2} ./ denom, 1), 'Color', [c, 1-c, 1], 'LineWidth', 1.5)
    end
    xlim([0 80]); ylim([0 1]); xticks([0 20 40 60 80])
    title(cfg.state_names{si}, 'FontSize', 8, 'FontWeight', 'normal')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 8, 'XTickLabelRotation', 0)
end
save_figure(cfg, 'figS8_immune_trend')

%% Results text: national population-weighted proportion with >= 1 prior infection
pop_2000 = zeros(S, A5);
for s = 1:S; pop_2000(s, :) = obs.pop((s-1)*T + 1, :); end
pop_2025 = obs.pop_2025;
p2000 = zeros(n_draw, 1);
p2025 = zeros(n_draw, 1);
for k = 1:n_draw
    p2000(k) = 1 - sum(squeeze(d.op_seroposall_pre_0(k, :, :)) .* pop_2000, 'all') / sum(pop_2000, 'all');
    p2025(k) = 1 - sum(squeeze(d.op_seroposall_0(k, :, T, :)) .* pop_2025, 'all') / sum(pop_2025, 'all');
end
fprintf('\nNational proportion with >= 1 prior infection (median [95%% CrI])\n');
fprintf('  start of 2000: %.1f%% [%.1f, %.1f]\n', 100 * [median(p2000), quantile(p2000, [0.025 0.975])]);
fprintf('  end of 2025:   %.1f%% [%.1f, %.1f]\n', 100 * [median(p2025), quantile(p2025, [0.025 0.975])]);
q2000 = quantile(p2000, [0.025 0.975]);
q2025 = quantile(p2025, [0.025 0.975]);
write_numbers(cfg, 'M02_catalytic', {
    'National % with >= 1 prior infection, start of 2000', 100*median(p2000), 100*q2000(1), 100*q2000(2), 'median [95% CrI], population-weighted'
    'National % with >= 1 prior infection, end of 2025',   100*median(p2025), 100*q2025(1), 100*q2025(2), 'median [95% CrI], population-weighted'});
