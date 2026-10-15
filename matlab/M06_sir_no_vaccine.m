%% Burden of disease without vaccination, 2026-2030: Fig. 3, Fig. S12,
%  Table S12, and the correlations reported in the Results.
%
%  Input : results/sir_inputs.mat (M05), data/processed (R pipeline)
%  Output: output/figures/fig3*.jpg, figS12_bod_age_by_region.jpg
%          output/tables/tableS12_bod_per100k.csv
%          output/sim/sim_no_vaccine.mat

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
obs    = load_observed_data(cfg);
sir_in = load(cfg.file_sir_inputs);
p      = sir_base_parameters(cfg, sir_in, 'shen');   % VE is irrelevant without vaccination

S = cfg.S; T = cfg.T; A5 = cfg.A5;
ny = cfg.n_sim_years;
[pop_state, prop_le1] = sir_state_population(sir_in);

%% Simulation for the three R0 scenarios
bod_names = {'case', 'hosp', 'dead', 'daly'};
bod = struct();
for m = 1:4; bod.(bod_names{m}) = zeros(S, A5, 3); end
for rc = 1:3
    fprintf('No vaccination, R0 %s\n', cfg.R0_names{rc});
    for si = 1:S
        out = sir_simulate(p, sir_in, si, rc);
        for m = 1:4; bod.(bod_names{m})(si, :, rc) = out.(bod_names{m}); end
    end
end
save(fullfile(cfg.sim_dir, 'sim_no_vaccine.mat'), 'bod', 'pop_state', 'prop_le1')

% annual burden per 100,000 by state: [S x 3]
annual_rate = @(M) squeeze(sum(M, 2)) ./ pop_state * 1e5 / ny;

%% Fig. 3A-D: annual burden per 100,000 by state (Q50, error bars Q25-Q75)
ylims   = [8000 250 15 400];
ylabels = {'Cases per 100K', 'Hospitalization per 100K', 'Deaths per 100K', 'DALYs per 100K'};
figure('Units', 'centimeters', 'Position', [0 0 18 10], 'PaperPositionMode', 'auto');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
for m = 1:4
    nexttile(m)
    v = annual_rate(bod.(bod_names{m}));
    v = v(cfg.reg_order, :);
    ymax = ylims(m);
    draw_region_background(cfg, ymax, 0.5)
    bar(v(:, 2), 0.7, 'FaceColor', 'k', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5)
    errorbar(1:S, v(:, 2), v(:, 2) - v(:, 1), v(:, 3) - v(:, 2), ...
             'LineStyle', 'none', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'CapSize', 3)
    xlim([0.25 S+0.75]); ylim([0 ymax])
    xticks(1:S); xticklabels(cfg.state_codes(cfg.reg_order))
    ylabel(ylabels{m})
    set(gca, 'FontName', 'SansSerif', 'FontSize', 10, 'XTickLabelRotation', 90)
end
save_figure(cfg, 'fig3A-D_bod_by_state')

%% Fig. 3E-F: annual DALYs vs past incidence and vs proportion with <= 1 prior infection
daly = annual_rate(bod.daly);                                         % [S x 3]
past_inc = sum(obs.cases_st(T-9:T, :), 1)' ./ pop_state * 1e5;        % 2016-2025
xvals  = {past_inc, prop_le1};
xlabs  = {'Past 10-year cumulative incidences (2016-2025)', ['Proportion ' char(8804) '1 infection as of 2026']};

figure('Units', 'centimeters', 'Position', [0 0 18 7], 'PaperPositionMode', 'auto');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
fprintf('\nCorrelation with annual DALYs per 100K (Q50)\n');
corr_rp = zeros(2, 2);   % [r p] for each indicator
for k = 1:2
    nexttile(k)
    x  = xvals{k};
    pf = polyfit(x, daly(:, 2), 1);
    xf = linspace(min(x) * 0.95, max(x) * 1.05, 100);
    plot(xf, polyval(pf, xf), 'k:', 'LineWidth', 1.5); hold on
    h = gobjects(5, 1);
    for r = 1:5
        idx = cfg.region_of == r;
        h(r) = plot(x(idx), daly(idx, 2), 'o', 'Color', cfg.region_colors(r, :), ...
                    'MarkerFaceColor', cfg.region_colors(r, :), 'MarkerSize', 4, 'LineWidth', 1);
        errorbar(x(idx), daly(idx, 2), daly(idx, 2) - daly(idx, 1), daly(idx, 3) - daly(idx, 2), ...
                 'LineStyle', 'none', 'Color', cfg.region_colors(r, :) * 0.7, 'LineWidth', 0.8, 'CapSize', 3);
    end
    [rv, pv] = corr(x, daly(:, 2));
    if pv < 0.001; ptxt = 'p < 0.001'; elseif pv < 0.01; ptxt = sprintf('p = %.3f', pv); else; ptxt = sprintf('p = %.2f', pv); end
    text(0.05, 0.75, {sprintf('r = %.2f', rv); ptxt}, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'FontName', 'SansSerif', 'FontSize', 9)
    fprintf('  %s: r = %.3f, p = %.4f\n', xlabs{k}, rv, pv);
    corr_rp(k, :) = [rv pv];
    xlabel(xlabs{k}, 'FontSize', 9)
    ylabel(['Annual DALY per 100K, 2026' char(8211) '2030'], 'FontSize', 9)
    set(gca, 'FontName', 'SansSerif', 'FontSize', 9); box on
    if k == 1
        legend(h, cfg.region_names, 'Location', 'northeast', 'FontSize', 8); legend('boxoff')
    end
end
save_figure(cfg, 'fig3E-F_daly_scatter')

write_numbers(cfg, 'M06_no_vaccine', {
    'Correlation r: past 10-year incidence vs annual DALYs per 100K', corr_rp(1, 1), NaN, NaN, sprintf('p = %.4g', corr_rp(1, 2))
    'Correlation r: proportion <= 1 prior infection (2026) vs annual DALYs per 100K', corr_rp(2, 1), NaN, NaN, sprintf('p = %.4g', corr_rp(2, 2))});
writetable(table(cfg.state_codes(cfg.reg_order)', 100 * prop_le1(cfg.reg_order), past_inc(cfg.reg_order), ...
                 'VariableNames', {'state', 'pct_le1_prior_infection_2026', 'past10y_incidence_per100k'}), ...
           fullfile(cfg.table_dir, 'state_indicators_2026.csv'))

plot_region_map_inset(cfg);
exportgraphics(gcf, fullfile(cfg.fig_dir, 'fig3F_inset_region_map.png'), 'Resolution', 600, ...
               'BackgroundColor', 'none');

%% Fig. S12: age distribution of the burden by region (Q50)
titles = {'Case proportion by age', 'Hospitalization proportion by age', ...
          'Death proportion by age', 'DALY proportion by age'};
figure('Units', 'centimeters', 'Position', [0 0 18 10], 'PaperPositionMode', 'auto');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
for m = 1:4
    nexttile(m)
    M = bod.(bod_names{m})(:, :, 2);
    reg_age = zeros(5, A5);
    for r = 1:5
        idx = cfg.region_of == r;
        reg_age(r, :) = sum(M(idx, :), 1) / sum(pop_state(idx)) * 1e5 / ny;
    end
    hp = plot(cfg.age_lower, (reg_age ./ sum(reg_age, 2))', 'LineWidth', 1.5);
    for r = 1:5; hp(r).Color = cfg.region_colors(r, :); end
    xticks(cfg.age_lower); xticklabels(cfg.age_labels); xlim([-0.5 80.5])
    title(titles{m}, 'FontWeight', 'normal')
    if m == 3
        lg = legend(cfg.region_names, 'Location', 'northwest'); legend('boxoff'); lg.NumColumns = 2;
    end
    set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
end
save_figure(cfg, 'figS12_bod_age_by_region')

%% Table S12: annual burden per 100,000 by state, median [Q25, Q75] over R0 scenarios
c = cell(S, 4);
for m = 1:4
    v = annual_rate(bod.(bod_names{m}));
    for k = 1:S
        si = cfg.reg_order(k);
        c{k, m} = sprintf('%.2f [%.2f, %.2f]', v(si, 2), v(si, 1), v(si, 3));
    end
end
writetable(cell2table(c, 'VariableNames', {'case_per100k', 'hosp_per100k', 'fatal_per100k', 'daly_per100k'}, ...
           'RowNames', cfg.state_codes(cfg.reg_order)), fullfile(cfg.table_dir, 'tableS12_bod_per100k.csv'), ...
           'WriteRowNames', true)
