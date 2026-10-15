%% Baseline vaccination scenario: 5 million children (5-14 y) per year,
%  distributed across states as in Brazil 2024/25. Fig. 4 (Shen VE),
%  Fig. S13 (Daniels vs Shen VE) and Table S10.
%
%  Input : results/sir_inputs.mat (M05), matlab/data/vaccine_doses_by_state_2024_2025.csv
%  Output: output/figures/fig4_bod_averted_baseline.jpg, figS13_*.jpg
%          output/tables/tableS10_vaccine_distribution.csv
%          output/sim/sim_baseline_vaccination.mat

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
obs    = load_observed_data(cfg);
sir_in = load(cfg.file_sir_inputs);

S = cfg.S; A5 = cfg.A5;
doses_annual = brazil_dose_allocation(cfg);
vactot = cfg.annual_doses * cfg.n_sim_years / 1000;   % thousands vaccinated over the horizon
bod_names = {'case', 'hosp', 'dead', 'daly'};

%% Simulations: no vaccine (once) and Brazil allocation for each VE setting
base = struct(); vacc = struct();
for m = 1:4
    base.(bod_names{m}) = zeros(S, A5, 3);
    vacc.(bod_names{m}) = zeros(S, A5, 3, 2);
end
for vi = 1:2
    p = sir_base_parameters(cfg, sir_in, cfg.ve_names{vi});
    for rc = 1:3
        fprintf('VE %s, R0 %s\n', cfg.ve_labels{vi}, cfg.R0_names{rc});
        for si = 1:S
            if vi == 1
                out = sir_simulate(p, sir_in, si, rc);
                for m = 1:4; base.(bod_names{m})(si, :, rc) = out.(bod_names{m}); end
            end
            vacc_in = struct('vacc', doses_annual(si), 'vacorder', cfg.vacorder.children);
            out = sir_simulate(p, sir_in, si, rc, vacc_in);
            for m = 1:4; vacc.(bod_names{m})(si, :, rc, vi) = out.(bod_names{m}); end
        end
    end
end
save(fullfile(cfg.sim_dir, 'sim_baseline_vaccination.mat'), 'base', 'vacc', 'doses_annual')

% averted per 1,000 vaccinated, national, by age: [3 (R0) x A5 x 2 (VE)]
averted = struct();
for m = 1:4
    averted.(bod_names{m}) = squeeze(sum(base.(bod_names{m}) - vacc.(bod_names{m}), 1)) / vactot;
    averted.(bod_names{m}) = permute(averted.(bod_names{m}), [2 1 3]);
end

%% Fig. 4: BoD averted per 1,000 vaccinated by age (Shen VE, Q50 with Q25-Q75)
dark_green  = [0.0 0.55 0.0];
light_green = [0.6 1.0 0.6];
shade_x     = [2.5 12.5 12.5 2.5];   % targeted ages 5-14
ylabels = {{'Cases averted', '(per 1000 vaccinated)'}, {'Hospitalizations averted', '(per 1000 vaccinated)'}, ...
           {'Deaths averted', '(per 1000 vaccinated)'}, {'DALYs averted', '(per 1000 vaccinated)'}};

figure('Units', 'centimeters', 'Position', [0 0 18 12], 'PaperPositionMode', 'auto');
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact')
for m = 1:4
    a = averted.(bod_names{m})(:, :, 1);   % Shen
    q25 = a(1, :); q50 = a(2, :); q75 = a(3, :);
    ymax = max(q75) * 1.3;
    nexttile(m)
    fill(shade_x, [0 0 ymax ymax], light_green, 'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on
    bar(cfg.age_lower, q50, 0.8, 'FaceColor', dark_green, 'EdgeColor', 'none')
    errorbar(cfg.age_lower, q50, q50 - q25, q75 - q50, 'LineStyle', 'none', ...
             'Color', [0.1 0.1 0.1], 'LineWidth', 1.0, 'CapSize', 3)
    xlim([-2.5 82.5]); ylim([0 ymax])
    xticks(cfg.age_lower); xticklabels(cfg.age_labels)
    set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
    ylabel(ylabels{m}, 'FontSize', 9)
    tot = sum(a, 2);
    if m == 3; fmt = 'Total: %.4f [%.4f, %.4f]'; else; fmt = 'Total: %.2f [%.2f, %.2f]'; end
    text(0.995, 0.95, sprintf(fmt, tot(2), min(tot([1 3])), max(tot([1 3]))), 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontName', 'SansSerif', 'FontSize', 9)
end
save_figure(cfg, 'fig4_bod_averted_baseline')

%% Fig. S13: Daniels vs Shen VE (one figure per outcome, by age and total)
col_ve = [0, 62, 116; 225, 36, 42] / 255;   % Daniels, Shen
ve_order = [2 1];                           % plot Daniels first
s13_ylabels = {'Cases averted', 'Hospitalizations averted', 'Fatality averted', 'DALY averted'};
fprintf('\nTotal BoD averted per 1,000 vaccinated (Q50 [Q25, Q75])\n');
for m = 1:4
    a = averted.(bod_names{m})(:, :, ve_order);   % [3 x A5 x 2]
    q50 = squeeze(a(2, :, :))'; q25 = squeeze(a(1, :, :))'; q75 = squeeze(a(3, :, :))';   % [2 x A5]
    tot = squeeze(sum(a, 2));                                                            % [3 x 2]
    ymax = max(q75(:)) * 1.15;

    figure('Units', 'centimeters', 'Position', [0 0 9 5], 'PaperPositionMode', 'auto');
    tiledlayout(1, 4, 'TileSpacing', 'compact')
    nexttile([1 3])
    fill(shade_x, [0 0 ymax ymax], light_green, 'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on
    bb = bar(cfg.age_lower, q50', 'EdgeColor', 'none');
    for k = 1:2
        bb(k).FaceColor = col_ve(k, :);
        xk = bb(k).XEndPoints;
        errorbar(xk, q50(k, :), q50(k, :) - q25(k, :), q75(k, :) - q50(k, :), '.', ...
                 'Color', 'k', 'LineWidth', 1.0, 'CapSize', 3, 'LineStyle', 'none')
    end
    yline(0, 'k:', 'LineWidth', 0.8)
    xticks(cfg.age_lower); xticklabels(cfg.age_labels); xlim([-2.5 82.5]); ylim([0 ymax])
    ylabel(s13_ylabels{m})
    set(gca, 'FontName', 'SansSerif', 'FontSize', 10)
    legend(bb, cfg.ve_labels(ve_order), 'Location', 'northeast', 'FontSize', 8); legend('boxoff')

    nexttile(4); hold on
    bt = bar(1, tot(2, :), 'EdgeColor', 'none');
    for k = 1:2
        bt(k).FaceColor = col_ve(k, :);
        errorbar(bt(k).XEndPoints(1), tot(2, k), tot(2, k) - tot(1, k), tot(3, k) - tot(2, k), '.', ...
                 'Color', 'k', 'LineWidth', 1.0, 'CapSize', 3, 'LineStyle', 'none')
    end
    xlim([0.5 1.5]); xticks(1); xticklabels({'Total'}); ylim([0 max(tot(3, :)) * 1.4])
    set(gca, 'FontName', 'SansSerif', 'FontSize', 10, 'XTickLabelRotation', 0)
    save_figure(cfg, ['figS13_bod_averted_ve_' bod_names{m}])

    for k = 1:2
        fprintf('  %-5s %-8s %.4f [%.4f, %.4f]\n', bod_names{m}, cfg.ve_labels{ve_order(k)}, ...
                tot(2, k), min(tot([1 3], k)), max(tot([1 3], k)));
    end
end

%% Table S10: vaccine distribution by state and per 1,000 children (5-14 y)
children_pop = sum(obs.pop_2025(:, cfg.target_ages.children), 2);
c = cell(S, 2);
for k = 1:S
    si = cfg.reg_order(k);
    c{k, 1} = sprintf('%.0f', doses_annual(si));
    c{k, 2} = sprintf('%.2f', doses_annual(si) / children_pop(si) * 1000);
end
writetable(cell2table(c, 'VariableNames', {'doses_total', 'doses_per1000_children'}, ...
           'RowNames', cfg.state_codes(cfg.reg_order)), ...
           fullfile(cfg.table_dir, 'tableS10_vaccine_distribution.csv'), 'WriteRowNames', true)

%% Numbers quoted in the text
rows = {};
for m = 1:4
    tot = squeeze(sum(averted.(bod_names{m}), 2));   % [3 (R0) x 2 (VE)]
    for vi = 1:2
        rows(end+1, :) = {sprintf('%s averted per 1,000 vaccinated, %s VE', bod_names{m}, cfg.ve_labels{vi}), ...
                          tot(2, vi), min(tot([1 3], vi)), max(tot([1 3], vi)), 'Q50 [range over R0 Q25-Q75]'}; %#ok<SAGROW>
    end
    rel = (tot(:, 2) ./ tot(:, 1) - 1) * 100;
    rows(end+1, :) = {sprintf('%s averted, Daniels relative to Shen (%%)', bod_names{m}), ...
                      rel(2), min(rel([1 3])), max(rel([1 3])), 'same R0 scenario, Q50 [range]'}; %#ok<SAGROW>
end
rows(end+1, :) = {'Doses over the horizon / population aged 5-14 in 2025 (%)', ...
                  100 * cfg.annual_doses * cfg.n_sim_years / sum(children_pop), NaN, NaN, ...
                  sprintf('%g doses per year x %d years', cfg.annual_doses, cfg.n_sim_years)};
write_numbers(cfg, 'M08_baseline_vaccination', rows)
