%% Fig. S9: model-based seroprevalence vs published serosurveys (refs in Appendix S3-4).
%  Model: posterior median and 95% CrI of the proportion with >= 1 prior
%  infection at the start of the survey year (pre-2000 surveys use the
%  baseline at the start of 2000). Survey values are drawn as
%    'step'     one segment per age interval
%    'multibin' min-max band and median per age interval (several sites)
%    'range'    a single min-max band over all ages
%
%  Input : results_big/catalytic_ppc_draws.mat
%  Output: output/figures/figS9_serosurvey_<region>.jpg

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
d   = load_catalytic_draws(cfg, 'ppc');

n_draw = size(d.op_seroposall_0, 1);
lo_idx = max(1, round(n_draw * 0.025));
hi_idx = min(n_draw, round(n_draw * 0.975));
col_model       = [0.15 0.15 0.15];
col_model_light = col_model + (1 - col_model) * 0.65;
col_survey      = [0 1 0];

% {state, survey year, age breaks, seroprevalence, style}
north = {
    'PA', 1997, [0 5 10 15 25 35 45 55 82], [26.7 37.5 43.3 41.3 45.6 43.3 37.4 52.5]/100, 'step'
    'AM', 2014, [0 17 48 80], [41.9 92.4 99.5]/100, 'step'};
northeast = {
    'CE', 1994, [0 80], [71 50 21 25 31 55 23 35 58]/100, 'range'
    'MA', 1996, [0 15 25 35 45 55 80], 1 - [76.4 60 59 54.4 57.3 44.8]/100, 'step'
    'BA', 1998, [0 80], [62.8 69.2 78.4]/100, 'range'
    'BA', 1998, [0 3], 26.6/100, 'step'
    'PB', 2005, [5 15 65], [84 96.8; 78.8 94.1; 69.9 77.7]/100, 'multibin'
    'BA', 2010, [5 10 15 80], [43.9 76.9 95.9]/100, 'step'
    'CE', 2012, [5 13], [46.6 46.6; 61.4 61.4]/100, 'multibin'
    'PE', 2015, [21.5-3, 21.5+3], [87.1 87.1]/100, 'step'
    'BA', 2015, [0 5 10 15 20 30 40 50], [57 67 77 85 90 92 95]/100, 'step'
    'RN', 2016, [0 17 48 80], [65.7 99 100]/100, 'step'
    'BA', 2017, [21 33], [87.44 87.44]/100, 'step'
    'PE', 2017, [5 15 25 35 45 55 66], [52.9 89.3 92.9 96.0 96.4 94.8]/100, 'step'
    'CE', 2018, [4 20 46 60 82], [33/44 161/177 68/69 111/114], 'step'
    'CE', 2019, [2 5 12], [3.4 32.0; 10.9 43.5]/100, 'multibin'};
centralwest = {
    'GO', 2002, [15 50], 37.3/100, 'step'
    'GO', 2009, [25.8-6.4 25.8+6.4], [49.4 49.4; 58.3 58.3]/100, 'multibin'
    'DF', 2020, [7 18 41 61 87], [40.3 57.9 59.9 71.1; 76.2 65.6 68.2 83.4]/100, 'multibin'
    'MS', 2024, [4 9 13 16], [13.4 22.6 27.7]/100, 'step'};
southeast = {
    'SP', 1998, [0 5 10 15 20 25 30 40 50 80], [19.63 15.48 11.61 15.93 17.01 15.30 14.62 11.23 15.59 14.80]/100, 'step'
    'MG', 2000, [0 16 60 80], [26 18 34.3; 17.02 26.7 18.5]/100, 'multibin'
    'RJ', 2006, [1 5 10 20], [50 71.4 83.3 20; 13 38.1 75 71; 38.9 41.5 61.8 82.1]/100, 'multibin'
    'RJ', 2014, [0 17 48 80], [45 96.1 98.2]/100, 'step'
    'SP', 2015, [10 20 40 60 92], [63.1 73.7 74.7 85.9]/100, 'step'
    'SP', 2017, [2 5 10 14 17], [33.1 33.1 33.1 33.1]/100, 'step'};

panels = {north, [0 0 5 8], [2 1], 'north'
          northeast, [0 0 15 20], [5 3], 'northeast'
          centralwest, [0 0 5 16], [4 1], 'centralwest'
          southeast, [0 0 15 8], [2 3], 'southeast'};

for f = 1:size(panels, 1)
    surveys = panels{f, 1};
    figure('Units', 'centimeters', 'Position', panels{f, 2}, 'PaperPositionMode', 'auto');
    tiledlayout(panels{f, 3}(1), panels{f, 3}(2), 'TileSpacing', 'tight')
    for k = 1:size(surveys, 1)
        nexttile(k)
        [code, yr, breaks, prev, style] = surveys{k, :};
        si = find(strcmp(cfg.state_codes, code));
        ti = max(0, yr - cfg.years(1));    % 0 = pre-2000 baseline; else end of the previous year
        draw_survey_panel(d, si, ti, cfg.age_lower, breaks, prev, style, ...
                          lo_idx, hi_idx, col_model, col_model_light, col_survey, ...
                          sprintf('%s (%d)', code, yr))
    end
    save_figure(cfg, ['figS9_serosurvey_' panels{f, 4}])
end

function draw_survey_panel(d, si, ti, ages, breaks, prev, style, lo_idx, hi_idx, ...
                           col_model, col_model_light, col_survey, title_str)
if ti == 0
    p_band = squeeze(1 - d.op_seroposall_pre_0(:, si, :));
else
    p_band = squeeze(1 - d.op_seroposall_0(:, si, ti, :));
end
p_s = sort(p_band);
fill([ages, fliplr(ages)], [p_s(hi_idx, :), fliplr(p_s(lo_idx, :))], ...
     col_model_light, 'EdgeColor', 'none', 'FaceAlpha', 0.7); hold on
plot(ages, median(p_s), '-', 'Color', col_model, 'LineWidth', 1.5)

switch style
    case 'step'
        for i = 1:numel(breaks) - 1
            plot(breaks(i:i+1), [prev(i) prev(i)], '-', 'Color', col_survey, 'LineWidth', 1.5)
        end
    case 'multibin'
        for i = 1:numel(breaks) - 1
            v = prev(:, i);
            fill([breaks(i) breaks(i+1) breaks(i+1) breaks(i)], [min(v) min(v) max(v) max(v)], ...
                 col_survey, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
            plot(breaks(i:i+1), [median(v) median(v)], ':', 'Color', col_survey, 'LineWidth', 1.5)
        end
    case 'range'
        fill([ages(1) ages(end) ages(end) ages(1)], [min(prev) min(prev) max(prev) max(prev)], ...
             col_survey, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
        plot(ages([1 end]), [median(prev) median(prev)], ':', 'Color', col_survey, 'LineWidth', 1.5)
end
xlim([0 80]); ylim([0 1])
title(title_str, 'FontSize', 9, 'FontWeight', 'normal')
set(gca, 'XTick', 0:10:80, 'FontName', 'SansSerif', 'FontSize', 12)
end
