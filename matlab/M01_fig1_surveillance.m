%% Fig. 1 (surveillance panels): cumulative reported incidence map and
%  age distribution of reported cases by year, 2000-2025.
%  The workflow diagram around these panels was drawn by hand.
%
%  Input : data/processed (R pipeline), matlab/data/shapefile
%  Output: output/figures/fig1_incidence_map.png, fig1_age_distribution.png

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
obs = load_observed_data(cfg);
S = cfg.S; T = cfg.T; A5 = cfg.A5;

%% Cumulative incidence per 100 population by state (choropleth, alpha = relative incidence)
inc = sum(obs.cases_st, 1)' ./ sum(obs.pop_st, 1)' * 100;   % [S x 1]
inc_norm = inc / max(inc);

shp = load_brazil_shapes(cfg);
figure('Units', 'centimeters', 'Position', [0 0 10 12], 'PaperPositionMode', 'auto');
hold on
for i = 1:numel(shp)
    for k = 1:numel(shp(i).parts)
        xy = shp(i).parts{k};
        patch(xy(:,1), xy(:,2), 'r', 'EdgeColor', 'k', 'LineWidth', 0.8, ...
              'FaceAlpha', inc_norm(shp(i).si));
    end
end
axis equal off
save_figure(cfg, 'fig1_incidence_map')

%% Proportion of reported cases by age group, one curve per year (cyan -> magenta)
cases_nat = zeros(T, A5);
for s = 1:S
    cases_nat = cases_nat + obs.cases((s-1)*T + (1:T), :);
end

figure('Units', 'centimeters', 'Position', [0 0 9 6], 'PaperPositionMode', 'auto');
hold on
for t = 1:T
    c = (t-1) / (T-1);
    plot(cases_nat(t, :) / sum(cases_nat(t, :)), 'Color', [c, 1-c, 1], 'LineWidth', 1)
end
xticks(1:A5); xticklabels(cfg.age_labels); xlim([1 A5])
xlabel('Age group'); ylabel('Proportion by age group')
set(gca, 'FontName', 'SansSerif', 'FontSize', 11)
save_figure(cfg, 'fig1_age_distribution')
