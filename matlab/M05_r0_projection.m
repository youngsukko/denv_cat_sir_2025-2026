%% Fig. S3: reproduction number under the endemic equilibrium and its
%  projection with Gaussian process regression; saves the inputs of the
%  SIR-type model.
%
%  R0 = lambda / mu* + 1 with lambda the per-serotype FOI (rescaled to the
%  0-4 age group) and 1/mu* the life expectancy. The posterior quartiles of
%  R0 in 2000-2025 are extrapolated separately for cfg.n_proj_years years.
%
%  Input : data/processed (R pipeline), results_big/catalytic_*_draws.mat
%  Output: output/figures/figS3_R0_projection.jpg
%          results/sir_inputs.mat (small; used by M06-M12)

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
obs = load_observed_data(cfg);
d   = load_catalytic_draws(cfg, 'all');

S = cfg.S; T = cfg.T; A5 = cfg.A5;
n_fut   = cfg.n_proj_years;
t_hist  = (1:T)';
t_fut   = (T+1 : T+n_fut)';
yr_all  = [cfg.years, cfg.years(end) + (1:n_fut)];
col_foi = [0.839 0.153 0.157];
nq      = numel(cfg.R0_quantiles);

%% R0 quartiles by state and year, with GPR projection
R0_proj = zeros(S, T + n_fut, nq);
for si = 1:S
    R_draws = squeeze(d.lambda_draws(:, si, :)) * cfg.life_expectancy + 1;   % [n x T]
    for qi = 1:nq
        r_hist = quantile(R_draws, cfg.R0_quantiles(qi), 1)';
        mdl    = fitrgp(t_hist, r_hist, 'Basis', 'linear', 'FitMethod', 'exact', 'PredictMethod', 'exact');
        R0_proj(si, :, qi) = [r_hist; max(1, predict(mdl, t_fut))]';
    end
end

%% Fig. S3
panel_map_layout(cfg, 18, 22);
for si = 1:S
    nexttile(cfg.state_tile(si))
    q25 = R0_proj(si, :, 1); q50 = R0_proj(si, :, 2); q75 = R0_proj(si, :, 3);
    h = 1:T; f = T:T+n_fut;
    fill([yr_all(h), fliplr(yr_all(h))], [q75(h), fliplr(q25(h))], col_foi, 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on
    fill([yr_all(f), fliplr(yr_all(f))], [q75(f), fliplr(q25(f))], col_foi, 'EdgeColor', 'none', 'FaceAlpha', 0.15)
    plot(yr_all(h), q50(h), '-', 'Color', col_foi, 'LineWidth', 1.5)
    plot(yr_all(f), q50(f), ':', 'Color', col_foi, 'LineWidth', 1.5)
    xline(cfg.years(end) + 0.5, 'k:', 'LineWidth', 0.8)
    xlim([yr_all(1), yr_all(end)]); ylim([1 inf])
    xticks([2000 2010 2020 2030]); xticklabels({'00', '10', '20', '30'})
    title(cfg.state_names{si}, 'FontSize', 8, 'FontWeight', 'normal')
    set(gca, 'FontName', 'SansSerif', 'FontSize', 8)
end
nexttile(cfg.legend_tile); axis off; hold on
h1 = plot(nan, nan, '-', 'Color', col_foi, 'LineWidth', 1.5);
h2 = plot(nan, nan, '--', 'Color', col_foi, 'LineWidth', 1.5);
h3 = fill(nan, nan, col_foi, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
legend([h1 h2 h3], {'Median R_0 (estimated)', 'Median R_0 (projected)', 'IQR'}, ...
       'Location', 'northwest', 'FontSize', 11); legend('boxoff')
save_figure(cfg, 'figS3_R0_projection')

%% Inputs of the SIR-type model (posterior medians, start of 2026)
immune_prop = zeros(S, A5, 5);   % proportion with 0..4 prior infections, end of 2025
for k = 0:4
    immune_prop(:, :, k+1) = squeeze(median(d.(sprintf('op_seroposall_%d', k))(:, :, T, :), 1));
end
pop_2025    = obs.pop_2025;
kappa_med   = median(d.kappa_draws, 1)';
omega_med   = median(d.omega_draws, 1)';
delta_med   = median(d.delta_draws, 1)';
state_codes = cfg.state_codes;
years_proj  = yr_all;
R0_quantiles = cfg.R0_quantiles;

save(cfg.file_sir_inputs, 'R0_proj', 'R0_quantiles', 'years_proj', 'T', 'immune_prop', ...
     'pop_2025', 'kappa_med', 'omega_med', 'delta_med', 'state_codes')
fprintf('Saved: %s\n', cfg.file_sir_inputs);
