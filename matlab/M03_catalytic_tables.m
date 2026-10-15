%% Catalytic model tables: Tables S1-S7 (CSV, states in regional order).
%
%  Input : data/processed (R pipeline), results_big/catalytic_*_draws.mat
%  Output: output/tables/tableS1_kappa.csv ... tableS7_immune_profile.csv

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
obs = load_observed_data(cfg);
d   = load_catalytic_draws(cfg, 'all');

S = cfg.S; T = cfg.T; A5 = cfg.A5;
codes = cfg.state_codes(cfg.reg_order);

fmt2   = '%.2f [%.2f, %.2f]';
fmt_nl = '%.2f\n[%.0f, %.0f]';
write_tbl = @(c, vars, rows, name) writetable(cell2table(c, 'VariableNames', vars, 'RowNames', rows), ...
    fullfile(cfg.table_dir, [name '.csv']), 'WriteRowNames', true);

%% Table S1: age-specific FOI multiplier relative to 0-4 years (95% CrI)
c = arrayfun(@(a) format_summary(d.kappa_draws(:, a), 1, fmt2), (1:A5)', 'UniformOutput', false);
write_tbl(c, {'kappa'}, cfg.age_labels_full, 'tableS1_kappa')

%% Table S2: FOI (4 x per-serotype FOI) x 100, pre-2000 and 2000-2025 (IQR)
periods = ['pre00', arrayfun(@(y) sprintf('y%02d', mod(y, 100)), cfg.years, 'UniformOutput', false)];
c = cell(S, T + 1);
for k = 1:S
    si = cfg.reg_order(k);
    c{k, 1} = format_summary(4 * d.lambda_pre_draws(:, si), 100, fmt_nl, 'iqr');
    for ti = 1:T
        c{k, ti+1} = format_summary(4 * d.lambda_draws(:, si, ti), 100, fmt_nl, 'iqr');
    end
end
write_tbl(c, periods, codes, 'tableS2_foi')

%% Table S3: infection-order-dependent apparent illness ratio (95% CrI)
c = arrayfun(@(k) format_summary(d.delta_draws(:, k), 1, fmt2), (1:4)', 'UniformOutput', false);
write_tbl(c, {'delta'}, {'Primary', 'Secondary', 'Tertiary', 'Quaternary'}, 'tableS3_delta')

%% Table S4: age-specific reporting multiplier (95% CrI)
c = arrayfun(@(a) format_summary(d.omega_draws(:, a), 1, fmt2), (1:A5)', 'UniformOutput', false);
write_tbl(c, {'omega'}, cfg.age_labels_full, 'tableS4_omega')

%% Table S5: apparent illness-to-case report ratio (%) by state and year (95% CrI)
c = cell(S, T);
for k = 1:S
    for ti = 1:T
        c{k, ti} = format_summary(d.rho_draws(:, cfg.reg_order(k), ti), 100, fmt_nl);
    end
end
write_tbl(c, periods(2:end), codes, 'tableS5_report_ratio')

%% Table S6: case/infection ratio, 2000-2025 (95% CrI)
ratio = case_infection_ratio(d);
c = arrayfun(@(k) format_summary(ratio(:, cfg.reg_order(k)), 1, fmt2), (1:S)', 'UniformOutput', false);
write_tbl(c, {'case_infection_ratio'}, codes, 'tableS6_case_infection_ratio')

%% Table S7: population share (%) by number of prior infections, 2000 and 2026 (95% CrI)
fmt_pct = '%.2f\n[%.2f, %.2f]';
c = cell(S, 10);
for k = 1:S
    si = cfg.reg_order(k);
    snapshots = {0, obs.pop((si-1)*T + 1, :); T, obs.pop((si-1)*T + T, :)};   % start 2000 / end 2025
    for j = 1:2
        Xk  = seropos_draws(d, si, snapshots{j, 1});
        pop = snapshots{j, 2};
        denom = Xk{1} + Xk{2} + Xk{3} + Xk{4} + Xk{5};
        for kx = 1:5
            share = (Xk{kx} ./ denom) * pop' / sum(pop);   % [n x 1], age-weighted within the state
            c{k, 5*(j-1) + kx} = format_summary(share, 100, fmt_pct);
        end
    end
end
vars = [strcat('y2000_', {'0','1','2','3','4'}), strcat('y2026_', {'0','1','2','3','4'})];
write_tbl(c, vars, codes, 'tableS7_immune_profile')

fprintf('Saved Tables S1-S7 to %s\n', cfg.table_dir);
