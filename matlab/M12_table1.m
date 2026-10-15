%% Table 1: burden of disease averted over 2026-2030 under each strategy.
%  Each cell: total averted Q50 [range Q25-Q75 R0 scenarios], then
%  (per 1,000 vaccinated Q50 [range]). One table per VE setting; the main
%  text uses Shen.
%
%  Input : results/sir_inputs.mat (M05), results/opt_allocation_*.mat (M09),
%          data/processed (R pipeline)
%  Output: output/tables/table1_bod_averted_{shen,daniels}.csv
%          output/sim/sim_table1.mat

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
obs    = load_observed_data(cfg);
sir_in = load(cfg.file_sir_inputs);

S  = cfg.S;
ny = cfg.n_sim_years;
doses_brazil = brazil_dose_allocation(cfg);
mass_pop  = sum(obs.pop_2025(:, cfg.target_ages.mass), 2);
mass_frac = mass_pop / sum(mass_pop);   % mass vaccination: same coverage of ages 5-59 in every state

strategy_names = {'Children, current distribution (baseline)', 'Children, DALY-optimized', ...
                  'Elderly, DALY-optimized', 'Mass vaccination 10M', 'Mass vaccination 20M', ...
                  'Mass vaccination 30M'};
strategy_vacorder = [cfg.vacorder.children, cfg.vacorder.children, cfg.vacorder.elderly, ...
                     cfg.vacorder.mass, cfg.vacorder.mass, cfg.vacorder.mass];
annual_total   = [cfg.annual_doses, cfg.annual_doses, cfg.annual_doses, cfg.mass_doses];
n_st = numel(strategy_names);
bod_names = {'case', 'hosp', 'dead', 'daly'};

%% No vaccination (independent of VE)
p0 = sir_base_parameters(cfg, sir_in, 'shen');
base = zeros(3, 4);   % national totals [R0 x outcome]
for rc = 1:3
    for si = 1:S
        out = sir_simulate(p0, sir_in, si, rc);
        base(rc, :) = base(rc, :) + cellfun(@(f) sum(out.(f)), bod_names);
    end
end

%% Strategies for each VE setting
averted = zeros(n_st, 3, 4, 2);   % [strategy x R0 x outcome x VE]
for vi = 1:2
    ve  = cfg.ve_names{vi};
    p   = sir_base_parameters(cfg, sir_in, ve);
    opt = load(cfg.file_opt(ve));
    for st = 1:n_st
        for rc = 1:3
            switch st
                case 1; doses = doses_brazil;
                case 2; doses = opt.doses_opt(:, rc, 1);
                case 3; doses = opt.doses_opt(:, rc, 2);
                otherwise; doses = annual_total(st) * mass_frac;
            end
            fprintf('[%s] %s, R0 %s\n', cfg.ve_labels{vi}, strategy_names{st}, cfg.R0_names{rc});
            total = zeros(1, 4);
            for si = 1:S
                out = sir_simulate(p, sir_in, si, rc, struct('vacc', doses(si), 'vacorder', strategy_vacorder(st)));
                total = total + cellfun(@(f) sum(out.(f)), bod_names);
            end
            averted(st, rc, :, vi) = base(rc, :) - total;
        end
    end
end
save(fullfile(cfg.sim_dir, 'sim_table1.mat'), 'base', 'averted', 'strategy_names', 'annual_total')

%% Tables
vactot = annual_total * ny / 1000;   % thousands vaccinated over the horizon
for vi = 1:2
    c = cell(n_st, 4);
    for st = 1:n_st
        for m = 1:4
            a  = averted(st, :, m, vi);   % [Q25 Q50 Q75]
            pk = a / vactot(st);
            c{st, m} = sprintf('%s\n[%s, %s]\n(%.2f\n[%.2f, %.2f])', format_thousands(a(2)), ...
                format_thousands(min(a([1 3]))), format_thousands(max(a([1 3]))), ...
                pk(2), min(pk([1 3])), max(pk([1 3])));
        end
    end
    tbl = cell2table(c, 'VariableNames', {'cases_averted', 'hospitalizations_averted', ...
                                          'deaths_averted', 'dalys_averted'}, 'RowNames', strategy_names);
    f = fullfile(cfg.table_dir, sprintf('table1_bod_averted_%s.csv', cfg.ve_names{vi}));
    writetable(tbl, f, 'WriteRowNames', true)
    fprintf('Saved: %s\n', f);
end

%% Long-format results for review, and the comparisons quoted in the text
long = {};
for vi = 1:2
    for st = 1:n_st
        for rc = 1:3
            for m = 1:4
                a   = averted(st, rc, m, vi);
                pk  = a / vactot(st);
                pk0 = averted(1, rc, m, vi) / vactot(1);   % baseline strategy, same VE and R0
                long(end+1, :) = {cfg.ve_labels{vi}, strategy_names{st}, cfg.R0_names{rc}, bod_names{m}, ...
                                  a, pk, 100 * (pk / pk0 - 1)}; %#ok<SAGROW>
            end
        end
    end
end
long_tbl = cell2table(long, 'VariableNames', {'VE', 'strategy', 'R0', 'outcome', 'averted_total', ...
                                              'averted_per1000_vaccinated', 'pct_change_vs_baseline_per1000'});
writetable(long_tbl, fullfile(cfg.table_dir, 'table1_long.csv'))

rows = {};
for st = 2:3
    for m = [4 3]   % DALYs, deaths
        sel  = strcmp(long_tbl.strategy, strategy_names{st}) & strcmp(long_tbl.outcome, bod_names{m});
        shen = sel & strcmp(long_tbl.VE, 'Shen');
        v_sh = long_tbl.pct_change_vs_baseline_per1000(shen);   % Q25, Q50, Q75
        v_all = long_tbl.pct_change_vs_baseline_per1000(sel);
        rows(end+1, :) = {sprintf('%s, %s per 1,000 vaccinated vs baseline (%%), Shen', strategy_names{st}, bod_names{m}), ...
                          v_sh(2), min(v_sh([1 3])), max(v_sh([1 3])), 'Q50 [range over R0]'}; %#ok<SAGROW>
        rows(end+1, :) = {sprintf('%s, %s per 1,000 vaccinated vs baseline (%%), Shen and Daniels', strategy_names{st}, bod_names{m}), ...
                          v_sh(2), min(v_all), max(v_all), 'Shen Q50 [range over R0 and VE]'}; %#ok<SAGROW>
    end
end
target_pop = [sum(obs.pop_2025(:, cfg.target_ages.children), 'all'), sum(obs.pop_2025(:, cfg.target_ages.elderly), 'all'), ...
              sum(obs.pop_2025(:, cfg.target_ages.mass), 'all')];
cov_cases  = {'children (5-14), 5M per year', annual_total(1), target_pop(1)
              'elderly (60+), 5M per year',   annual_total(3), target_pop(2)
              'ages 5-59, 10M per year',       annual_total(4), target_pop(3)
              'ages 5-59, 20M per year',       annual_total(5), target_pop(3)
              'ages 5-59, 30M per year',       annual_total(6), target_pop(3)};
for k = 1:size(cov_cases, 1)
    rows(end+1, :) = {sprintf('Coverage: doses per year / 2025 population, %s (%%)', cov_cases{k, 1}), ...
                      100 * cov_cases{k, 2} / cov_cases{k, 3}, NaN, NaN, ...
                      sprintf('over %d years: %.1f%%', ny, 100 * ny * cov_cases{k, 2} / cov_cases{k, 3})}; %#ok<SAGROW>
end
write_numbers(cfg, 'M12_table1', rows)
