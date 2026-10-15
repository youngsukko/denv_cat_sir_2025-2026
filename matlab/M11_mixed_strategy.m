%% Fig. 5: burden-of-disease averting ratios under linear combinations of the
%  children- and elderly-targeting optimized allocations (Shen VE), and the
%  nationwide ratios quoted in the Results (100% and 60% to children).
%
%  Input : results/sir_inputs.mat (M05), results/opt_allocation_shen.mat (M09)
%  Output: output/figures/fig5_mixed_strategy.png (R0 Q50)
%          output/sim/sim_mixed_strategy.mat

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
sir_in = load(cfg.file_sir_inputs);

VE_NAME = 'shen';
FIG_R0  = 2;            % R0 scenario shown in Fig. 5 (Q50)
w_grid  = 0:0.05:1;     % share of the supply given to the children-targeting allocation

S = cfg.S;
p   = sir_base_parameters(cfg, sir_in, VE_NAME);
opt = load(cfg.file_opt(VE_NAME));
nW  = numel(w_grid);
bod_names = {'case', 'hosp', 'dead', 'daly'};

%% Simulations: no vaccine and each mix, totals over age by state
base = zeros(S, 3, 4);        % [state x R0 x outcome]
mix  = zeros(S, 3, nW, 4);    % [state x R0 x w x outcome]
for rc = 1:3
    for si = 1:S
        out = sir_simulate(p, sir_in, si, rc);
        base(si, rc, :) = cellfun(@(f) sum(out.(f)), bod_names);
    end
    for wi = 1:nW
        w = w_grid(wi);
        fprintf('R0 %s, children share %.2f\n', cfg.R0_names{rc}, w);
        for si = 1:S
            groups = struct('idx',  {cfg.target_ages.children, cfg.target_ages.elderly}, ...
                            'dose', {w * opt.doses_opt(si, rc, 1), (1 - w) * opt.doses_opt(si, rc, 2)});
            out = sir_simulate(p, sir_in, si, rc, groups);
            mix(si, rc, wi, :) = cellfun(@(f) sum(out.(f)), bod_names);
        end
    end
end
save(fullfile(cfg.sim_dir, 'sim_mixed_strategy.mat'), 'base', 'mix', 'w_grid')

%% Fig. 5
ylab = {'Case', 'Hospitalization', 'Fatality', 'DALY'};
x_axis = 1 - w_grid;   % 0 = children-targeting only, 1 = elderly-targeting only
figure('Units', 'centimeters', 'Position', [0 0 18 12], 'PaperPositionMode', 'auto');
tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'loose')
for m = 1:4
    nexttile(m); hold on
    h = gobjects(6, 1);
    for r = 1:5
        a = averting_ratio(base, mix, cfg.region_of == r, m);
        h(r) = plot(x_axis, a(FIG_R0, :), '-', 'Color', cfg.region_colors(r, :), 'LineWidth', 1.5);
    end
    a = averting_ratio(base, mix, true(S, 1), m);
    h(6) = plot(x_axis, a(FIG_R0, :), 'k--', 'LineWidth', 1.5);
    ylabel(sprintf('%s averting ratio (%%)', ylab{m}))
    xlim([0 1]); xticks([0 0.5 1]); xticklabels({'Only\newlinechildren', '50%', 'Only\newlineelderly'})
    set(gca, 'FontName', 'SansSerif', 'FontSize', 10)
    if m == 1
        lg = legend(h, [cfg.region_short, {'Nationwide'}], 'Location', 'best', 'FontSize', 8);
        legend('boxoff'); lg.NumColumns = 3;
    end
end
save_figure(cfg, 'fig5_mixed_strategy')

%% Results text: nationwide averting ratios, Q50 [range over Q25-Q75]
fprintf('\nNationwide averting ratio (%%), Q50 [range over R0 Q25-Q75]\n');
for w_target = [1.0 0.6]
    [~, wi] = min(abs(w_grid - w_target));
    fprintf('Children share %.0f%%\n', 100 * w_grid(wi));
    for m = 1:4
        a = averting_ratio(base, mix, true(S, 1), m);
        fprintf('  %-16s %.2f [%.2f, %.2f]\n', ylab{m}, a(2, wi), min(a([1 3], wi)), max(a([1 3], wi)));
    end
end

function a = averting_ratio(base, mix, idx, m)
% percentage of outcome m averted in the states idx: [3 (R0) x numel(w_grid)]
b  = sum(base(idx, :, m), 1)';                    % [3 x 1]
mx = reshape(sum(mix(idx, :, :, m), 1), 3, []);   % [3 x nW]
a  = (b - mx) ./ b * 100;
end
