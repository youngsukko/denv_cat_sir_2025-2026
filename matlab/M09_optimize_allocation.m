%% DALY-minimizing allocation of 5 million vaccinations per year across states,
%  for children (5-14 y) or elderly (60+ y), each R0 scenario and VE setting.
%  fmincon (interior point, BFGS Hessian, forward finite differences).
%
%  This step is slow and its result can differ slightly between runs, so the
%  allocations are saved in results/ (tracked by git) and read by M10-M12.
%  Each case is written as soon as it finishes; set the *_TO_RUN lists to
%  run a subset.
%
%  Input : results/sir_inputs.mat (M05)
%  Output: results/opt_allocation_shen.mat, results/opt_allocation_daniels.mat

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
sir_in = load(cfg.file_sir_inputs);

VE_TO_RUN      = {'shen', 'daniels'};
TARGETS_TO_RUN = 1:2;      % 1 = children (5-14 y), 2 = elderly (60+ y)
R0_TO_RUN      = 1:3;      % Q25, Q50, Q75
USE_PARALLEL   = false;    % parallel finite differences (needs Parallel Computing Toolbox)

S = cfg.S;
target_idx = {cfg.target_ages.children, cfg.target_ages.elderly};
[pop_state, ~] = sir_state_population(sir_in);
doses_annual_brazil = brazil_dose_allocation(cfg);

opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 300, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-8, ...
    'Algorithm', 'interior-point', ...
    'HessianApproximation', 'bfgs', ...
    'FiniteDifferenceType', 'forward', ...
    'UseParallel', USE_PARALLEL);

for vi = 1:numel(VE_TO_RUN)
    ve_name  = VE_TO_RUN{vi};
    out_file = cfg.file_opt(ve_name);
    p = sir_base_parameters(cfg, sir_in, ve_name);

    if isfile(out_file)
        res = load(out_file);
    else
        res.doses_opt    = nan(S, 3, 2);   % doses per year by state [S x R0 x target]
        res.daly_opt     = nan(3, 2);      % national DALYs over the horizon, optimal allocation
        res.daly_uniform = nan(3, 2);      % same with equal coverage of the target group
    end

    for tc = TARGETS_TO_RUN
        for rc = R0_TO_RUN
            fprintf('\n=== VE %s | %s | R0 %s ===\n', ve_name, cfg.target_names{tc}, cfg.R0_names{rc});
            % x: relative coverage of the state population; doses sum to the national supply
            to_doses = @(x) cfg.annual_doses * x .* pop_state / sum(x .* pop_state);
            obj = @(x) national_daly(to_doses(max(x, 0)), p, sir_in, rc, target_idx{tc});

            x0 = ones(S, 1);
            [x_opt, daly_opt] = fmincon(obj, x0, [], [], [], [], zeros(S, 1), [], [], opts);

            res.doses_opt(:, rc, tc)  = to_doses(x_opt);
            res.daly_opt(rc, tc)      = daly_opt;
            res.daly_uniform(rc, tc)  = obj(x0);
            fprintf('DALYs: uniform %.1f, optimal %.1f (%.2f%% lower)\n', res.daly_uniform(rc, tc), ...
                    daly_opt, 100 * (1 - daly_opt / res.daly_uniform(rc, tc)));

            res.pop_state           = pop_state;
            res.doses_annual_brazil = doses_annual_brazil;
            res.state_codes         = cfg.state_codes;
            res.target_names        = cfg.target_names;
            res.R0_names            = cfg.R0_names;
            res.ve_name             = ve_name;
            res.n_sim_years         = cfg.n_sim_years;
            save(out_file, '-struct', 'res')
            fprintf('Saved: %s\n', out_file);
        end
    end
end

function total = national_daly(doses, p, sir_in, rc, idx)
total = 0;
for si = 1:numel(doses)
    out   = sir_simulate(p, sir_in, si, rc, struct('idx', idx, 'dose', doses(si)));
    total = total + sum(out.daly);
end
end
