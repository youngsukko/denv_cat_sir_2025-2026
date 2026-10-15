function inputs = sir_base_parameters(cfg, sir_in, ve_name)
%SIR_BASE_PARAMETERS inputs struct of func_ode_alt_v5 / func_ode_alt_v6
%   (Tables S8, S9, S11), without the scenario-specific fields beta1,
%   vacc/vacorder (v5) or vacc_targets (v6).
%   sir_in : contents of results/sir_inputs.mat (M05_r0_projection.m)
%   ve_name: 'shen' or 'daniels'

lifeexp = load(fullfile(cfg.data_dir, 'data_lifeexp.mat'), 'data_lifeexp');
surv    = load(fullfile(cfg.data_dir, 'data_survrate.mat'), 'survrate');
birth   = load(fullfile(cfg.data_dir, 'data_birthrate.mat'), 'data_birthper1K_since1960');
data_birthper1K_since1960 = birth.data_birthper1K_since1960;
data_birthper1K_since1960(end:end+50) = data_birthper1K_since1960(end);   % hold the last year
survrate = surv.survrate(1:5:end)';

inputs = struct();
inputs.agenum    = cfg.A5;
inputs.beta0     = sir_in.kappa_med(:);   % age-specific FOI multiplier (catalytic model)
inputs.cageadj   = sir_in.omega_med(:);   % age-specific reporting multiplier
inputs.cprop_1   = sir_in.delta_med(1);   % apparent illness ratio, primary ... quaternary
inputs.cprop_2   = sir_in.delta_med(2);
inputs.cprop_3   = sir_in.delta_med(3);
inputs.cprop_4   = sir_in.delta_med(4);
inputs.gamma     = 1/(1/52);              % recovery, 1-week infectious period
inputs.ohm       = 1/1;                   % waning of cross-protection, 1 year
inputs.i1 = 3/4; inputs.i2 = 2/4; inputs.i3 = 1/4;
inputs.mu        = 1/(5*1);               % ageing between five-year groups
inputs.mu_final  = 1/(8.9*1);             % death rate of the 80+ group
inputs.q         = survrate(1:16)';
inputs.birthrate = data_birthper1K_since1960(64);
inputs.ohm_vac   = 1/(1);                 % <1 y -> 1-5 y after vaccination
inputs.ohm_vac2  = 1/(4);                 % 1-5 y -> >5 y
inputs.ohm_vac3  = 0;                     % no return to unvaccinated

% case hospitalization and case fatality ratios by age (Table S11)
h_ref = [0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
         0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
f_ref = [0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
         0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';
inputs.h1 = h_ref; inputs.h2 = h_ref; inputs.h3 = h_ref; inputs.h4 = h_ref;
inputs.f1 = f_ref; inputs.f2 = f_ref; inputs.f3 = f_ref; inputs.f4 = f_ref;

% vaccine effectiveness: e_* (<1 y), ew_* (1-5 y), eww_* (>5 y); n/p = sero-negative/positive;
% i/d/h/m = infection/disease/hospitalization/death (death = hospitalization)
ve = ve_settings(ve_name);
pre  = {'e', 'ew', 'eww'};
sero = {'n', 'p'};
for st = 1:3
    for sg = 1:2
        tag = sprintf('%s_%s', pre{st}, sero{sg});
        inputs.([tag '_i']) = ve.inf(st, sg);
        inputs.([tag '_d']) = ve.dis(st, sg);
        inputs.([tag '_h']) = ve.hsp(st, sg);
        inputs.([tag '_m']) = ve.hsp(st, sg);
    end
end

% not used by the ODE: projection horizon and DALY weights
inputs.n_years     = cfg.n_sim_years;
inputs.lifeexp_age = lifeexp.data_lifeexp(1:5:end)';
end
