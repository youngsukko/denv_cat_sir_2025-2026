function p = sir_base_parameters(cfg, sir_in, ve_name)
%SIR_BASE_PARAMETERS Parameters of the SIR-type model (Table S8, S9, S11).
%   sir_in : contents of results/sir_inputs.mat (see M05_r0_projection.m)
%   ve_name: 'shen' or 'daniels'

lifeexp = load(fullfile(cfg.data_dir, 'data_lifeexp.mat'), 'data_lifeexp');      % IBGE, ages 0..80
surv    = load(fullfile(cfg.data_dir, 'data_survrate.mat'), 'survrate');         % IBGE, ages 0..84
birth   = load(fullfile(cfg.data_dir, 'data_birthrate.mat'), 'data_birthper1K_since1960');  % World Bank, 1960-2022

survrate = surv.survrate(1:5:end);

p.agenum   = cfg.A5;
p.beta0    = sir_in.kappa_med(:);     % age-specific FOI multiplier (catalytic model)
p.cageadj  = sir_in.omega_med(:);     % age-specific reporting multiplier
p.delta    = sir_in.delta_med(:);     % apparent illness ratio by infection order
p.gamma    = 52;                      % recovery rate, 1-week infectious period
p.ohm      = 1;                       % waning of cross-protection, 1 year
p.mu       = 1/5;                     % ageing between five-year groups
p.mu_final = 1/8.9;                   % death rate in the 80+ group
p.q        = survrate(1:16);          % survival into the next age group
p.birthrate = birth.data_birthper1K_since1960(end);   % per 1,000, latest year
p.mu_star  = 1 / cfg.life_expectancy;

% vaccine-derived protection: <1 y -> 1-5 y -> >5 y (no return to unvaccinated)
p.ohm_vac  = 1;
p.ohm_vac2 = 1/4;
p.ohm_vac3 = 0;

% case hospitalization and case fatality ratios by age (Table S11)
p.h = [0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
       0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
p.f = [0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
       0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';

% DALY components: years lived with disability per case / hospitalization,
% years of life lost per death by age group
p.yld_case    = 0.00084;
p.yld_hosp    = 0.0051;
p.lifeexp_age = lifeexp.data_lifeexp(1:5:end)';

p.n_years     = cfg.n_sim_years;
p.vacc_groups = struct('idx', {}, 'dose', {});
p = apply_ve_settings(p, ve_settings(ve_name));
end

function p = apply_ve_settings(p, ve)
% rows: unvaccinated, <1 y, 1-5 y, >5 y after vaccination; columns: seronegative, seropositive
p.ve_inf  = [0 0; ve.inf];
p.ve_dis  = [0 0; ve.dis];
p.ve_hsp  = [0 0; ve.hsp];
p.ve_mort = [0 0; ve.mort];
p.ve_name = ve.name;
end
