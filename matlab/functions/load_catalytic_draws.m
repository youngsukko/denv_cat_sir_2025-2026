function d = load_catalytic_draws(cfg, part)
%LOAD_CATALYTIC_DRAWS Posterior draws of the catalytic model (results_big).
%   part = 'param' : lambda, lambda_pre, kappa, omega, delta, rho draws
%   part = 'ppc'   : immune profiles, cumulative predicted cases, infections
%   part = 'all'   : both
%
%   Array layouts (n = number of draws):
%     lambda_draws, rho_draws    [n x S x T]   (lambda rescaled by kappa(0-4))
%     lambda_pre_draws           [n x S]
%     kappa_draws, omega_draws   [n x A5]      (kappa relative to 0-4 years)
%     delta_draws                [n x 4]
%     op_seroposall_k            [n x S x T x A5]  proportion with k prior
%                                infections at the end of year t, k = 0..4
%     op_seroposall_pre_k        [n x S x A5]  same at the start of 2000
%     inf_total                  [n x S x T x A5]  expected new infections
%     cum_cases_age              [n x S x A5]  predicted cases, 2000-2025

if nargin < 2; part = 'all'; end
d = struct();

if any(strcmp(part, {'param', 'all'}))
    check_file(cfg.file_param);
    d = merge(d, load(cfg.file_param, 'lambda_draws', 'lambda_pre_draws', 'kappa_draws', ...
                      'omega_draws', 'delta_draws', 'rho_draws'));
end
if any(strcmp(part, {'ppc', 'all'}))
    check_file(cfg.file_ppc);
    vars = [arrayfun(@(k) sprintf('op_seroposall_%d', k), 0:4, 'UniformOutput', false), ...
            arrayfun(@(k) sprintf('op_seroposall_pre_%d', k), 0:4, 'UniformOutput', false), ...
            {'inf_total', 'cum_cases_age'}];
    d = merge(d, load(cfg.file_ppc, vars{:}));
end
end

function check_file(f)
if ~isfile(f)
    error(['Missing %s.\nRun R/04_fit_catalytic_model.R and R/05_export_posterior.R, ' ...
           'or copy the existing posterior files into results_big/ (see README).'], f);
end
end

function a = merge(a, b)
fn = fieldnames(b);
for k = 1:numel(fn); a.(fn{k}) = b.(fn{k}); end
end
