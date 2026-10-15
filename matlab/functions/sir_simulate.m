function out = sir_simulate(p, sir_in, si, rc, vacc_groups)
%SIR_SIMULATE Simulate one state over the projection horizon (2026-2030).
%   p          : parameters (sir_base_parameters)
%   sir_in     : contents of results/sir_inputs.mat
%   si         : state index (cfg.state_codes order)
%   rc         : reproduction number scenario (1 = Q25, 2 = Q50, 3 = Q75)
%   vacc_groups: struct array with fields idx (age groups) and dose
%                (people vaccinated per year); empty for no vaccination
%
%   out.case, out.hosp, out.dead, out.daly: [1 x A] totals over the horizon

if nargin < 5; vacc_groups = struct('idx', {}, 'dose', {}); end

A    = p.agenum;
T    = sir_in.T;
R0   = sir_in.R0_proj(si, :, rc);
Xk   = state_immune_counts(sir_in, si);

p.beta1       = R0(T + (1:p.n_years)) * p.gamma;   % beta = R0 * alpha
p.vacc_groups = vacc_groups;
y0 = sir_initial_state(Xk, R0(T), p);

[~, y] = ode23(@(t, yy) sir_ode(t, yy, p), 0:p.n_years, y0);

out.case = y(end, A*48 + (1:A));
out.hosp = y(end, A*49 + (1:A));
out.dead = y(end, A*50 + (1:A));
out.daly = out.case * p.yld_case + out.hosp * p.yld_hosp + out.dead .* p.lifeexp_age;
end
