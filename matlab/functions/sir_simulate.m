function out = sir_simulate(inputs, sir_in, si, rc, vacc_input)
%SIR_SIMULATE Simulate one state over 2026-2030 and return the burden by age.
%   inputs    : sir_base_parameters
%   sir_in    : contents of results/sir_inputs.mat
%   si        : state index (cfg.state_codes order)
%   rc        : R0 scenario (1 = Q25, 2 = Q50, 3 = Q75)
%   vacc_input: struct with fields vacc (doses per year) and vacorder
%               (target ages, func_ode_alt_v5), or a struct array with fields
%               idx and dose (one element per target group, func_ode_alt_v6).
%               Omit for no vaccination.
%
%   out.case, out.hosp, out.dead, out.daly: [1 x A5] totals over the horizon

if nargin < 5; vacc_input = struct('vacc', 0, 'vacorder', 1); end

A5 = inputs.agenum;
T  = sir_in.T;
R0pred = sir_in.R0_proj(:, :, rc);
Xk = state_immune_counts(sir_in, si);

% transmission rates from the projected R0 (beta = R0 * alpha); the ODE uses
% beta1(floor(t) + 1), so all projected years are passed
inputs.beta1 = R0pred(si, T+1:end) * 52;
y0 = sir_initial_state(Xk(:,1), Xk(:,2), Xk(:,3), Xk(:,4), Xk(:,5), R0pred(si, T), A5);
tspan = 0:inputs.n_years;

if isfield(vacc_input, 'vacorder')
    inputs.vacc     = vacc_input.vacc;
    inputs.vacorder = vacc_input.vacorder;
    [~, y] = ode23(@(t, yy) func_ode_alt_v5(t, yy, inputs), tspan, y0);
else
    inputs.vacc_targets = vacc_input;
    [~, y] = ode23(@(t, yy) func_ode_alt_v6(t, yy, inputs), tspan, y0);
end

Casen = y(:, A5*48+(1:A5));
Hosp  = y(:, A5*49+(1:A5));
Dead  = y(:, A5*50+(1:A5));
YLD_age = Casen*0.00084 + Hosp*0.0051;
YLL_age = Dead .* inputs.lifeexp_age;

out.case = Casen(end, :);
out.hosp = Hosp(end, :);
out.dead = Dead(end, :);
out.daly = YLD_age(end, :) + YLL_age(end, :);
end
