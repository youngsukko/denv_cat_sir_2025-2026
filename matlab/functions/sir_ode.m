function dydt = sir_ode(t, y, p)
%SIR_ODE Age-structured SIR-type dengue model with vaccination (time in years).
%   The state vector holds four vaccination strata (unvaccinated, <1 year,
%   1-5 years and >5 years after vaccination). Each stratum stores the
%   blocks [S0 I1 R1 S1 I2 R2 S2 I3 R3 S3 I4 R4], one value per age group:
%   Sk are susceptible with k prior infections, Ik and Rk are infectious and
%   recovered after the k-th infection. The last three blocks accumulate
%   reported cases, hospitalizations and deaths (51 x A elements in total).
%
%   p: parameter struct from sir_base_parameters + apply_ve_settings, with
%      p.beta1       transmission rate for each simulated year
%      p.vacc_groups struct array (idx: age groups, dose: doses per year)

y(y < 0) = 0;
A  = p.agenum;
X  = reshape(y(1:48*A), A, 12, 4);   % [age x block x stratum]
iS = [1 4 7 10];                     % S0..S3
iI = [2 5 8 11];                     % I1..I4
iR = [3 6 9 12];                     % R1..R4
S  = X(:, iS, :);
I  = X(:, iI, :);
R  = X(:, iR, :);

%% Force of infection
N     = sum(X(:));
beta1 = p.beta1(min(floor(t) + 1, numel(p.beta1)));
foi   = beta1 * p.beta0 * (sum(p.beta0 .* sum(sum(I, 3), 2)) / N);   % [A x 1]

susc = [1, 3/4, 2/4, 1/4];   % (1 - k/4): serotypes not yet encountered, k = 0..3
sero = [1 2 2 2];            % serostatus by k: 1 seronegative, 2 seropositive
inf_mult = 1 - p.ve_inf(:, sero)';   % [k x stratum]

new_inf = zeros(A, 4, 4);            % S_k -> I_(k+1), by stratum
for v = 1:4
    new_inf(:, :, v) = foi .* (susc .* inf_mult(:, v)') .* S(:, :, v);
end

%% Infection, recovery and loss of temporary cross-protection
dS = -new_inf;
dI = new_inf - p.gamma * I;
dR = p.gamma * I;
dR(:, 1:3, :) = dR(:, 1:3, :) - p.ohm * R(:, 1:3, :);
dS(:, 2:4, 1) = dS(:, 2:4, 1) + p.ohm * R(:, 1:3, 1);
% breakthrough infections reset vaccine-derived protection to its maximum
dS(:, 2:4, 2) = dS(:, 2:4, 2) + p.ohm * sum(R(:, 1:3, 2:4), 3);

%% Vaccination of unvaccinated S and R in each target group
for g = 1:numel(p.vacc_groups)
    idx  = p.vacc_groups(g).idx;
    dose = p.vacc_groups(g).dose;
    if dose <= 0; continue; end
    pool = sum(sum(S(idx, :, 1))) + sum(sum(R(idx, :, 1)));
    rate = min(1, dose / max(pool, eps));
    dS(idx, :, 1) = dS(idx, :, 1) - rate * S(idx, :, 1);
    dS(idx, :, 2) = dS(idx, :, 2) + rate * S(idx, :, 1);
    dR(idx, :, 1) = dR(idx, :, 1) - rate * R(idx, :, 1);
    dR(idx, :, 2) = dR(idx, :, 2) + rate * R(idx, :, 1);
end

%% Waning of vaccine-derived protection (susceptible classes)
dS(:, :, 2) = dS(:, :, 2) - p.ohm_vac  * S(:, :, 2);
dS(:, :, 3) = dS(:, :, 3) + p.ohm_vac  * S(:, :, 2) - p.ohm_vac2 * S(:, :, 3);
dS(:, :, 4) = dS(:, :, 4) + p.ohm_vac2 * S(:, :, 3) - p.ohm_vac3 * S(:, :, 4);
dS(:, :, 1) = dS(:, :, 1) + p.ohm_vac3 * S(:, :, 4);

%% Demography (infectious individuals are not aged)
dS = dS + aging(S, p);
dR = dR + aging(R, p);
dS(1, 1, 1) = dS(1, 1, 1) + p.birthrate * N / 1000;

%% Burden of disease
dis_mult  = 1 - p.ve_dis(:, sero)';
hsp_mult  = 1 - p.ve_hsp(:, sero)';
mort_mult = 1 - p.ve_mort(:, sero)';
cum_case  = zeros(A, 1);
cum_hosp  = zeros(A, 1);
cum_death = zeros(A, 1);
for v = 1:4
    apparent  = p.cageadj .* new_inf(:, :, v) .* p.delta(:)';   % [A x k]
    cum_case  = cum_case  + apparent * dis_mult(:, v);
    cum_hosp  = cum_hosp  + (p.h .* apparent) * hsp_mult(:, v);
    cum_death = cum_death + (p.f .* apparent) * mort_mult(:, v);
end

dX = zeros(A, 12, 4);
dX(:, iS, :) = dS;
dX(:, iI, :) = dI;
dX(:, iR, :) = dR;
dydt = [dX(:); cum_case; cum_hosp; cum_death];
end

function d = aging(X, p)
% survivors move to the next age group; the last group leaves at mu_final
A = size(X, 1);
d = zeros(size(X));
d(2:A, :, :)   = p.mu * p.q .* X(1:A-1, :, :);
d(1:A-1, :, :) = d(1:A-1, :, :) - p.mu * X(1:A-1, :, :);
d(A, :, :)     = d(A, :, :) - p.mu_final * X(A, :, :);
end
