function ratio = case_infection_ratio(d)
%CASE_INFECTION_RATIO Predicted reported cases per infection, 2000-2025.
%   ratio: [n_draw x S], cumulative over years and age groups.

infections = sum(sum(d.inf_total, 3), 4);   % [n x S]
ratio = sum(d.cum_cases_age, 3) ./ infections;
end
