function obs = load_observed_data(cfg)
%LOAD_OBSERVED_DATA Reported cases and population from data/processed.
%   obs.cases, obs.pop : [S*T x A5]; row (s-1)*T + t is state
%                        cfg.state_codes{s} in year cfg.years(t), columns
%                        follow cfg.age_labels.
%   obs.cases_st, obs.pop_st : [T x S] totals over age groups.
%   obs.pop_2025 : [S x A5] population in the last data year.

obs.cases = read_long_table(cfg.file_cases, 'cases', cfg);
obs.pop   = read_long_table(cfg.file_pop,   'pop',   cfg);

obs.cases_st = zeros(cfg.T, cfg.S);
obs.pop_st   = zeros(cfg.T, cfg.S);
obs.pop_2025 = zeros(cfg.S, cfg.A5);
for s = 1:cfg.S
    rows = (s-1)*cfg.T + (1:cfg.T);
    obs.cases_st(:, s) = sum(obs.cases(rows, :), 2);
    obs.pop_st(:, s)   = sum(obs.pop(rows, :), 2);
    obs.pop_2025(s, :) = obs.pop(rows(end), :);
end
end

function M = read_long_table(file, value_col, cfg)
if ~isfile(file)
    error('Missing %s. Run the R pipeline (R/02, R/03) first.', file);
end
opts = detectImportOptions(file);
opts = setvartype(opts, {'uf', 'age_group'}, 'char');
tbl  = readtable(file, opts);

[~, s] = ismember(tbl.uf, cfg.state_codes);
[~, a] = ismember(tbl.age_group, cfg.age_labels_full);
t = tbl.year - cfg.years(1) + 1;
keep = t >= 1 & t <= cfg.T;
assert(all(s > 0) && all(a > 0), 'Unexpected state or age group in %s', file);
assert(nnz(keep) == cfg.S * cfg.T * cfg.A5, 'Incomplete state x year x age grid in %s', file);

M = zeros(cfg.S * cfg.T, cfg.A5);
M(sub2ind(size(M), (s(keep)-1)*cfg.T + t(keep), a(keep))) = tbl.(value_col)(keep);
end
