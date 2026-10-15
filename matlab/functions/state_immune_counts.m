function Xk = state_immune_counts(sir_in, si)
%STATE_IMMUNE_COUNTS [A x 5] people with 0..4 prior infections, start of 2026.
%   Posterior median proportions (catalytic model) times the 2025 population.

Xk = squeeze(sir_in.immune_prop(si, :, :)) .* sir_in.pop_2025(si, :)';
end
