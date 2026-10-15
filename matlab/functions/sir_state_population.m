function [pop_total, prop_le1] = sir_state_population(sir_in)
%SIR_STATE_POPULATION Population at the start of 2026 from the immune profile.
%   pop_total: [S x 1] total over age and infection history
%   prop_le1 : [S x 1] proportion with zero or one prior infection

S = size(sir_in.immune_prop, 1);
pop_total = zeros(S, 1);
prop_le1  = zeros(S, 1);
for si = 1:S
    Xk = state_immune_counts(sir_in, si);
    pop_total(si) = sum(Xk, 'all');
    prop_le1(si)  = sum(Xk(:, 1:2), 'all') / pop_total(si);
end
end
