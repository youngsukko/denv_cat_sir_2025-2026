function y0 = sir_initial_state(Xk, R0, p)
%SIR_INITIAL_STATE Initial condition at the start of 2026.
%   Xk: [A x 5] number of people with k = 0..4 prior infections (catalytic model).
%   Under the endemic equilibrium a fraction mu*/(mu* + alpha) * (1 - 1/R0)
%   of each class with k < 4 prior infections is infectious; recovered
%   compartments start empty and everyone else is susceptible.

A  = p.agenum;
sc = p.mu_star / (p.mu_star + p.gamma) * (1 - 1/R0);
y0 = zeros(51*A, 1);
for k = 0:3
    y0(A*(3*k)     + (1:A)) = Xk(:, k+1) * (1 - sc);   % S_k
    y0(A*(3*k + 1) + (1:A)) = Xk(:, k+1) * sc;         % I_(k+1)
end
y0(A*11 + (1:A)) = Xk(:, 5);                           % R4: immune to all serotypes
end
