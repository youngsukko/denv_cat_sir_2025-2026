function Xk = seropos_draws(d, si, ti)
%SEROPOS_DRAWS Immune profile draws of one state as a 1x5 cell of [n_draw x A5].
%   ti = 0 gives the start of 2000 (pre-2000 baseline); ti = 1..T gives the
%   end of year cfg.years(ti).

Xk = cell(1, 5);
for k = 0:4
    if ti == 0
        Xk{k+1} = squeeze(d.(sprintf('op_seroposall_pre_%d', k))(:, si, :));
    else
        Xk{k+1} = squeeze(d.(sprintf('op_seroposall_%d', k))(:, si, ti, :));
    end
end
end
