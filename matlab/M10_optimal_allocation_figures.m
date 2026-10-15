%% Optimized state allocations: Figs. S14-S17 (scatter vs past incidence and
%  choropleth maps, Shen and Daniels VE), Fig. S18 (Daniels vs Shen) and the
%  correlations reported in Appendix S5-3.
%
%  Input : results/opt_allocation_{shen,daniels}.mat (M09), results/sir_inputs.mat (M05),
%          data/processed
%  Output: output/figures/figS14-S18_*.png
%          output/tables/allocation_correlations.csv, daniels_shen_similarity.csv

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg    = project_config();
obs    = load_observed_data(cfg);
sir_in = load(cfg.file_sir_inputs);
S = cfg.S; T = cfg.T;

opt = struct();
for vi = 1:2
    opt.(cfg.ve_names{vi}) = load(cfg.file_opt(cfg.ve_names{vi}), 'doses_opt');
end
pop_state    = sir_state_population(sir_in);
doses_brazil = brazil_dose_allocation(cfg);

past_inc   = sum(obs.cases_st(T-9:T, :), 1)' ./ pop_state * 1e5;   % 2016-2025 per 100K
target_idx = {cfg.target_ages.children, cfg.target_ages.elderly};
target_pop = [sum(obs.pop_2025(:, target_idx{1}), 2), sum(obs.pop_2025(:, target_idx{2}), 2)];
target_tag = {'children', 'elderly'};

dark_green = [0.0 0.8 0.0];
target_col = [0.8 0.1 0.1; 0.0 0.2 0.8];   % children red, elderly blue
lighten    = @(c, frac) c + (1 - c) * frac;
fig_no     = struct('scatter', [14 15], 'map', [16 17]);   % [Shen Daniels]

shp = load_brazil_shapes(cfg);
label_offset = containers.Map({'CE', 'DF', 'PE', 'RJ', 'RN', 'SE', 'GO'}, ...
                              {[-1 0], [0 0.8], [0 -0.8], [0 -0.8], [0 1], [0 -1], [-1 0]});

corr_rows = {};
for vi = 1:2
    ve = cfg.ve_names{vi};
    doses_opt = opt.(ve).doses_opt;   % [S x R0 x target]

    for tc = 1:2
        col = target_col(tc, :);
        per_target = doses_opt(:, :, tc) ./ target_pop(:, tc) * 1e5;   % [S x 3]
        dmin = min(per_target(:)); dmax = max(per_target(:));          % shared colour range

        for rc = 1:3
            tag = sprintf('%s_%s_%s', ve, target_tag{tc}, cfg.R0_names{rc});

            %% Figs. S14/S15: doses per 100K population vs past incidence
            bra = doses_brazil ./ pop_state * 1e5;
            dop = doses_opt(:, rc, tc) ./ pop_state * 1e5;
            xf  = linspace(min(past_inc) * 0.95, max(past_inc) * 1.05, 100);
            figure('Units', 'centimeters', 'Position', [0 0 7 5.5], 'PaperPositionMode', 'auto');
            plot(past_inc, bra, 'o', 'Color', dark_green, 'LineWidth', 2); hold on
            plot(past_inc, dop, 'o', 'Color', col, 'LineWidth', 2)
            plot(xf, polyval(polyfit(past_inc, bra, 1), xf), '--', 'Color', dark_green, 'LineWidth', 2)
            plot(xf, polyval(polyfit(past_inc, dop, 1), xf), '--', 'Color', col, 'LineWidth', 2)
            ylim([0 max(dop) * 1.1])
            if rc == 1
                legend('Brazil 2024/25', sprintf('Optimal (%s)', cfg.target_names{tc}), 'Location', 'north')
                legend('boxoff')
            end
            set(gca, 'FontName', 'SansSerif', 'FontSize', 10)
            save_figure(cfg, sprintf('figS%d_scatter_%s', fig_no.scatter(vi), tag))

            %% Figs. S16/S17: doses per 100K target population, choropleth
            nC = 256;
            pale = lighten(col, 0.85);
            cmap = [linspace(pale(1), col(1), nC)', linspace(pale(2), col(2), nC)', linspace(pale(3), col(3), nC)'];
            v_norm = (per_target(:, rc) - dmin) / max(dmax - dmin, eps);

            figure('Units', 'centimeters', 'Position', [0 0 11 12], 'PaperPositionMode', 'auto');
            hold on
            for i = 1:numel(shp)
                ci = max(1, min(nC, round(v_norm(shp(i).si) * (nC - 1)) + 1));
                for k = 1:numel(shp(i).parts)
                    xy = shp(i).parts{k};
                    patch(xy(:,1), xy(:,2), cmap(ci, :), 'EdgeColor', 'none', 'FaceAlpha', 1)
                end
            end
            for i = 1:numel(shp)
                off = [0 0];
                if isKey(label_offset, shp(i).code); off = label_offset(shp(i).code); end
                text(shp(i).cx + off(1), shp(i).cy + off(2), shp(i).code, ...
                     'HorizontalAlignment', 'center', 'FontName', 'SansSerif')
            end
            axis equal off
            colormap(gca, cmap); caxis([dmin dmax])
            cb = colorbar;
            cb.Label.String = 'Optimal vaccination per 100K (target pop)';
            cb.FontName = 'SansSerif'; cb.FontSize = 12;
            save_figure(cfg, sprintf('figS%d_map_%s', fig_no.map(vi), tag))

            %% correlations with past incidence
            [r_t, p_t] = corr(past_inc, per_target(:, rc));
            [r_p, p_p] = corr(past_inc, dop);
            corr_rows(end+1, :) = {cfg.ve_labels{vi}, cfg.target_names{tc}, cfg.R0_names{rc}, ...
                                   r_t, p_t, r_p, p_p}; %#ok<SAGROW>
        end
    end
end
[r_b, p_b] = corr(past_inc, doses_brazil ./ pop_state * 1e5);
corr_rows(end+1, :) = {'-', 'Brazil 2024/25', '-', nan, nan, r_b, p_b};
corr_tbl = cell2table(corr_rows, 'VariableNames', {'VE', 'target', 'R0', ...
    'r_per100k_target', 'p_per100k_target', 'r_per100k_population', 'p_per100k_population'});
writetable(corr_tbl, fullfile(cfg.table_dir, 'allocation_correlations.csv'))
disp(corr_tbl)

%% Fig. S18: Daniels vs Shen allocations per 100K target population
col_dan  = [0, 62, 116] / 255;
col_shen = [225, 36, 42] / 255;
darken   = @(c, frac) c .* (1 - frac);
col_q    = {[lighten(col_dan, 0.6); col_dan; darken(col_dan, 0.4)], ...
            [lighten(col_shen, 0.6); col_shen; darken(col_shen, 0.4)]};   % Q25 light .. Q75 dark
q_offset = [-0.28 0 0.28];

for tc = 1:2
    dan  = opt.daniels.doses_opt(cfg.reg_order, :, tc) ./ target_pop(cfg.reg_order, tc) * 1e5;
    shen = opt.shen.doses_opt(cfg.reg_order, :, tc)    ./ target_pop(cfg.reg_order, tc) * 1e5;
    ymax = max([dan(:); shen(:)]) * 1.1;

    figure('Units', 'centimeters', 'Position', [0 0 16 6], 'PaperPositionMode', 'auto');
    draw_region_background(cfg, ymax, 0.15)
    hd = gobjects(3, 1); hs = gobjects(3, 1);
    for k = 1:S
        for q = 1:3
            x0 = k + q_offset(q);
            plot([x0 x0], [dan(k, q) shen(k, q)], 'k:', 'LineWidth', 0.8)
            hd(q) = plot(x0, dan(k, q), 'o', 'MarkerFaceColor', col_q{1}(q, :), ...
                         'MarkerEdgeColor', col_q{1}(q, :), 'MarkerSize', 3.5);
            hs(q) = plot(x0, shen(k, q), 'o', 'MarkerFaceColor', col_q{2}(q, :), ...
                         'MarkerEdgeColor', col_q{2}(q, :), 'MarkerSize', 3.5);
        end
    end
    xlim([0.25 S+0.75]); ylim([0 ymax])
    xticks(1:S); xticklabels(cfg.state_codes(cfg.reg_order))
    set(gca, 'XTickLabelRotation', 90, 'FontName', 'SansSerif', 'FontSize', 13)
    if tc == 2
        legend([hd; hs], {'Daniels Q25', 'Daniels Q50', 'Daniels Q75', 'Shen Q25', 'Shen Q50', 'Shen Q75'}, ...
               'Location', 'northeast', 'FontSize', 13, 'NumColumns', 2); legend('boxoff')
    end
    save_figure(cfg, ['figS18_daniels_vs_shen_' target_tag{tc}])
end

%% Similarity of Daniels- and Shen-derived allocations
sim_rows = {};
for tc = 1:2
    for rc = 1:3
        dd = opt.daniels.doses_opt(:, rc, tc);
        ss = opt.shen.doses_opt(:, rc, tc);
        tvd = 0.5 * sum(abs(dd / sum(dd) - ss / sum(ss)));   % total variation distance
        sim_rows(end+1, :) = {cfg.target_names{tc}, cfg.R0_names{rc}, corr(dd, ss, 'Type', 'Pearson'), ...
                              corr(dd, ss, 'Type', 'Spearman'), 100 * (1 - tvd)}; %#ok<SAGROW>
    end
end
sim_tbl = cell2table(sim_rows, 'VariableNames', {'target', 'R0', 'pearson_r', 'spearman_rho', 'similarity_pct'});
writetable(sim_tbl, fullfile(cfg.table_dir, 'daniels_shen_similarity.csv'))
disp(sim_tbl)
