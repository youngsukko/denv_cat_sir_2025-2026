%% Appendix: Daniels vs Shen 비교 — 주별, q25/50/75 마커 (막대 없음)
close all; clc;

Sd = load('results/SIR_opt_result_dan.mat');   % Daniels
Ss = load('results/SIR_opt_result.mat');        % Shen
% 둘 다: doses_opt [S x 3(R0) x 2(target)], target_names, R0_names, state_codes, popbystate

state_codes  = Ss.state_codes;
target_names = Ss.target_names;
S = 27; T = 26; A5 = 17;

col_dan  = [0, 62, 116]/255;
col_shen = [225, 36, 42]/255;

lighten = @(c, frac) c + (1-c)*frac;
darken  = @(c, frac) c .* (1-frac);

% q25=옅게, q50=원색, q75=진하게
col_dan_q  = [lighten(col_dan,0.6);  col_dan;  darken(col_dan,0.4)];
col_shen_q = [lighten(col_shen,0.6); col_shen; darken(col_shen,0.4)];

ibge_of_state = [12,27,13,16,29,23,53,32,52,21, ...
                 31,50,51,15,25,26,22,41,33,24, ...
                 11,14,43,42,28,35,17]; %#ok<NASGU>

north_idx     = find(ismember(state_codes,{'AC','AM','AP','PA','RO','RR','TO'}));
northeast_idx = find(ismember(state_codes,{'AL','BA','CE','MA','PB','PE','PI','RN','SE'}));
centwest_idx  = find(ismember(state_codes,{'DF','GO','MS','MT'}));
southeast_idx = find(ismember(state_codes,{'ES','MG','RJ','SP'}));
south_idx     = find(ismember(state_codes,{'PR','RS','SC'})); %#ok<NASGU>

north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];
reg_order = [north_order,northeast_order,centwest_order,southeast_order,south_order];
codes_reord = state_codes(reg_order);

col_region = [1.000 0.498 0.055; 0.596 0.306 0.639; 0.090 0.745 0.812; ...
              0.173 0.627 0.173; 0.549 0.341 0.271];
reg_bound = [0, cumsum([numel(north_order),numel(northeast_order), ...
    numel(centwest_order),numel(southeast_order),numel(south_order)])];

% ── 연령별 인구 (타겟연령대 per-100K 계산용) ──────────────────
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
pop_mat = table2array(readtable('results/pop_age_state_year.csv'));
pop_mat = pop_mat(:, age_reorder_idx);
target_age_idx = {[2,3], 13:17};   % tc=1: children(5-14), tc=2: elderly(60+)
target_pop = zeros(S,2);
for si = 1:S
    pop_2025_si = pop_mat((si-1)*T+T, :);
    target_pop(si,1) = sum(pop_2025_si(target_age_idx{1}));
    target_pop(si,2) = sum(pop_2025_si(target_age_idx{2}));
end

mkdir('figures')

q_offset = [-0.28, 0, 0.28];   % q25/q50/q75 서브 위치

for tc = 1:2
    dan_per100k  = squeeze(Sd.doses_opt(reg_order,:,tc)) ./ target_pop(reg_order,tc) * 1e5;  % [S x 3]
    shen_per100k = squeeze(Ss.doses_opt(reg_order,:,tc)) ./ target_pop(reg_order,tc) * 1e5;  % [S x 3]

    figure('Units','centimeters','Position',[0 0 16 6],'PaperPositionMode','auto');
    hold on
    ymax_bar = max([dan_per100k(:); shen_per100k(:)])*1.1;
    for r = 1:5
        fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
             [0 0 1 1]*ymax_bar, ...
             col_region(r,:),'EdgeColor','none','FaceAlpha',0.15)
    end

    hdan  = gobjects(3,1);
    hshen = gobjects(3,1);
    for si = 1:S
        for q = 1:3
            x0 = si + q_offset(q);
            y_dan  = dan_per100k(si,q);
            y_shen = shen_per100k(si,q);
            plot([x0, x0], [y_dan, y_shen], 'k:', 'LineWidth', 0.8)
            hdan(q)  = plot(x0, y_dan,  'o', 'MarkerFaceColor', col_dan_q(q,:), ...
                'MarkerEdgeColor', col_dan_q(q,:),  'MarkerSize', 3.5);
            hshen(q) = plot(x0, y_shen, 'o', 'MarkerFaceColor', col_shen_q(q,:), ...
                'MarkerEdgeColor', col_shen_q(q,:), 'MarkerSize', 3.5);
        end
    end

    xlim([0.25 S+0.75]); ylim([0 ymax_bar])
    xticks(1:S); xticklabels(codes_reord)
    set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',13)
    % ylabel('Doses per 100K (target age group)')
    if tc == 2
    legend([hdan; hshen], {'Daniels Q25','Daniels Q50','Daniels Q75', ...
        'Shen Q25','Shen Q50','Shen Q75'}, ...
        'Location','northeast','FontSize',13,'NumColumns',2); legend('boxoff')
    end
    exportgraphics(gcf, sprintf('figures/fig_appendix_dan_shen_%s.jpg', ...
        lower(strrep(target_names{tc},' ','_'))), 'Resolution',600)
end
fprintf('Done.\n')

%% Daniels vs Shen 배분전략 유사도 (q25/q50/q75 × children/elderly)
qvals = [25 50 75];
fprintf('\n%-10s %-6s %10s %10s %10s\n', 'target','q','Pearson_r','Spearman','1-TVD(%%)')
sim_rows = {};
for tc = 1:2
    for q = 1:3
        d = Sd.doses_opt(:,q,tc);
        s = Ss.doses_opt(:,q,tc);

        r_pearson  = corr(d, s, 'Type','Pearson');
        r_spearman = corr(d, s, 'Type','Spearman');
        tvd = 0.5 * sum(abs(d/sum(d) - s/sum(s)));   % Total Variation Distance (분포로 봤을 때)
        sim_pct = (1 - tvd) * 100;                    % 얼마나 안 옮겨도 되는지 (%)

        fprintf('%-10s q%-5d %10.3f %10.3f %10.2f\n', ...
            target_names{tc}, qvals(q), r_pearson, r_spearman, sim_pct)
        sim_rows(end+1,:) = {target_names{tc}, sprintf('q%d',qvals(q)), ...
            r_pearson, r_spearman, sim_pct}; %#ok<SAGROW>
    end
end
T_sim = cell2table(sim_rows, 'VariableNames', ...
    {'target','R0_quantile','pearson_r','spearman_rho','similarity_pct_1minusTVD'});
mkdir('result_csv')
writetable(T_sim, 'result_csv/dan_shen_similarity.csv')
disp(T_sim)