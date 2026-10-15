%% Optimal allocation visualization
% Requires: results/SIR_opt_result.mat (from C_sir_opt.m)
%           results/cases_age_state_year.csv
%           results/pop_age_state_year.csv (age-reordered)
%
% 색상 규칙: 고령층(elderly) 타겟 관련 플롯은 파란색(dark_blue)으로 표시.
%           그 외 타겟은 기존 색(dark_red / magenta / black) 유지.
close all; clc;

load('results/SIR_opt_result_dan.mat')
% doses_opt [S × 3(R0) × 2(target)]
% daly_opt, daly_brazil [3 × 2]
% doses_annual_brazil [S × 1]
% popbystate [S × 1]
% target_names, R0_names, state_codes

S = 27; T = 26; A5 = 17;
annual_national = 5e6;

dark_green  = [0.0 0.8 0.0];
dark_red    = [0.8 0.1 0.1];
dark_blue   = [0.0 0.2 0.8];

lighten = @(c, frac) c + (1-c)*frac;

% 고령층 타겟 판별 키워드 (target_names 표기에 맞게 필요시 조정)
elderly_keywords = {'elderly','older','60+','고령'};

ibge_of_state = [12,27,13,16,29,23,53,32,52,21, ...
                 31,50,51,15,25,26,22,41,33,24, ...
                 11,14,43,42,28,35,17];

col_region = [1.000 0.498 0.055; 0.596 0.306 0.639; 0.090 0.745 0.812; ...
              0.173 0.627 0.173; 0.549 0.341 0.271];

north_idx     = find(ismember(state_codes,{'AC','AM','AP','PA','RO','RR','TO'}));
northeast_idx = find(ismember(state_codes,{'AL','BA','CE','MA','PB','PE','PI','RN','SE'}));
centwest_idx  = find(ismember(state_codes,{'DF','GO','MS','MT'}));
southeast_idx = find(ismember(state_codes,{'ES','MG','RJ','SP'}));
south_idx     = find(ismember(state_codes,{'PR','RS','SC'}));

north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];
reg_order = [north_order,northeast_order,centwest_order,southeast_order,south_order];
codes_reord = state_codes(reg_order);

region_of = zeros(S,1);
region_of(north_idx)=1; region_of(northeast_idx)=2; region_of(centwest_idx)=3;
region_of(southeast_idx)=4; region_of(south_idx)=5;

reg_bound = [0, cumsum([numel(north_order),numel(northeast_order), ...
    numel(centwest_order),numel(southeast_order),numel(south_order)])];

% ── 누적 인시던스 계산 (2014-2023, T-12:T-3) ──────────────────
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
cases_mat = table2array(readtable('results/cases_age_state_year.csv'));
cases_mat = cases_mat(:, age_reorder_idx);

past_inc = zeros(S,1);
for si = 1:S
    idx = (si-1)*T + (T-9:T);   % 2016-2025 
    past_inc(si) = sum(cases_mat(idx,:),'all') / popbystate(si) * 1e5;
end

% ── 연령별 인구 (타겟연령대 per-100K 계산용) ──────────────────
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
Shfile = shaperead('BR_UF_2023.shp');
cxvec = zeros(1,length(Shfile));
cyvec = zeros(1,length(Shfile));
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    valid = ~isnan(x) & ~isnan(y);
    pg = polyshape(x(valid),y(valid));
    [cx,cy] = centroid(pg);
    cxvec(i) = cx; cyvec(i) = cy;
end
cxvec(26) = cxvec(26) - 1;

label_offset = containers.Map( ...
    {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
     'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
     'RO','RR','RS','SC','SE','SP','TO'}, ...
    {[0,0],[0,0],[0,0],[0,0],[0,0],[-1,0],[0,0.8],[0,0],[0,0],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,-0.8],[0,0],[0,0],[0,-0.8],[0,1], ...
     [0,0],[0,0],[0,0],[0,0],[0,-1],[0,0],[0,0]});

for tc = 1:2
    % ── 타겟별 색상 결정: 고령층 타겟이면 파란색, 아니면 기존 색 ──
    is_elderly = any(contains(lower(target_names{tc}), elderly_keywords));
    if is_elderly
        opt_color   = dark_blue;
        patch_color = dark_blue;
    else
        opt_color   = dark_red;
        patch_color = dark_red;
    end

    %% ── Fig 3: 주별 bar — q25/50/75 한 그림에, Brazil alloc은 초록 마커 ──
    bra_r   = doses_annual_brazil(reg_order) ./ target_pop(reg_order,tc) * 1e5;
    opt_r25 = doses_opt(reg_order,1,tc)      ./ target_pop(reg_order,tc) * 1e5;
    opt_r50 = doses_opt(reg_order,2,tc)      ./ target_pop(reg_order,tc) * 1e5;
    opt_r75 = doses_opt(reg_order,3,tc)      ./ target_pop(reg_order,tc) * 1e5;

    col_q25 = lighten(opt_color, 0.65);
    col_q50 = lighten(opt_color, 0.30);
    col_q75 = opt_color;

    % figure('Units','centimeters','Position',[0 0 18 6],'PaperPositionMode','auto');
    % hold on
    % ymax_bar = max([bra_r; opt_r25; opt_r50; opt_r75])*1.2;
    % for r = 1:5
    %     fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
    %          [0 0 1 1]*ymax_bar, ...
    %          col_region(r,:),'EdgeColor','none','FaceAlpha',0.15)
    % end
    % bb = bar(1:S, [opt_r25, opt_r50, opt_r75], 0.85, 'grouped', 'EdgeColor','none');
    % bb(1).FaceColor = col_q25;
    % bb(2).FaceColor = col_q50;
    % bb(3).FaceColor = col_q75;
    % hbra = plot(1:S, bra_r, 'o', 'MarkerEdgeColor', dark_green, ...
    %     'MarkerFaceColor', dark_green, 'MarkerSize', 4, 'LineStyle','none');
    % 
    % xlim([0.25 S+0.75]); ylim([0 ymax_bar])
    % xticks(1:S); xticklabels(codes_reord)
    % set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',8)
    % ylabel('Doses per 100K (target age group)')
    % title(sprintf('%s', target_names{tc}), 'FontName','SansSerif','FontWeight','normal')
    % legend([bb, hbra], {'Optimal q25','Optimal q50','Optimal q75','Brazil alloc'}, ...
    %     'Location','northeast','FontSize',8); legend('boxoff')
    % exportgraphics(gcf, sprintf('figures/fig_bar_opt_%s.jpg', ...
    %     lower(strrep(target_names{tc},' ','_'))), 'Resolution',600)

    % ── 컬러바 범위를 R0(q25/50/75) 3개 지도 간 통일 ──────────────
    doses_per100k_target_all = doses_opt(:,:,tc) ./ target_pop(:,tc) * 1e5;  % [S x 3]
    dmin_tc = min(doses_per100k_target_all(:));
    dmax_tc = max(doses_per100k_target_all(:));

    for rc = 1:3
        doses_opt_si = doses_opt(:, rc, tc);

        %% ── Fig 1: Scatter (인시던스 vs 접종량) ─────────────────
        doses_bra_per100k = doses_annual_brazil ./ popbystate * 1e5;
        doses_opt_per100k = doses_opt_si        ./ popbystate * 1e5;

        ppp_bra = polyfit(past_inc, doses_bra_per100k, 1);
        ppp_opt = polyfit(past_inc, doses_opt_per100k, 1);
        xfit = linspace(min(past_inc)*0.95, max(past_inc)*1.05, 100);

        figure('Units','centimeters','Position',[0 0 7 5.5],'PaperPositionMode','auto');
        plot(past_inc, doses_bra_per100k, 'o', 'color', dark_green, 'LineWidth', 2)
        hold on
        plot(past_inc, doses_opt_per100k, 'o', 'color', opt_color, 'LineWidth', 2)
        plot(xfit, polyval(ppp_bra, xfit), '--', 'color', dark_green, 'LineWidth', 2)
        plot(xfit, polyval(ppp_opt, xfit), '--', 'color', opt_color, 'LineWidth', 2)
        % xlabel('Cumulative incidence (2016-2025)')
        % ylabel('Vaccination per 100K population')
        ylim([0 max([doses_opt_per100k;doses_opt_per100k]*1.1)])
        if rc == 1
        legend('Brazil 2024/25', sprintf('Optimal (%s)', target_names{tc}), ...
            'Location','north'); legend('boxoff')
        end
        set(gca,'FontName','SansSerif','FontSize',10)
        exportgraphics(gcf, sprintf('figures/fig_scatter_opt_%s_%s.jpg', ...
            lower(strrep(target_names{tc},' ','_')), R0_names{rc}), 'Resolution',600)

        %% ── Fig 2: Choropleth (타겟연령대 인구 10만명당, 자기색 옅음→진함, 컬러바) ──
        doses_opt_per100k_target = doses_opt_si ./ target_pop(:,tc) * 1e5;

        dmin = dmin_tc; dmax = dmax_tc;
        opt_norm = (doses_opt_per100k_target - dmin) / max(dmax - dmin, eps);   % 0~1

        nC = 256;
        pale_color = lighten(patch_color, 0.85);   % 타겟 색 자체의 옅은 버전 (흰색 아님, 다른 색 안섞임)
        cmap_local = [linspace(pale_color(1), patch_color(1), nC)', ...
                      linspace(pale_color(2), patch_color(2), nC)', ...
                      linspace(pale_color(3), patch_color(3), nC)'];

        figure('Units','centimeters','Position',[0 0 11 12],'PaperPositionMode','auto');
        hold on
        for i = 1:length(Shfile)
            x = Shfile(i).X; y = Shfile(i).Y;
            cd_uf = str2double(Shfile(i).CD_UF);
            si = find(ibge_of_state == cd_uf);
            nan_idx = find(isnan(x));
            starts = [1, nan_idx+1]; ends = [nan_idx-1, length(x)];
            for j = 1:length(starts)
                seg = starts(j):ends(j);
                if length(seg) >= 3
                    if ~isempty(si)
                        ci = max(1, min(nC, round(opt_norm(si)*(nC-1))+1));
                        patch(x(seg),y(seg),cmap_local(ci,:),'EdgeColor','none','FaceAlpha',1)
                    else
                        patch(x(seg),y(seg),[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.5)
                    end
                end
            end
        end
        for i = 1:length(Shfile)
            cd_uf = str2double(Shfile(i).CD_UF);
            si = find(ibge_of_state == cd_uf);
            if ~isempty(si)
                off = label_offset(state_codes{si});
                text(cxvec(i)+off(1), cyvec(i)+off(2), state_codes{si}, ...
                    'HorizontalAlignment','center','FontName','SansSerif')
            end
        end
        axis equal off
        colormap(gca, cmap_local)
        caxis([dmin dmax])
        cb = colorbar;
        cb.Label.String = 'Optimal vaccination per 100K (target pop)';
        cb.FontName = 'sansserif';
        cb.FontSize = 12;
        % title(sprintf('Optimal allocation: %s | R0=%s', target_names{tc}, R0_names{rc}), ...
        %     'FontName','SansSerif')
        exportgraphics(gcf, sprintf('figures/fig_choropleth_opt_%s_%s.jpg', ...
            lower(strrep(target_names{tc},' ','_')), R0_names{rc}), 'Resolution',600)

        fprintf('[%s / R0=%s] DALY: Brazil=%.1f  Optimal=%.1f  (%.2f%% reduction)\n', ...
            target_names{tc}, R0_names{rc}, daly_brazil(rc,tc), daly_opt(rc,tc), ...
            (daly_brazil(rc,tc)-daly_opt(rc,tc))/daly_brazil(rc,tc)*100)
    end
end
fprintf('Done.\n')

%% 누적 인시던스 vs 접종량(per100K, 타겟인구) 상관관계 — q25/50/75 × children/elderly
fprintf('\n%-12s %-8s %10s %10s\n', 'target', 'R0', 'corr(r)', 'p-value')
for tc = 1:2
    for rc = 1:3
        doses_per100k = doses_opt(:,rc,tc) ./ target_pop(:,tc) * 1e5;
        [r, p] = corr(past_inc, doses_per100k);
        fprintf('%-12s %-8s %10.3f %10.4f\n', target_names{tc}, R0_names{rc}, r, p)
    end
end