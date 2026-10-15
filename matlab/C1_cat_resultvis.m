%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Catalytic model viz (Stan PPC 기반, MATLAB 재계산 없음)
%  Input:
%    results/catalytic_ppc_draws.mat          (op_seroposall_*, model_inc, cum_cases_age, ...)
%    results/catalytic_lambda_kappa_draws.mat (lambda/kappa/omega/delta/rho/lambda_pre draws)
%    results/cases_age_state_year.csv
%    results/pop_age_state_year.csv
%    results/incidence_per100.csv
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;
%%
fsz = 11;
fsz_pan = 8;
%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;
years       = 2000:2025;
agerangeidx = 0:5:80;

agelabels = arrayfun(@(i) sprintf('%d-',agerangeidx(i)), 1:16, 'UniformOutput',false);
agelabels{17} = '80+';

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};
state_names = {'Acre','Alagoas','Amazonas','Amapá','Bahia','Ceará', ...
               'Distrito Federal','Espírito Santo','Goiás','Maranhão', ...
               'Minas Gerais','Mato Grosso do Sul','Mato Grosso','Pará', ...
               'Paraíba','Pernambuco','Piauí','Paraná','Rio de Janeiro', ...
               'Rio Grande do Norte','Rondônia','Roraima','Rio Grande do Sul', ...
               'Santa Catarina','Sergipe','São Paulo','Tocantins'};

old_to_new = [21,1,3,22,14,4,27,10,17,6,20,15,16,2,25,5,11,19,26,18,24,23,12,13,9,7];

north_codes     = {'AC','AM','AP','PA','RO','RR','TO'};
northeast_codes = {'AL','BA','CE','MA','PB','PE','PI','RN','SE'};
centwest_codes  = {'DF','GO','MS','MT'};
southeast_codes = {'ES','MG','RJ','SP'};
south_codes     = {'PR','RS','SC'};

north_idx     = find(ismember(state_codes, north_codes));
northeast_idx = find(ismember(state_codes, northeast_codes));
centwest_idx  = find(ismember(state_codes, centwest_codes));
southeast_idx = find(ismember(state_codes, southeast_codes));
south_idx     = find(ismember(state_codes, south_codes));

% 지역 매핑 (1=North, 2=Northeast, 3=Central-west, 4=Southeast, 5=South)
region_of = containers.Map(state_codes, ...
    [1,2,1,1,2,2,3,4,3,2, ...
     4,3,3,1,2,2,2,5,4,2, ...
     1,1,5,5,2,4,1]);

col_region = [1.000 0.498 0.055; 0.596 0.306 0.639; 0.090 0.745 0.812; ...
              0.173 0.627 0.173; 0.549 0.341 0.271];

cmap = [150 150 150; 255 221 102; 255 153 51; 255 102 51; 204 0 0] / 255;

col_foi    = [0.839 0.153 0.157];
col_report = [0.580 0.404 0.741];
col_model  = [0.15 0.15 0.15];
col_data   = [0.85 0.10 0.10];
lighten    = @(c, frac) c + (1-c)*frac;
col_foi_light    = lighten(col_foi, 0.65);
col_report_light = lighten(col_report, 0.65);
col_model_light  = lighten(col_model, 0.65);

NCOLS = 6;
state_to_tile = zeros(1,S);
tile_def = {'RR',1,2; 'AP',1,3; 'AM',2,2; 'PA',2,3; 'MA',2,4; 'CE',2,5;
            'AC',3,1; 'RO',3,2; 'TO',3,3; 'PI',3,4; 'PB',3,5; 'RN',3,6;
            'MT',4,2; 'GO',4,3; 'BA',4,4; 'PE',4,5; 'AL',4,6;
            'MS',5,2; 'DF',5,3; 'MG',5,4; 'SE',5,5;
            'SP',6,3; 'RJ',6,4; 'ES',6,5;
            'PR',7,3; 'SC',7,4; 'RS',8,3};
for k = 1:size(tile_def,1)
    si = find(strcmp(state_codes, tile_def{k,1}));
    state_to_tile(si) = (tile_def{k,2}-1)*NCOLS + tile_def{k,3};
end
legend_tile = 5;

pct025 = @(x) quantile(x, 0.025);
pct975 = @(x) quantile(x, 0.975);
pct25  = @(x) quantile(x, 0.25);
pct75  = @(x) quantile(x, 0.75);
pmed   = @(x) median(x);

mkdir('figures')

%% ── helper: incidence CSV 로더 (구포맷/신포맷 자동 인식) ────────
function inc_mat = load_incidence(path, T, S, cases_mat, pop_mat)
    tbl = readtable(path);
    vn  = tbl.Properties.VariableNames;
    is_year = strcmpi(vn, 'year');
    raw = table2array(tbl(:, ~is_year));
    if isequal(size(raw), [T, S])
        inc_mat = raw;
    elseif size(raw,1) == S*T
        warning('incidence_per100.csv가 연령별 포맷입니다. cases/pop로부터 재계산합니다.');
        inc_mat = zeros(T, S);
        for s = 1:S
            idx = (s-1)*T + (1:T);
            c = sum(cases_mat(idx,:), 2);
            p = sum(pop_mat(idx,:), 2);
            inc_mat(:,s) = c ./ p * 100;
        end
    else
        error('incidence_per100.csv 크기를 인식할 수 없습니다: [%d x %d]', size(raw,1), size(raw,2));
    end
end

%% ── helper: validation 점 플롯 (나이값/구간경계 자동 인식) ──────
function draw_vali_pre2000(seropos_pre_0, si, agerangeidx, ...
                           tmp_dat_ages, tmp_prop, minidx, maxidx, title_str)
    p_band = squeeze(1 - seropos_pre_0(:,si,:));   % [n_draw x A5]
    p_s = sort(p_band);
    fill([agerangeidx, fliplr(agerangeidx)], [p_s(maxidx,:), fliplr(p_s(minidx,:))], ...
         [0.8 0.4 0.7],'EdgeColor','none','FaceAlpha',0.3); hold on
    plot(agerangeidx, median(p_s), '-','color',[0.6 0 0.5],'LineWidth',1.5)

    tmp_prop = tmp_prop(:)';
    if numel(tmp_dat_ages) == numel(tmp_prop)
        obs_x = tmp_dat_ages;
    else
        obs_x = (tmp_dat_ages(1:end-1)+tmp_dat_ages(2:end))/2;
        obs_x(end) = min(obs_x(end), 80);
    end
    plot(obs_x, tmp_prop, 'ko','MarkerFaceColor','k','MarkerSize',4)

    xlabel('Age'); ylabel('Seroprevalence')
    xlim([0 80]); ylim([0 1])
    title(title_str,'FontSize',8,'FontWeight','normal')
    set(gca,'FontName','SansSerif','FontSize',8)
end

function draw_vali_post2000(op_seroposall_0, si, ti, agerangeidx, ...
                             tmp_dat_ages, tmp_prop, minidx, maxidx, title_str)
    X0b = squeeze(op_seroposall_0(:,si,ti,:));
    sero = sort(1 - X0b);
    fill([agerangeidx, fliplr(agerangeidx)], [sero(maxidx,:), fliplr(sero(minidx,:))], ...
         [0.3 0.5 0.8],'EdgeColor','none','FaceAlpha',0.3); hold on
    plot(agerangeidx, median(sero), '-','color',[0 0.2 0.7],'LineWidth',1.5)

    tmp_prop = tmp_prop(:)';
    if numel(tmp_dat_ages) == numel(tmp_prop)
        obs_x = tmp_dat_ages;
    else
        obs_x = (tmp_dat_ages(1:end-1)+tmp_dat_ages(2:end))/2;
        obs_x(end) = min(obs_x(end), 80);
    end
    plot(obs_x, tmp_prop, 'ro','MarkerFaceColor','r','MarkerSize',4)

    xlabel('Age'); ylabel('Seroprevalence')
    xlim([0 80]); ylim([0 1])
    title(sprintf('%s (t=%d, %d)', title_str, ti, ti+1999),'FontSize',8,'FontWeight','normal')
    set(gca,'FontName','SansSerif','FontSize',8)
end

%% ── 1. 데이터 로드 (전부 .mat / .csv, 재계산 없음) ──────────────
fprintf('Loading catalytic_ppc_draws.mat...\n');
load('results/catalytic_ppc_draws.mat', ...
    'op_seroposall_0','op_seroposall_1','op_seroposall_2', ...
    'op_seroposall_3','op_seroposall_4', ...
    'op_seroposall_pre_0','op_seroposall_pre_1','op_seroposall_pre_2', ...
    'op_seroposall_pre_3','op_seroposall_pre_4', ...
    'model_inc','cum_cases_age','cum_case_age_prop_model');

fprintf('Loading catalytic_lambda_kappa_draws.mat...\n');
load('results/catalytic_lambda_kappa_draws.mat', ...
    'lambda_draws','kappa_draws','omega_draws','delta_draws', ...
    'rho_draws','lambda_pre_draws');
% lambda_draws,rho_draws: [n_draw x S x T]
% kappa_draws,omega_draws: [n_draw x A5]
% delta_draws: [n_draw x 4]
% lambda_pre_draws: [n_draw x S]

n_draw = size(op_seroposall_0,1);
minidx = max(1,   round(n_draw*0.025));
maxidx = min(n_draw, round(n_draw*0.975));
q25idx = max(1,   round(n_draw*0.25));
q75idx = min(n_draw, round(n_draw*0.75));

fprintf('Loading CSV files (observed data only)...\n');
cases_mat = table2array(readtable('results/cases_age_state_year.csv'));   % [S*T x A5]
pop_mat   = table2array(readtable('results/pop_age_state_year.csv'));
% ── CSV 컬럼이 R의 sort(unique(age_group))로 잘못 정렬되어 있음 ──
% (알파벳순 정렬이라 "5-9"가 "45-49"와 "50-54" 사이에 끼어듦)
% 올바른 연령순으로 재배열
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order = agelabels;   % {'0-','5-','10-',...,'80+'} — 이미 메인 스크립트에 있음

% agelabels는 '0-' 형태라 비교를 위해 같은 포맷으로 변환
wrong_order_short = cellfun(@(s) strsplit(s,'-'), wrong_order, 'UniformOutput', false);
wrong_order_short = cellfun(@(c) [c{1} '-'], wrong_order_short, 'UniformOutput', false);
wrong_order_short{strcmp(wrong_order,'80+')} = '80+';

[~, reorder_idx] = ismember(correct_order, wrong_order_short);
assert(all(reorder_idx > 0), 'CSV 컬럼 매핑 실패 — wrong_order 확인 필요');

cases_mat = cases_mat(:, reorder_idx);
pop_mat   = pop_mat(:, reorder_idx);



inc_mat   = load_incidence('results/incidence_per100.csv', T, S, cases_mat, pop_mat);  % [T x S]
inc_max   = max(inc_mat(:)) * 1.05;

% 주별 연도별 케이스/인구 [T x S] (연령 합산, 관측 데이터)
data_case_bystate = zeros(T, S);
data_pop_bystate  = zeros(T, S);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_case_bystate(:,s) = sum(cases_mat(idx,:), 2);
    data_pop_bystate(:,s)  = sum(pop_mat(idx,:),   2);
end
data_inc_bystate = data_case_bystate ./ data_pop_bystate * 100;

% 누적 케이스(전 기간 합) 연령분포 [S x A5] (관측 데이터)
data_case_age_bystate = zeros(S, A5);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_case_age_bystate(s,:) = sum(cases_mat(idx,:), 1);
end
data_case_age_prop = data_case_age_bystate ./ sum(data_case_age_bystate,2);

% 지역별 케이스/인구/연령분포
data_case_reg = [sum(data_case_bystate(:,north_idx),2), sum(data_case_bystate(:,northeast_idx),2), ...
                 sum(data_case_bystate(:,centwest_idx),2), sum(data_case_bystate(:,southeast_idx),2), ...
                 sum(data_case_bystate(:,south_idx),2)];
data_pop_reg  = [sum(data_pop_bystate(:,north_idx),2), sum(data_pop_bystate(:,northeast_idx),2), ...
                 sum(data_pop_bystate(:,centwest_idx),2), sum(data_pop_bystate(:,southeast_idx),2), ...
                 sum(data_pop_bystate(:,south_idx),2)];
data_case_age_reg = zeros(5, A5);
for s = 1:S
    cas = data_case_age_bystate(s,:);
    if     ismember(state_codes{s}, north_codes);     data_case_age_reg(1,:) = data_case_age_reg(1,:) + cas;
    elseif ismember(state_codes{s}, northeast_codes); data_case_age_reg(2,:) = data_case_age_reg(2,:) + cas;
    elseif ismember(state_codes{s}, centwest_codes);  data_case_age_reg(3,:) = data_case_age_reg(3,:) + cas;
    elseif ismember(state_codes{s}, southeast_codes); data_case_age_reg(4,:) = data_case_age_reg(4,:) + cas;
    elseif ismember(state_codes{s}, south_codes);     data_case_age_reg(5,:) = data_case_age_reg(5,:) + cas;
    end
end

fprintf('Loading catalytic_lambda_kappa_draws.mat...\n');
load('results/catalytic_lambda_kappa_draws.mat', ...
    'lambda_draws','kappa_draws','omega_draws','delta_draws', ...
    'rho_draws','lambda_pre_draws');

% ── kappa_age 평균 정규화 (해석 가능한 FOI 스케일) ──────────────
% kappa_adj: 평균 1로 정규화된 연령별 FOI 배율
% lambda_adj/lambda_pre_adj: kappa 평균만큼 보정되어 실제 FOI 스케일을 가짐
% mean_kappa = mean(kappa_draws, 2);             % [n_draw x 1]
mean_kappa = 1;             % [n_draw x 1]
kappa_adj  = kappa_draws ./ mean_kappa;        % [n_draw x A5], 평균=1

lambda_adj = lambda_draws .* mean_kappa;       % [n_draw x S x T] (implicit expansion)
% lambda_adj = lambda_draws*0;
% for i = 1:S
%     lambda_adj(:,i,:) = lambda_draws(:,i,:).*mean_kappa;
% end
lambda_pre_adj = lambda_pre_draws .* mean_kappa;  % [n_draw x S]

% 누적(2000-2025) 인구 by 연령대, by 주 [S x A5]  (분모용)
data_pop_age_bystate = zeros(S, A5);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_pop_age_bystate(s,:) = sum(pop_mat(idx,:), 1);
end

% 누적(2000-2025) 총인구, by 주 [S x 1]  (분모 고정값)
data_pop_total_bystate = sum(data_pop_age_bystate, 2);   % [S x 1]
%% ── 2. 지역별 인시던스 시계열 ─────────────────────────────────
% figure('Units','centimeters','Position',[0 0 12 6],'PaperPositionMode','auto');
% ppp1 = plot(years, data_case_reg ./ data_pop_reg * 100, 'LineWidth', 1.5);
% for r=1:5; ppp1(r).Color = col_region(r,:); end
% legend('North','Northeast','Central-west','Southeast','South', 'Location','northwest'); legend('boxoff');
% xlim([years(1)-0.5, years(end)+0.5]); xticks(2000:4:2025)
% ylabel('Incidence (per 100)'); set(gca,'FontName','SansSerif','FontSize',fsz)
% exportgraphics(gcf,'figures/fig_incidence_region.jpg','Resolution',600)

%% ── 3. 지역별 연령분포 ────────────────────────────────────────
% figure('Units','centimeters','Position',[0 0 8 6],'PaperPositionMode','auto');
% ppp1 = plot(agerangeidx, (data_case_age_reg ./ sum(data_case_age_reg,2))', 'LineWidth',1.5);
% for r=1:5; ppp1(r).Color = col_region(r,:); end
% ylabel('Proportion by age group'); xlabel('Age group')
% xticklabels(agelabels); xticks(agerangeidx); xlim([-0.5, 80.5])
% set(gca,'FontName','SansSerif','FontSize',10)
% exportgraphics(gcf,'figures/fig_age_dist_region.jpg','Resolution',600)

%% ── 4. 연령별 FOI 보정인자 (kappa_age) ──────────────────────
figure('Units','centimeters','Position',[0 0 7 6],'PaperPositionMode','auto');
for i = 1:A5
    vals = kappa_adj(:, i);
    plot([i i], [pct025(vals), pct975(vals)], 'color', col_foi, 'LineWidth', 1.5);
    hold on
    % plot([i i], [pct25(vals), pct75(vals)], 'color', col_foi, 'LineWidth', 1.5*2.2);
    plot(i, pmed(vals), 'ko', 'LineWidth', 1.5, 'MarkerFaceColor', col_foi, 'MarkerSize', 3);
end
ylabel('FOI multiplier by age'); xlabel('Age group');
xlim([0 18]); 
h95 = plot(nan,nan,'-','color',col_foi,'LineWidth',1.5);
% hiqr = plot(nan,nan,'-','color',col_foi,'LineWidth',1.5*2.2);
hmed = plot(nan,nan,'ko','MarkerFaceColor',col_foi,'MarkerSize',3);
% legend([hmed hiqr h95], {'Median','IQR (25-75%)','95% CrI'}, 'Location','south'); legend('boxoff')
% legend([hmed h95], {'Median','IQR (25-75%)','95% CrI'}, 'Location','south'); legend('boxoff')
legend([hmed h95], {'Median','95% CrI'}, 'Location','south'); legend('boxoff')
set(gca,'xtick',1:17,'XTickLabel',agelabels,'FontName','SansSerif','FontSize',11)
% exportgraphics(gcf,'figures/fig_kappa_age.jpg','Resolution',600)
%% ── 5. 연령별 보고율 보정인자 (omega_age) ────────────────────
figure('Units','centimeters','Position',[0 0 7 6],'PaperPositionMode','auto');
for i = 1:A5
    vals = omega_draws(:, i);
    plot([i i], [pct025(vals), pct975(vals)], 'color', col_report, 'LineWidth', 1.5);
    hold on
    % plot([i i], [pct25(vals), pct75(vals)], 'color', col_report, 'LineWidth', 1.5*2.2);
    plot(i, pmed(vals), 'ko', 'LineWidth', 1.5, 'MarkerFaceColor', col_report, 'MarkerSize', 3);
end
ylabel('Report rate multiplier'); xlabel('Age group');
xlim([0 18]); ylim([0 1])
h95 = plot(nan,nan,'-','color',col_report,'LineWidth',1.5);
% hiqr = plot(nan,nan,'-','color',col_report,'LineWidth',1.5*2.2);
hmed = plot(nan,nan,'ko','MarkerFaceColor',col_report,'MarkerSize',3);
% legend([hmed hiqr h95], {'Median','IQR (25-75%)','95% CrI'}, 'Location','south'); legend('boxoff')
legend([hmed h95], {'Median','95% CrI'}, 'Location','south'); legend('boxoff')
set(gca,'xtick',1:17,'XTickLabel',agelabels,'FontName','SansSerif','FontSize',fsz)
exportgraphics(gcf,'figures/fig_omega_age.jpg','Resolution',600)

%% ── 6. 감염 순서별 증상발현율 (delta) ───────────────────────
figure('Units','centimeters','Position',[0 0 7 8],'PaperPositionMode','auto');
for k = 1:4
    vals = delta_draws(:, k);
    boxchart(k*ones(size(vals)), vals, 'BoxFaceColor', col_report, 'MarkerStyle', 'none');
    hold on
end
xlim([0.5 4.5]); xticks(1:4)
xticklabels({'Primary','Secondary','Tertiary','Quaternary'})
ylabel('Infection to apparent illness ratio')
set(gca,'FontName','SansSerif','FontSize',fsz)
exportgraphics(gcf,'figures/fig_delta.jpg','Resolution',600)

%% ── 7. Pre-2000 FOI (lambda_pre) ────────────────────────────
state_order = {'RO','AC','AM','RR','PA','AP','TO', ...
               'MA','PI','CE','RN','PB','PE','AL','SE','BA', ...
               'MS','MT','GO','DF', ...
               'MG','RJ','SP','ES', ...
               'PR','SC','RS'};
region_labels7 = [ones(1,7), 2*ones(1,9), 3*ones(1,4), 4*ones(1,4), 5*ones(1,3)];
reorder_idx = cellfun(@(c) find(strcmp(state_codes, c)), state_order);

vals_reordered = 4*lambda_pre_adj(:, reorder_idx);   % [n_draw x 27]
group_s = repmat(1:S, n_draw, 1);
group_r = repmat(region_labels7, n_draw, 1);

figure('Units','centimeters','Position',[0 0 9 12],'PaperPositionMode','auto');
h = boxchart(group_s(:), vals_reordered(:), 'GroupByColor', group_r(:), ...
    'MarkerStyle','none','BoxWidth',3,'BoxFaceAlpha',1);
for r = 1:5
    set(h(r),'BoxFaceColor',col_region(r,:));
end
xticks(linspace(0.55,S+0.5,S)); xticklabels(state_order); xlim([0 S+1])
ylabel('Pre-2000 FOI')
view([90 90])
ylim([0 0.05])
set(gca,'FontName','SansSerif','FontSize',fsz)
exportgraphics(gcf,'figures/fig_lambda_pre.jpg','Resolution',600)

%% ── 8. 주별/연도별 FOI 패널맵 ───────────────────────────────
figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8, NCOLS, 'TileSpacing','tight','Padding','tight');
for si = 1:S
    tno = state_to_tile(si); if tno == 0; continue; end
    nexttile(tno)

    foi_med = zeros(1,T); foi_lo95 = zeros(1,T); foi_hi95 = zeros(1,T);
    foi_lo50 = zeros(1,T); foi_hi50 = zeros(1,T);
    for ti = 1:T
        vals = 4 * squeeze(lambda_adj(:, si, ti));
        foi_lo95(ti) = pct025(vals); foi_hi95(ti) = pct975(vals);
        foi_lo50(ti) = pct25(vals);  foi_hi50(ti) = pct75(vals);
        foi_med(ti)  = pmed(vals);
    end

    yyaxis right
    plot(years, inc_mat(:,si), 'k:', 'LineWidth', 1.2)
    ylim([0 inc_max]); ax = gca; ax.YAxis(2).Color = 'k'; ax.YAxis(1).Color = col_foi;

    yyaxis left
    % fill([years, fliplr(years)], [foi_hi95, fliplr(foi_lo95)], ...
    %      col_foi_light, 'EdgeColor','none','FaceAlpha',0.5); hold on
    fill([years, fliplr(years)], [foi_hi50, fliplr(foi_lo50)], ...
         col_foi, 'EdgeColor','none','FaceAlpha',0.5); hold on
    plot(years, foi_med, 'color', col_foi, 'LineWidth', 1.2)
    ylim([0 0.3]); xlim([years(1)-0.5 years(end)+0.5])
    xticks([2000 2010 2020]); xticklabels({'00','10','20'})
    title(state_names{si},'FontSize',8,'FontWeight','normal')
    set(gca,'FontName','SansSerif','FontSize',8)
end
nexttile(legend_tile); axis off; hold on
% h1 = fill(nan,nan,col_foi_light,'FaceAlpha',0.5,'EdgeColor','none');
h2 = fill(nan,nan,col_foi,'FaceAlpha',0.5,'EdgeColor','none');
h3 = plot(nan,nan,'-','color',col_foi,'LineWidth',1.5);
h4 = plot(nan,nan,'k:','LineWidth',1.5);
% legend([h3 h2 h1 h4], {'Median FOI','IQR (25-75%)','95% CrI','Incidence (data, per 100)'}, ...
%     'Location','north','FontName','SansSerif','FontSize',11); legend('boxoff')
legend([h3 h2 h4], {'Median FOI','IQR (25-75%)','Incidence (data, per 100)'}, ...
    'Location','north','FontName','SansSerif','FontSize',11); legend('boxoff')
exportgraphics(gcf,'figures/fig_foi_panelmap.jpg','Resolution',600)

%% ── 9. 주별/연도별 보고율 패널맵 ────────────────────────────
figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8, NCOLS, 'TileSpacing','tight','Padding','tight');
for si = 1:S
    tno = state_to_tile(si); if tno == 0; continue; end
    nexttile(tno)

    rho_med = zeros(1,T); rho_lo95 = zeros(1,T); rho_hi95 = zeros(1,T);
    rho_lo50 = zeros(1,T); rho_hi50 = zeros(1,T);
    for ti = 1:T
        vals = squeeze(rho_draws(:, si, ti));
        rho_lo95(ti) = pct025(vals); rho_hi95(ti) = pct975(vals);
        rho_lo50(ti) = pct25(vals);  rho_hi50(ti) = pct75(vals);
        rho_med(ti)  = pmed(vals);
    end

    yyaxis right
    plot(years, inc_mat(:,si), 'k:', 'LineWidth', 1.2)
    ylim([0 inc_max]); ax = gca; ax.YAxis(2).Color = 'k'; ax.YAxis(1).Color = col_report;

    yyaxis left
    fill([years, fliplr(years)], [rho_hi95, fliplr(rho_lo95)], ...
         col_report_light, 'EdgeColor','none','FaceAlpha',0.5); hold on
    fill([years, fliplr(years)], [rho_hi50, fliplr(rho_lo50)], ...
         col_report, 'EdgeColor','none','FaceAlpha',0.5);
    plot(years, rho_med, 'color', col_report, 'LineWidth', 1.2)
    ylim([0 1]); xlim([years(1)-0.5 years(end)+0.5])
    xticks([2000 2010 2020]); xticklabels({'00','10','20'})
    title(state_names{si},'FontSize',8,'FontWeight','normal')
    set(gca,'FontName','SansSerif','FontSize',8)
end
nexttile(legend_tile); axis off; hold on
h1 = fill(nan,nan,col_report_light,'FaceAlpha',0.5,'EdgeColor','none');
h2 = fill(nan,nan,col_report,'FaceAlpha',0.5,'EdgeColor','none');
h3 = plot(nan,nan,'-','color',col_report,'LineWidth',1.5);
h4 = plot(nan,nan,'k:','LineWidth',1.5);
legend([h3 h2 h1 h4], {'Median report rate','IQR (25-75%)','95% CrI','Incidence (data, per 100)'}, ...
    'Location','north','FontName','SansSerif','FontSize',11); legend('boxoff')
exportgraphics(gcf,'figures/fig_rho_panelmap.jpg','Resolution',600)

%% ── 10. helper: 면역 프로파일 패널 ───────────────────────────────
function plot_immune_panel(X0r, X1r, X2r, X3r, X4r, agerangeidx, cmap, sname)
    denom = X0r + X1r + X2r + X3r + X4r;
    X0n = X0r./denom; X1n = X1r./denom; X2n = X2r./denom;
    X3n = X3r./denom; X4n = X4r./denom;

    bb = bar(agerangeidx, [mean(X0n);mean(X1n);mean(X2n);mean(X3n);mean(X4n)], ...
        1.1,'stacked','EdgeColor','none');
    for k=1:5; bb(k).FaceColor=cmap(k,:); end
    hold on
    plot(agerangeidx, mean(X0n)+mean(X1n), 'k-', 'LineWidth',1.5)
    xlim([0 80]); ylim([0 1]); xticks([0 20 40 60 80])
    set(gca,'FontName','SansSerif','FontSize',8)
    title(sname,'FontSize',11,'FontWeight','normal')
end

function add_legend_tile(tile_no, cmap)
    nexttile(tile_no); axis off; hold on
    bb = bar(nan*ones(1,17), nan*ones(5,17), 0.95,'stacked','EdgeColor','none');
    for k=1:5; bb(k).FaceColor=cmap(k,:); end
    bb(1).BaseLine.Visible='off';
    lg = legend({'0','1','2','3','4'},'Location','northwest','FontSize',11,'NumColumns',3);
    legend('boxoff'); title(lg,'Prior exposure','FontWeight','normal')
end

%% ── 11. 2000년 기준 면역 프로파일 (Stan에서 직접 받은 pre-2000) ──
fprintf('Fig: 2000 baseline immune profile...\n');
figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','tight')
for si = 1:S
    tno = state_to_tile(si); if tno==0; continue; end
    nexttile(tno)
    X0b = squeeze(op_seroposall_pre_0(:,si,:));
    X1b = squeeze(op_seroposall_pre_1(:,si,:));
    X2b = squeeze(op_seroposall_pre_2(:,si,:));
    X3b = squeeze(op_seroposall_pre_3(:,si,:));
    X4b = squeeze(op_seroposall_pre_4(:,si,:));
    plot_immune_panel(X0b,X1b,X2b,X3b,X4b,agerangeidx,cmap,state_codes{si})
end
add_legend_tile(legend_tile, cmap)
% sgtitle('Immune profile — Year 2000 (pre-2000 baseline)','FontSize',11,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_immune_2000.jpg','Resolution',600)

%% ── 12. 2025년 기준 면역 프로파일 ───────────────────────────────
fprintf('Fig: 2025 immune profile...\n');
figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','tight')
for si = 1:S
    tno = state_to_tile(si); if tno==0; continue; end
    nexttile(tno)
    X0b=squeeze(op_seroposall_0(:,si,T,:));
    X1b=squeeze(op_seroposall_1(:,si,T,:));
    X2b=squeeze(op_seroposall_2(:,si,T,:));
    X3b=squeeze(op_seroposall_3(:,si,T,:));
    X4b=squeeze(op_seroposall_4(:,si,T,:));
    plot_immune_panel(X0b,X1b,X2b,X3b,X4b,agerangeidx,cmap,state_codes{si})
end
add_legend_tile(legend_tile, cmap)
% sgtitle('Immune profile — Year 2025','FontSize',11,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_immune_2025.jpg','Resolution',600)

%% ── 13. 연도별 스냅샷 ───────────────────────────────────────────
% for yearidx = [5 10 15 20]
%     fprintf('Fig: year snapshot %d...\n', yearidx+2000);
%     figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto','Color','white');
%     tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','compact')
%     for si = 1:S
%         tno = state_to_tile(si); if tno==0; continue; end
%         nexttile(tno)
%         X0b=squeeze(op_seroposall_0(:,si,yearidx,:));
%         X1b=squeeze(op_seroposall_1(:,si,yearidx,:));
%         X2b=squeeze(op_seroposall_2(:,si,yearidx,:));
%         X3b=squeeze(op_seroposall_3(:,si,yearidx,:));
%         X4b=squeeze(op_seroposall_4(:,si,yearidx,:));
%         denom = X0b+X1b+X2b+X3b+X4b;
%         bb = bar(agerangeidx, [mean(X0b./denom);mean(X1b./denom);mean(X2b./denom); ...
%             mean(X3b./denom);mean(X4b./denom)], 0.95,'stacked','EdgeColor','none');
%         for k=1:5; bb(k).FaceColor=cmap(k,:); end
%         hold on
%         plot(agerangeidx, mean(X0b./denom)+mean(X1b./denom), 'c-','LineWidth',2)
%         xlim([0 80]); ylim([0 1]); xticks([0 20 40 60 80])
%         set(gca,'FontName','SansSerif','FontSize',7)
%         title(state_names{si},'FontSize',8,'FontWeight','normal')
%     end
%     add_legend_tile(legend_tile, cmap)
%     sgtitle(sprintf('Immune profile — Year %d', yearidx+2000),'FontSize',11,'FontName','SansSerif')
%     exportgraphics(gcf,sprintf('figures/fig_immune_%d.jpg',yearidx+2000),'Resolution',600)
% end

%% ── 14. Colorful seroprevalence trend ──────────────────────────
fprintf('Fig: colorful seroprevalence trend...\n');
figure('Units','centimeters','Position',[0 0 20 23],'PaperPositionMode','auto','Color','white');
tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','compact')
for si = 1:S
    tno = state_to_tile(si); if tno==0; continue; end
    nexttile(tno); hold on
    for yearidx = 1:T
        X0b = squeeze(op_seroposall_0(:,si,yearidx,:));
        X1b = squeeze(op_seroposall_1(:,si,yearidx,:));
        X2b = squeeze(op_seroposall_2(:,si,yearidx,:));
        X3b = squeeze(op_seroposall_3(:,si,yearidx,:));
        X4b = squeeze(op_seroposall_4(:,si,yearidx,:));
        denom = X0b+X1b+X2b+X3b+X4b;
        sero1 = mean(X0b./denom) + mean(X1b./denom);
        plot(agerangeidx, sero1, 'color',[(yearidx-1)/(T-1), 1-(yearidx-1)/(T-1), 1],'LineWidth',1.5)
    end
    xlim([0 80]); ylim([0 1]); xticks([0 20 40 60 80])
    set(gca,'FontName','SansSerif','FontSize',8,'XTickLabelRotation',0)
    title(state_names{si},'FontSize',8,'FontWeight','normal')
end
exportgraphics(gcf,'figures/fig_serotrend_colorful.jpg','Resolution',600)

%% ── 15. Seroprevalence time series panel map (IQR + 95% CI) ────
% fprintf('Fig: seroprevalence time series...\n');
% figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
% tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','tight')
% for si = 1:S
%     tno = state_to_tile(si); if tno==0; continue; end
%     nexttile(tno)
%     seropos = squeeze(1 - mean(op_seroposall_0(:,si,:,:),4));  % [n_draw x T]
%     sero_s  = sort(seropos);
% 
%     % yyaxis right
%     % plot(years, inc_mat(:,si), 'k:', 'LineWidth',1.2)
%     % ylim([0 inc_max]); ax=gca; ax.YAxis(2).Color='k';
%     % yyaxis left
%     fill([years, fliplr(years)], [sero_s(maxidx,:), fliplr(sero_s(minidx,:))], ...
%          col_report_light,'EdgeColor','none','FaceAlpha',0.5); hold on
%     fill([years, fliplr(years)], [sero_s(q75idx,:), fliplr(sero_s(q25idx,:))], ...
%          col_report,'EdgeColor','none','FaceAlpha',0.5);
%     plot(years, median(sero_s), 'color',col_report,'LineWidth',1.2)
%     ylim([0 1]); ax.YAxis(1).Color=col_report;
%     xlim([years(1)-0.5, years(end)+0.5])
%     xticks([2000 2010 2020]); xticklabels({'00','10','20'})
%     set(gca,'FontName','SansSerif','FontSize',8)
%     title(state_names{si},'FontSize',8,'FontWeight','normal')
% end
% nexttile(legend_tile); axis off; hold on
% h1=fill(nan,nan,col_report_light,'FaceAlpha',0.5,'EdgeColor','none');
% h2=fill(nan,nan,col_report,'FaceAlpha',0.5,'EdgeColor','none');
% h3=plot(nan,nan,'-','color',col_report,'LineWidth',1.5);
% % h4=plot(nan,nan,'k:','LineWidth',1.5);
% legend([h3 h2 h1],{'Median seropos','IQR (25-75%)','95% CrI'}, ...
%     'Location','northwest','FontSize',11); legend('boxoff')
% exportgraphics(gcf,'figures/fig_seropos_map.jpg','Resolution',600)

%% ── 16. 주별 연도별 케이스 인시던스 패널맵 (PPC model vs data, IQR+95%CI) ──
% fprintf('Fig: state-year case incidence panel map (model vs data)...\n');
% model_med  = squeeze(median(model_inc, 1));               % [T x S]
% model_lo50 = squeeze(quantile(model_inc, 0.25, 1));
% model_hi50 = squeeze(quantile(model_inc, 0.75, 1));
% model_lo95 = squeeze(quantile(model_inc, 0.025, 1));
% model_hi95 = squeeze(quantile(model_inc, 0.975, 1));
% 
% figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
% tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','tight')
% for si = 1:S
%     tno = state_to_tile(si); if tno==0; continue; end
%     nexttile(tno)
%     % 95% CI - 연하게
%     fill([years, fliplr(years)], [model_hi95(:,si)', fliplr(model_lo95(:,si)')], ...
%          col_model_light, 'EdgeColor','none','FaceAlpha',0.5); hold on
%     % IQR - 진하게
%     fill([years, fliplr(years)], [model_hi50(:,si)', fliplr(model_lo50(:,si)')], ...
%          col_model, 'EdgeColor','none','FaceAlpha',0.5);
%     plot(years, model_med(:,si), 'color', col_model, 'LineWidth', 1.3)
%     scatter(years, data_inc_bystate(:,si), 16, col_data, 'square', 'filled')
% 
%     xlim([years(1)-0.5, years(end)+0.5])
%     ylim([0 max([data_inc_bystate(:); model_hi95(:)])*1.05])
%     xticks([2000 2010 2020]); xticklabels({'00','10','20'})
%     set(gca,'FontName','SansSerif','FontSize',8)
%     title(state_names{si},'FontSize',8,'FontWeight','normal')
% end
% nexttile(legend_tile); axis off; hold on
% h1 = plot(nan,nan,'-','color',col_model,'LineWidth',1.5);
% h2 = fill(nan,nan,col_model,'FaceAlpha',0.5,'EdgeColor','none');
% h3 = fill(nan,nan,col_model_light,'FaceAlpha',0.5,'EdgeColor','none');
% h4 = scatter(nan,nan,16,col_data,'square','filled');
% legend([h1 h2 h3 h4], {'Model median (PPC)','IQR (25-75%)','95% CI','Observed (data)'}, ...
%     'Location','northwest','FontSize',11); legend('boxoff')
% sgtitle('Reported case incidence by state and year: model (PPC) vs data','FontSize',11,'FontName','SansSerif')
% exportgraphics(gcf,'figures/fig_case_incidence_bystate.jpg','Resolution',600)

%% ── 17. 주별 누적 케이스 인시던스(연령별, 분모=총인구, per 100) 패널맵 ──
fprintf('Fig: cumulative case incidence per 100 total population by age group (model vs data)...\n');

% 모델: cum_cases_age [n_draw x S x A5] (절대 케이스수) / 주 전체 누적 인구 * 100
cum_inc_age_model = zeros(n_draw, S, A5);
for s = 1:S
    cum_inc_age_model(:,s,:) = squeeze(cum_cases_age(:,s,:)) / data_pop_total_bystate(s) * 100000;
end

model_age_med  = squeeze(median(cum_inc_age_model, 1));      % [S x A5]
model_age_lo50 = squeeze(quantile(cum_inc_age_model, 0.25, 1));
model_age_hi50 = squeeze(quantile(cum_inc_age_model, 0.75, 1));
model_age_lo95 = squeeze(quantile(cum_inc_age_model, 0.025, 1));
model_age_hi95 = squeeze(quantile(cum_inc_age_model, 0.975, 1));

% 데이터: 절대 케이스수 / 주 전체 누적 인구 * 100
data_inc_age_bystate = data_case_age_bystate ./ data_pop_total_bystate * 100000;  % [S x A5]

figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8,NCOLS,'TileSpacing','tight','Padding','tight')
for si = 1:S
    tno = state_to_tile(si); if tno==0; continue; end
    nexttile(tno)
    fill([agerangeidx, fliplr(agerangeidx)], ...
         [model_age_hi95(si,:), fliplr(model_age_lo95(si,:))], ...
         col_model_light, 'EdgeColor','none','FaceAlpha',0.5); hold on
    fill([agerangeidx, fliplr(agerangeidx)], ...
         [model_age_hi50(si,:), fliplr(model_age_lo50(si,:))], ...
         col_model, 'EdgeColor','none','FaceAlpha',0.5);
    plot(agerangeidx, model_age_med(si,:), 'o-', 'color', col_model, ...
         'MarkerFaceColor', col_model, 'MarkerSize', 2.5, 'LineWidth', 1)
    bar(agerangeidx, data_inc_age_bystate(si,:), 0.6, 'FaceColor', col_data, ...
        'EdgeColor','none', 'FaceAlpha', 0.55)

    xlim([-2.5 82.5])
    % ylim([0 max([data_inc_age_bystate(:); model_age_hi95(:)])*1.1])
    xticks([0 20 40 60 80])
    set(gca,'FontName','SansSerif','FontSize',8)
    title(state_names{si},'FontSize',8,'FontWeight','normal')
end
nexttile(legend_tile); axis off; hold on
h1 = plot(nan,nan,'o-','color',col_model,'MarkerFaceColor',col_model,'LineWidth',1);
h2 = fill(nan,nan,col_model,'FaceAlpha',0.5,'EdgeColor','none');
h3 = fill(nan,nan,col_model_light,'FaceAlpha',0.5,'EdgeColor','none');
h4 = bar(nan, nan, 'FaceColor', col_data, 'EdgeColor','none','FaceAlpha',0.55);
legend([h1 h2 h3 h4], {'Model median','IQR (25-75%)','95% CrI','Observed (data)'}, ...
    'Location','northwest','FontSize',11); legend('boxoff')
% sgtitle('Cumulative case contribution per 100 total population by age group, by state (2000–2025): model vs data', ...
%     'FontSize',11,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_case_age_dist_bystate.jpg','Resolution',600)

%% ── 19. FOI → R₀ 변환 + GPR 미래 추정 (IQR+95%CI) ──────────────
% lambda_adj는 이미 kappa 평균으로 정규화되어 있으므로 추가 보정 불필요
L_exp    = 76.4;
t_hist   = 1:T;
t_fut    = (T+1):(T+10);
xlab_all = [years, (years(end)+1):(years(end)+10)];

gpr_opts_lin = {'Basis','linear','FitMethod','exact','PredictMethod','exact'};

R0pred_q25 = []; R0pred_q50 = []; R0pred_q75 = [];

figure('Units','centimeters','Position',[0 0 18 22],'PaperPositionMode','auto');
tiledlayout(8, NCOLS, 'TileSpacing','tight','Padding','tight')
for si = 1:S
    tno = state_to_tile(si); if tno == 0; continue; end

    lam_si  = squeeze(lambda_adj(:,si,:));        % [n_draw x T], 이미 보정됨
    R_draws = lam_si * L_exp + 1;                  % [n_draw x T]

    q25_hist = quantile(R_draws, 0.25, 1);
    q50_hist = quantile(R_draws, 0.50, 1);
    q75_hist = quantile(R_draws, 0.75, 1);

    try
        mdl50 = fitrgp(t_hist', q50_hist', gpr_opts_lin{:});
        mdl25 = fitrgp(t_hist', q25_hist', gpr_opts_lin{:});
        mdl75 = fitrgp(t_hist', q75_hist', gpr_opts_lin{:});
        q50_fut = max(1, predict(mdl50, t_fut'));
        q25_fut = max(1, predict(mdl25, t_fut'));
        q75_fut = max(1, predict(mdl75, t_fut'));
    catch
        p50 = polyfit(t_hist(end-5:end), q50_hist(end-5:end), 1);
        p25 = polyfit(t_hist(end-5:end), q25_hist(end-5:end), 1);
        p75 = polyfit(t_hist(end-5:end), q75_hist(end-5:end), 1);
        q50_fut = max(1, polyval(p50, t_fut));
        q25_fut = max(1, polyval(p25, t_fut));
        q75_fut = max(1, polyval(p75, t_fut));
    end

    q50_all = [q50_hist, q50_fut'];
    q25_all = [q25_hist, q25_fut'];
    q75_all = [q75_hist, q75_fut'];
    R0pred_q50 = [R0pred_q50; q50_all];
    R0pred_q25 = [R0pred_q25; q25_all];
    R0pred_q75 = [R0pred_q75; q75_all];

    nexttile(tno)
    fill([xlab_all(1:T), fliplr(xlab_all(1:T))], [q75_hist, fliplr(q25_hist)], ...
         col_foi, 'EdgeColor','none','FaceAlpha',0.3); hold on
    fill([xlab_all(T:end), fliplr(xlab_all(T:end))], ...
         [[q75_hist(end), q75_fut'], fliplr([q25_hist(end), q25_fut'])], ...
         col_foi, 'EdgeColor','none','FaceAlpha',0.15)
    plot(xlab_all(1:T),   q50_hist, '-',  'color',col_foi, 'LineWidth',1.5)
    plot(xlab_all(T:end), [q50_hist(end), q50_fut'], ':', 'color',col_foi, 'LineWidth',1.5)
    xline(years(end)+0.5, 'k:', 'LineWidth',0.8)

    xlim([years(1), xlab_all(end)]); ylim([1 inf])
    xticks([2000 2010 2020 2030]); xticklabels({'00','10','20','30'})
    set(gca,'FontName','SansSerif','FontSize',8)
    title(state_names{si},'FontSize',8,'FontWeight','normal')
end
nexttile(legend_tile); axis off; hold on
h1 = plot(nan,nan,'-', 'color',col_foi,'LineWidth',1.5);
h2 = plot(nan,nan,'--','color',col_foi,'LineWidth',1.5);
h3 = fill(nan,nan,col_foi,'FaceAlpha',0.3,'EdgeColor','none');
legend([h1 h2 h3], {'Median R₀ (estimated)','Median R₀ (projected)','IQR'}, ...
    'Location','northwest','FontSize',11); legend('boxoff')
% sgtitle('R₀ projection (FOI × 76.4 + 1)', 'FontSize',11,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_R0_projection.jpg','Resolution',600)

save('results/R0predict_new.mat','R0pred_q25','R0pred_q50','R0pred_q75', ...
     'state_codes','xlab_all','T')

%% Case/Infection ratio by state (catalytic model draws 기반)
%% Case/Infection ratio by state (catalytic model draws 기반, inf_total 직접 사용)
%% 0. 기본 설정
S = 27; T = 26; A5 = 17;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};

north_idx     = find(ismember(state_codes,{'AC','AM','AP','PA','RO','RR','TO'}));
northeast_idx = find(ismember(state_codes,{'AL','BA','CE','MA','PB','PE','PI','RN','SE'}));
centwest_idx  = find(ismember(state_codes,{'DF','GO','MS','MT'}));
southeast_idx = find(ismember(state_codes,{'ES','MG','RJ','SP'}));
south_idx     = find(ismember(state_codes,{'PR','RS','SC'}));

col_region = [1.000 0.498 0.055;   % North
              0.596 0.306 0.639;   % Northeast
              0.090 0.745 0.812;   % Central-west
              0.173 0.627 0.173;   % Southeast
              0.549 0.341 0.271];  % South

% 지역순 재배열 (첨부 이미지와 동일한 순서)
north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];
reg_order = [north_order, northeast_order, centwest_order, ...
             southeast_order, south_order];
region_of = zeros(S,1);
region_of(north_idx)=1; region_of(northeast_idx)=2; region_of(centwest_idx)=3;
region_of(southeast_idx)=4; region_of(south_idx)=5;

fprintf('Loading data...\n')
load('results/catalytic_ppc_draws.mat', 'cum_cases_age', 'inf_total');
% inf_total: [n_draw x S x T x A5] (to_matlab_order 적용됨, rho/omega/phi 미적용 총감염)
n_draw = size(cum_cases_age,1);

cum_inf_age = squeeze(sum(inf_total, 3));   % [n_draw x S x A5]

% cum_cases_age: [n_draw x S x A5] (PPC 기반 기대 보고 케이스, 누적)
% cum_inf_age:   [n_draw x S x A5] (모델 총감염 기대값, 누적, 보고율 미적용)
ratio_draws = sum(cum_cases_age, 3) ./ sum(cum_inf_age, 3);  % [n_draw x S]

ratio_med = median(ratio_draws, 1);          % [1 x S]
ratio_lo  = quantile(ratio_draws, 0.025, 1);
ratio_hi  = quantile(ratio_draws, 0.975, 1);
ratio_q25 = quantile(ratio_draws, 0.25, 1);
ratio_q75 = quantile(ratio_draws, 0.75, 1);

% 콘솔 출력
fprintf('\n--- Case/Infection ratio by state (q50, 95%% CI) ---\n')
for si = 1:S
    fprintf('  %s: %.3f  [%.3f, %.3f]\n', state_codes{si}, ...
        ratio_med(si), ratio_lo(si), ratio_hi(si))
end
fprintf('  National (median of state medians): %.3f\n', median(ratio_med))

mkdir('figures')
figure('Units','centimeters','Position',[0 0 20 7],'PaperPositionMode','auto');
hold on

% 지역 배경
reg_bound = [0, cumsum([numel(north_order), numel(northeast_order), ...
    numel(centwest_order), numel(southeast_order), numel(south_order)])];
ymax_bg = max(ratio_hi)*1.3;
for r = 1:5
    fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
         [0 0 ymax_bg ymax_bg], col_region(r,:),'EdgeColor','none','FaceAlpha',0.15)
end

% 지역순 재배열
r_med = ratio_med(reg_order);
r_lo  = ratio_lo(reg_order);
r_hi  = ratio_hi(reg_order);
r_q25 = ratio_q25(reg_order);
r_q75 = ratio_q75(reg_order);
reg_c = region_of(reg_order);

% IQR 캡 없는 에러바 (95% CI, 연하게)
for si = 1:S
    c = col_region(reg_c(si),:);
    % 95% CI 선
    % plot([si si], [r_lo(si) r_hi(si)], '-', 'Color', c*0.6+0.4, 'LineWidth', 1.0)
    % IQR 굵게
    % plot([si si], [r_q25(si) r_q75(si)], '-', 'Color', c, 'LineWidth', 1.5)
    plot([si si], [r_lo(si) r_hi(si)], '-', 'Color', c, 'LineWidth', 1.5)
    % 중앙값 마커
    plot(si, r_med(si), 's', 'MarkerFaceColor', c, 'MarkerEdgeColor', c*0.7, ...
         'MarkerSize', 7, 'LineWidth', 0.5)
end

codes_reord = state_codes(reg_order);
xlim([0.25 S+0.75]); ylim([0 ymax_bg])
ylim([0, 0.15])
xticks(1:S); xticklabels(codes_reord)
set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',12)
ylabel('Case / Infection ratio')
% yline(median(ratio_med), 'k:', 'LineWidth', 1.0)   % 전국 중앙값 참고선
box on

exportgraphics(gcf,'figures/fig_case_infection_ratio.jpg','Resolution',600)
fprintf('Saved: figures/fig_case_infection_ratio.jpg\n')
%% ── FOI 시계열 (선택 3개 주, pre-2000 포함) ────────────────────
% ▼ 여기서 주 번호 선택 (알파벳순 인덱스, 1~27)
%   5=BA(Bahia), 26=SP(São Paulo), 24=SC(Santa Catarina)
sel_idx   = [5, 26, 24];
sel_names = state_names(sel_idx);
sel_cols  = {'r', 'b', 'k'};
sel_lsty  = {'-', '-', '--'};

figure('Units','centimeters','Position',[0 0 7 6],'PaperPositionMode','auto');
hold on

xlab_all = ['P00', cellstr(num2str((0:T-1)', '%02d'))'];  % P00, 00, 01, ..., 25

for k = 1:3
    si = sel_idx(k);

    % pre-2000 FOI median (lambda_pre_adj)
    foi_pre = median(lambda_pre_adj(:, si)) * 4;   % ×4: 4 serotypes

    % 연도별 FOI median [1×T]
    foi_t = zeros(1, T);
    for ti = 1:T
        foi_t(ti) = median(squeeze(lambda_adj(:, si, ti))) * 4;
    end

    % x=1: P00, x=2..T+1: 연도
    % plot(1, foi_pre, 'o', 'Color', sel_cols{k}, 'MarkerFaceColor', sel_cols{k}, 'MarkerSize', 4)
    % plot(2:T+1, foi_t, sel_lsty{k}, 'Color', sel_cols{k}, 'LineWidth', 1.5, ...
    %     'DisplayName', sel_names{k})
    plot(1:T+1, 4*[foi_pre foi_t], sel_lsty{k}, 'Color', sel_cols{k}, 'LineWidth', 1.5, ...
        'DisplayName', sel_names{k})
end

% P00 구분선
xline(1.5, 'k:', 'LineWidth', 0.8)

xticks(1:2:T+1)
xticklabels(xlab_all(1:2:end))
xlim([0.5, T+1.5])
ylim([0 inf])
ylabel('Force of infection (FOI)')

% legend (연속선 마커용)
hleg = gobjects(3,1);
for k = 1:3
    hleg(k) = plot(nan, nan, sel_lsty{k}, 'Color', sel_cols{k}, 'LineWidth', 1.5);
end
legend(hleg, sel_names, 'Location','northeast','FontSize',9); legend('boxoff')
set(gca,'FontName','SansSerif','FontSize',11)

exportgraphics(gcf,'figures/fig_FOI_selected_states.jpg','Resolution',600)

%% ── 18. National summary statistics for manuscript text (Table S7 discussion) ──
% Fully posterior-based, population-weighted national estimates.
% Addresses reviewer request (Henrik, Comment 26): proportion of population
% with >=1 prior infection in 2000 vs 2025, and average annual infections
% per capita. All quantities computed per posterior draw, then summarized
% as median + 95% CrI (no reliance on point estimates/medians pre-aggregation).

fprintf('Computing national summary statistics (posterior-based)...\n');

% inf_total: [n_draw x S x T x A5], modeled NEW infections per year (not cumulative)
if ~exist('inf_total','var')
    load('results/catalytic_ppc_draws.mat','inf_total');
end

% ---- (a) Population by state, year 2000 (t=1) and year 2025 (t=T) ----
pop_age_2000 = zeros(S, A5);
pop_age_2025 = zeros(S, A5);
for s = 1:S
    idx0 = (s-1)*T + 1;    % t = 2000
    idxE = (s-1)*T + T;    % t = 2025
    pop_age_2000(s,:) = pop_mat(idx0,:);
    pop_age_2025(s,:) = pop_mat(idxE,:);
end
pop_state_2000 = sum(pop_age_2000, 2);
pop_state_2025 = sum(pop_age_2025, 2);
natl_pop_2000  = sum(pop_state_2000);
natl_pop_2025  = sum(pop_state_2025);

% ---- (b) Population-weighted proportion with >=1 prior infection, per draw ----
prop_atleast1_2000_draws = zeros(n_draw,1);
prop_atleast1_2025_draws = zeros(n_draw,1);

for d = 1:n_draw
    % --- 2000 baseline ---
    x0_2000 = squeeze(op_seroposall_pre_0(d,:,:));        % [S x A5], proportion naive
    naive_count_2000 = sum(x0_2000 .* pop_age_2000, 2);   % [S x 1]
    prop_atleast1_2000_draws(d) = 1 - sum(naive_count_2000) / natl_pop_2000;

    % --- 2025 (ti = T) ---
    x0_2025 = squeeze(op_seroposall_0(d,:,T,:));          % [S x A5]
    naive_count_2025 = sum(x0_2025 .* pop_age_2025, 2);
    prop_atleast1_2025_draws(d) = 1 - sum(naive_count_2025) / natl_pop_2025;
end

fprintf('\n--- National %% with >=1 prior infection ---\n');
fprintf('2000: median %.1f%%  [95%% CrI %.1f, %.1f]\n', ...
    100*median(prop_atleast1_2000_draws), 100*pct025(prop_atleast1_2000_draws), 100*pct975(prop_atleast1_2000_draws));
fprintf('2025: median %.1f%%  [95%% CrI %.1f, %.1f]\n', ...
    100*median(prop_atleast1_2025_draws), 100*pct025(prop_atleast1_2025_draws), 100*pct975(prop_atleast1_2025_draws));

% ---- (c) Average annual infections per capita, 2000-2025 (model-based, not finite-difference) ----
annual_infections_percapita_draws = zeros(n_draw,1);
pop_yr = zeros(T,1);
for t = 1:T
    pop_yr(t) = sum(pop_mat((0:S-1)'*T + t, :), 'all');  % total national pop that year
end
for d = 1:n_draw
    inf_yr = squeeze(sum(sum(inf_total(d,:,:,:),2),4));  % [T x 1], national infections per year
    annual_infections_percapita_draws(d) = mean(inf_yr ./ pop_yr);
end

fprintf('\n--- Average annual infections per capita (2000-2025) ---\n');
fprintf('median %.4f  [95%% CrI %.4f, %.4f]\n', ...
    median(annual_infections_percapita_draws), pct025(annual_infections_percapita_draws), pct975(annual_infections_percapita_draws));

% ---- (d) Save for manuscript text pull ----
summary_stats = struct( ...
    'prop_atleast1_2000_median', median(prop_atleast1_2000_draws), ...
    'prop_atleast1_2000_lo95',   pct025(prop_atleast1_2000_draws), ...
    'prop_atleast1_2000_hi95',   pct975(prop_atleast1_2000_draws), ...
    'prop_atleast1_2025_median', median(prop_atleast1_2025_draws), ...
    'prop_atleast1_2025_lo95',   pct025(prop_atleast1_2025_draws), ...
    'prop_atleast1_2025_hi95',   pct975(prop_atleast1_2025_draws), ...
    'annual_inf_percapita_median', median(annual_infections_percapita_draws), ...
    'annual_inf_percapita_lo95',   pct025(annual_infections_percapita_draws), ...
    'annual_inf_percapita_hi95',   pct975(annual_infections_percapita_draws));
save('results/national_summary_stats.mat', 'summary_stats');
fprintf('\nSaved: results/national_summary_stats.mat\n');
%% ── Model calibration (incidence version): observed vs. predicted, by state ──
% Same as above, but divided by cumulative 2000-2025 population per state
% (data_pop_total_bystate, from Section 1) to express as incidence per
% 100,000 — consistent with the scale already used in Section 17
% (cum_inc_age_model / data_inc_age_bystate).

% ---- Observed: total cases per state (sum over age groups) ----
obs_cases = sum(data_case_age_bystate, 2);   % [S x 1]
fprintf('Computing state-level observed vs. model INCIDENCE (per 100,000, 2000-2025)...\n');
model_cases_draws = sum(cum_cases_age, 3);   % [n_draw x S]


% obs_cases, model_cases_draws, reg_order, region_of, col_region, state_codes
% already computed/defined above.

obs_incidence_pc   = obs_cases ./ data_pop_total_bystate * 100000;          % [S x 1]
model_incidence_pc = model_cases_draws ./ data_pop_total_bystate' * 100000; % [n_draw x S]

model_inc_med  = mean(model_incidence_pc, 1);
model_inc_lo95 = quantile(model_incidence_pc, 0.025, 1);
model_inc_hi95 = quantile(model_incidence_pc, 0.975, 1);

% reorder by region
obs_inc_reord   = obs_incidence_pc(reg_order);
model_inc_med_r = model_inc_med(reg_order);
model_inc_lo95_r = model_inc_lo95(reg_order);
model_inc_hi95_r = model_inc_hi95(reg_order);

figure('Units','centimeters','Position',[0 0 17 5],'PaperPositionMode','auto');
hold on

ymax_bg = max([obs_inc_reord(:); model_inc_hi95_r(:)]) * 1.15;

for r = 1:5
    fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
         [0 0 ymax_bg ymax_bg], col_region(r,:), 'EdgeColor','none','FaceAlpha',0.15)
end

bar(1:S, obs_inc_reord, 0.6, 'FaceColor', col_data, 'EdgeColor','none', 'FaceAlpha', 0.75)

for si = 1:S
    plot([si si], [model_inc_lo95_r(si) model_inc_hi95_r(si)], 'k-', 'LineWidth', 1.2)
    plot(si, model_inc_med_r(si), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 5, 'LineWidth', 0.5)
end

xlim([0.25 S+0.75]); ylim([0 ymax_bg])
xticks(1:S); xticklabels(codes_reord)
set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',11)
% ylabel('Cumulative incidence')
% ylabel({'Average annual incidence', ['(per 100,000, 2000' char(8211) '2025)']})
ylabel({'Annual incidence', ['(2000' char(8211) '2025)']})
box on

h1 = fill(nan,nan,col_data,'FaceAlpha',0.75,'EdgeColor','none');
h2 = plot(nan,nan,'ko','MarkerFaceColor','k','MarkerSize',5,'LineWidth',0.5);
h3 = plot(nan,nan,'k-','LineWidth',1.2);
legend([h1 h2 h3], {'Observed (data)','Model mean','Model 95% CrI'}, ...
    'Location','north','FontSize',9); legend('boxoff')

exportgraphics(gcf,'figures/fig_model_calibration_bystate_incidence.jpg','Resolution',600)
copygraphics(gcf, 'ContentType', 'vector')
fprintf('Saved: figures/fig_model_calibration_bystate_incidence.jpg\n')

%% ── Brazil region map (Fig. 3: state colors = region colors) ──
% Section 0의 col_region, *_codes를 그대로 써서 다른 모든 그림과 색/지역 분류가 동일.
% 필요: 현재 폴더에 BR_UF_2023.shp (IBGE), Section 0 실행된 상태.

label_offset = containers.Map( ...
    {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
     'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
     'RO','RR','RS','SC','SE','SP','TO'}, ...
    {[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0.8],[0,0],[0,0],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,-0.8],[0,0],[0,0],[0,-0.8],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]});

Shfile = shaperead("BR_UF_2023.shp");

% 라벨용 centroid
cxvec = zeros(1,length(Shfile));
cyvec = zeros(1,length(Shfile));
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    valid = ~isnan(x) & ~isnan(y);
    pg = polyshape(x(valid),y(valid));
    [cx,cy] = centroid(pg);
    cxvec(i) = cx; cyvec(i) = cy;
end
cxvec(26) = cxvec(26) - 1;   % 기존 코드 그대로 (Shfile 26번째 주 라벨 보정)

figure('Units','centimeters','Position',[0 0 10 12],'PaperPositionMode','auto');
hold on
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    code = Shfile(i).SIGLA_UF;
    if     ismember(code, north_codes);     fc = col_region(1,:);
    elseif ismember(code, northeast_codes); fc = col_region(2,:);
    elseif ismember(code, centwest_codes);  fc = col_region(3,:);
    elseif ismember(code, southeast_codes); fc = col_region(4,:);
    elseif ismember(code, south_codes);     fc = col_region(5,:);
    else;                                   fc = [0.8 0.8 0.8];
    end
    nan_idx = find(isnan(x));
    starts = [1, nan_idx+1]; ends = [nan_idx-1, length(x)];
    for j = 1:length(starts)
        seg = starts(j):ends(j);
        if length(seg) >= 3
            patch(x(seg),y(seg),fc,'EdgeColor','w','LineWidth',0.8,'FaceAlpha',1)
        end
    end
end
for i = 1:length(Shfile)
    off = label_offset(Shfile(i).SIGLA_UF);
    text(cxvec(i)+off(1), cyvec(i)+off(2), Shfile(i).SIGLA_UF, ...
        'HorizontalAlignment','center','FontName','SansSerif')
end
axis equal off
set(gca,'FontName','SansSerif','FontSize',14)
exportgraphics(gcf,'figures/fig_map_region.jpg','Resolution',600)