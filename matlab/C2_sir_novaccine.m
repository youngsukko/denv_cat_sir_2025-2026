%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — SIR simulation + BoD (no vaccine baseline)
%  Requires: func_ode_alt_v5.m (unchanged)
%             results/catalytic_ppc_draws.mat
%             results/catalytic_lambda_kappa_draws.mat
%             results/R0predict_new.mat
%             results/pop_age_state_year.csv
%             results/cases_age_state_year.csv
%             data_lifeexp.mat, data_survrate.mat, data_birthrate.mat
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};
state_names = {'Acre','Alagoas','Amazonas','Amapá','Bahia','Ceará', ...
               'Distrito Federal','Espírito Santo','Goiás','Maranhão', ...
               'Minas Gerais','Mato Grosso do Sul','Mato Grosso','Pará', ...
               'Paraíba','Pernambuco','Piauí','Paraná','Rio de Janeiro', ...
               'Rio Grande do Norte','Rondônia','Roraima','Rio Grande do Sul', ...
               'Santa Catarina','Sergipe','São Paulo','Tocantins'};

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

age_labels = {'0-','5-','10-','15-','20-','25-','30-','35-', ...
              '40-','45-','50-','55-','60-','65-','70-','75-','80+'};

%% ── 1. 데이터 로드 (전부 새 소스, 재계산 없음) ──────────────────
fprintf('Loading data...\n')
load('results/catalytic_ppc_draws.mat', ...
    'op_seroposall_0','op_seroposall_1','op_seroposall_2', ...
    'op_seroposall_3','op_seroposall_4');
n_draw = size(op_seroposall_0,1);

load('results/R0predict_new.mat','R0pred_q50','R0pred_q25','R0pred_q75');
% R0pred_q50 [S × (T+10)]

load('data_lifeexp.mat','data_lifeexp')
data_lifeexp_age = data_lifeexp(1:5:end)';  % [1×17], 5세 밴드

load('results/catalytic_lambda_kappa_draws.mat', 'kappa_draws','omega_draws','delta_draws');
% kappa_age 평균 정규화 (R0 projection 때와 동일한 물리량 기준 통일)
mean_kappa  = 1;          % [n_draw x 1]
kappa_adj   = kappa_draws ./ mean_kappa;     % [n_draw x A5], 평균=1

% ── pop_age_state_year.csv 연령 컬럼 순서 버그 보정 ──────────────
% (R의 sort(unique(age_group)) 알파벳정렬 때문에 "5-9"가 45-49와 50-54 사이에 끼어듦)
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
assert(all(age_reorder_idx > 0), '연령 컬럼 매핑 실패');

pop_mat = table2array(readtable('results/pop_age_state_year.csv'));  % [S*T × A5]
pop_mat = pop_mat(:, age_reorder_idx);   % 올바른 연령순으로 재배열

load('data_survrate.mat'); load('data_birthrate.mat')
data_birthper1K_since1960(end:end+50) = data_birthper1K_since1960(end);
survrate = survrate(1:5:end)';

%% ── 2. 공통 inputs 설정 ──────────────────────────────────────
inputs = struct();
inputs.agenum   = A5;
inputs.beta0    = median(kappa_adj)';          % [A5×1] age-FOI multiplier (정규화됨, 평균=1)
inputs.cageadj  = median(omega_draws)';        % [A5×1] age-reporting multiplier

inputs.cprop_1  = median(delta_draws(:,1));
inputs.cprop_2  = median(delta_draws(:,2));
inputs.cprop_3  = median(delta_draws(:,3));
inputs.cprop_4  = median(delta_draws(:,4));

inputs.gamma    = 1/(1/52);    % 1주 감염기간 (연 단위)
inputs.ohm      = 1/1;         % 1년 일시면역
inputs.i1       = 3/4;
inputs.i2       = 2/4;
inputs.i3       = 1/4;
inputs.mu       = 1/(5*1);     % 5년 연령이동 (연 단위)
inputs.mu_final = 1/(8.9*1);   % 80+ 사망률

inputs.q        = survrate(1:16)';
inputs.birthrate = data_birthper1K_since1960(64);

% 백신 효과 파라미터는 유지하되, 이번 시나리오에서는 vacc=0이라 사용되지 않음
inputs.e_n_d=0.562;  inputs.e_n_i=0.281;  inputs.e_n_h=0.844; inputs.e_n_m=0.844;
inputs.e_p_d=0.720;  inputs.e_p_i=0.360;  inputs.e_p_h=0.844; inputs.e_p_m=0.844;
inputs.ew_n_d=0.405; inputs.ew_n_i=0.2025; inputs.ew_n_h=0.844; inputs.ew_n_m=0.844;
inputs.ew_p_d=0.531; inputs.ew_p_i=0.2655; inputs.ew_p_h=0.844; inputs.ew_p_m=0.844;
inputs.eww_n_d=0.233; inputs.eww_n_i=0.1165; inputs.eww_n_h=0.844; inputs.eww_n_m=0.844;
inputs.eww_p_d=0.319; inputs.eww_p_i=0.1595; inputs.eww_p_h=0.844; inputs.eww_p_m=0.844;
inputs.ohm_vac=1/(1); inputs.ohm_vac2=1/(4); inputs.ohm_vac3=0;

h_ref = [0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
         0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
f_ref = [0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
         0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';
inputs.h1=h_ref; inputs.h2=h_ref; inputs.h3=h_ref; inputs.h4=h_ref;
inputs.f1=f_ref; inputs.f2=f_ref; inputs.f3=f_ref; inputs.f4=f_ref;

tspan = 0:4;   % 5년 시뮬레이션 (연 단위)

%% ── 3. 시나리오 설정: NO VACCINE BASELINE ───────────────────
vacc_scenario = 0;   % 0 = no vaccination
vacorder      = 1;   % vacc=0일 때는 사용되지 않음 (placeholder)

%% ── 4. 시뮬레이션 (q25/q50/q75 각각) ───────────────────────
R0_scenarios   = {R0pred_q25, R0pred_q50, R0pred_q75};
scenario_names = {'q25','q50','q75'};

Case_age_3  = zeros(S, A5, 3);
Hosp_age_3  = zeros(S, A5, 3);
Dead_age_3  = zeros(S, A5, 3);
DALY_age_3  = zeros(S, A5, 3);
popbystate  = zeros(S, 1);
riskpop     = zeros(S, 1);

for sc = 1:3
    R0pred_sc = R0_scenarios{sc};
    fprintf('Scenario %s...\n', scenario_names{sc})

    for si = 1:S
        pop_2025 = pop_mat((si-1)*T + T, :)';

        A0 = median(squeeze(op_seroposall_0(:,si,T,:)),1)' .* pop_2025;
        A1 = median(squeeze(op_seroposall_1(:,si,T,:)),1)' .* pop_2025;
        A2 = median(squeeze(op_seroposall_2(:,si,T,:)),1)' .* pop_2025;
        A3 = median(squeeze(op_seroposall_3(:,si,T,:)),1)' .* pop_2025;
        A4 = median(squeeze(op_seroposall_4(:,si,T,:)),1)' .* pop_2025;
        popset = [A0,A1,A2,A3,A4];

        if sc == 2   % q50일 때만 저장
            popbystate(si) = sum(popset,'all');
            riskpop(si)    = sum(popset(:,1)+popset(:,2),'all')/sum(popset,'all');
        end

        R0mat_si      = R0pred_sc(si, T);
        inputs.beta1  = R0pred_sc(si, T+1:end) * 52;
        inputs.vacc   = vacc_scenario;
        inputs.vacorder = vacorder;

        mmu = 1/76.4;
        y0  = zeros(51*A5, 1);
        y0(A5*0+(1:A5)) = A0 - mmu/(mmu+52)*(1-1/R0mat_si)*A0;
        y0(A5*3+(1:A5)) = A1 - mmu/(mmu+52)*(1-1/R0mat_si)*A1;
        y0(A5*6+(1:A5)) = A2 - mmu/(mmu+52)*(1-1/R0mat_si)*A2;
        y0(A5*9+(1:A5)) = A3 - mmu/(mmu+52)*(1-1/R0mat_si)*A3;
        y0(A5*11+(1:A5))= A4;
        y0(A5*1+(1:A5)) = mmu/(mmu+52)*(1-1/R0mat_si)*A0;
        y0(A5*4+(1:A5)) = mmu/(mmu+52)*(1-1/R0mat_si)*A1;
        y0(A5*7+(1:A5)) = mmu/(mmu+52)*(1-1/R0mat_si)*A2;
        y0(A5*10+(1:A5))= mmu/(mmu+52)*(1-1/R0mat_si)*A3;

        [~, y] = ode23(@(t,yy) func_ode_alt_v5(t,yy,inputs), tspan, y0);

        Casen = y(:, A5*48+(1:A5));
        Hosp  = y(:, A5*49+(1:A5));
        Dead  = y(:, A5*50+(1:A5));
        YLD_age = Casen*0.00084 + Hosp*0.0051;
        YLL_age = Dead .* data_lifeexp_age;

        Case_age_3(si,:,sc) = Casen(end,:);
        Hosp_age_3(si,:,sc) = Hosp(end,:);
        Dead_age_3(si,:,sc) = Dead(end,:);
        DALY_age_3(si,:,sc) = YLD_age(end,:) + YLL_age(end,:);
    end
end

% 편의 변수 (q50 = sc=2)
Case_0_age = Case_age_3(:,:,2);
Hosp_0_age = Hosp_age_3(:,:,2);
Dead_0_age = Dead_age_3(:,:,2);
daly_0_age = DALY_age_3(:,:,2);

save('results/SIR_BoD_baseline_novaccine.mat', ...
    'Case_age_3','Hosp_age_3','Dead_age_3','DALY_age_3', ...
    'Case_0_age','Hosp_0_age','Dead_0_age','daly_0_age', ...
    'popbystate','riskpop','state_codes','state_names')
fprintf('Done.\n')

%% ── 5. BoD 시각화 (연령 분포, 지역별) ────────────────────────
north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];

reg_order = [north_order, northeast_order, centwest_order, ...
             southeast_order, south_order];

reg_labels = [repmat(1,1,numel(north_idx)), repmat(2,1,numel(northeast_idx)), ...
              repmat(3,1,numel(centwest_idx)), repmat(4,1,numel(southeast_idx)), ...
              repmat(5,1,numel(south_idx))];
codes_reord = state_codes(reg_order);
pop_reord   = popbystate(reg_order)';

reg_idx  = {north_idx, northeast_idx, centwest_idx, southeast_idx, south_idx};
pop_reg  = cellfun(@(idx) sum(popbystate(idx)), reg_idx);
case_reg = cell2mat(cellfun(@(idx) sum(Case_0_age(idx,:),1)/sum(popbystate(idx))*1e5/5, reg_idx,'UniformOutput',false)');
hosp_reg = cell2mat(cellfun(@(idx) sum(Hosp_0_age(idx,:),1)/sum(popbystate(idx))*1e5/5, reg_idx,'UniformOutput',false)');
dead_reg = cell2mat(cellfun(@(idx) sum(Dead_0_age(idx,:),1)/sum(popbystate(idx))*1e5/5, reg_idx,'UniformOutput',false)');
daly_reg = cell2mat(cellfun(@(idx) sum(daly_0_age(idx,:),1)/sum(popbystate(idx))*1e5/5, reg_idx,'UniformOutput',false)');

fav = 0.5;
reg_bound = [0, cumsum([numel(north_idx), numel(northeast_idx), ...
    numel(centwest_idx), numel(southeast_idx), numel(south_idx)])];

function draw_region_bg(reg_bound, col_region, ymax, fav)
    for r = 1:5
        fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
             [0 0 ymax ymax], col_region(r,:),'EdgeColor','none','FaceAlpha',fav)
        hold on
    end
end

figure('Units','centimeters','Position',[0 0 18 10],'PaperPositionMode','auto');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

nexttile(1)
ppp=plot(0:5:80,(case_reg./sum(case_reg,2))','LineWidth',1.5);
for r=1:5; ppp(r).Color=col_region(r,:); end
xticks(0:5:80); xticklabels(age_labels); xlim([-0.5 80.5])
title('Case proportion by age','FontWeight','normal')
set(gca,'FontName','SansSerif','FontSize',11)

nexttile(2)
ppp=plot(0:5:80,(hosp_reg./sum(hosp_reg,2))','LineWidth',1.5);
for r=1:5; ppp(r).Color=col_region(r,:); end
xticks(0:5:80); xticklabels(age_labels); xlim([-0.5 80.5])
title('Hospitalization proportion by age','FontWeight','normal')
set(gca,'FontName','SansSerif','FontSize',11)

nexttile(3)
ppp=plot(0:5:80,(dead_reg./sum(dead_reg,2))','LineWidth',1.5);
for r=1:5; ppp(r).Color=col_region(r,:); end
xticks(0:5:80); xticklabels(age_labels); xlim([-0.5 80.5])
title('Death proportion by age','FontWeight','normal')
legg=legend('North','Northeast','Central-west','Southeast','South','Location','northwest');
legend('boxoff'); legg.NumColumns=2;
set(gca,'FontName','SansSerif','FontSize',11)

nexttile(4)
ppp=plot(0:5:80,(daly_reg./sum(daly_reg,2))','LineWidth',1.5);
for r=1:5; ppp(r).Color=col_region(r,:); end
xticks(0:5:80); xticklabels(age_labels); xlim([-0.5 80.5])
title('DALY proportion by age','FontWeight','normal')
set(gca,'FontName','SansSerif','FontSize',11)

exportgraphics(gcf,'figures/fig_BoD_age_proportion_novaccine.jpg','Resolution',600)

%% ── 6. 주별 BoD + IQR 범위 ──────────────────────────────────
ylims   = [8000 250 15 400];
ylabels = {'Cases per 100K','Hospitalization per 100K','Deaths per 100K','DALYs per 100K'};
mats_3  = {Case_age_3, Hosp_age_3, Dead_age_3, DALY_age_3};

figure('Units','centimeters','Position',[0 0 18 10],'PaperPositionMode','auto');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for m = 1:4
    nexttile(m)
    for r = 1:5
        fill([reg_bound(r) reg_bound(r+1) reg_bound(r+1) reg_bound(r)]+0.5, ...
             [0 0 ylims(m) ylims(m)], col_region(r,:), ...
             'EdgeColor','none','FaceAlpha',fav); hold on
    end

    v50 = sum(mats_3{m}(:,:,2),2)./popbystate*1e5/5;
    v25 = sum(mats_3{m}(:,:,1),2)./popbystate*1e5/5;
    v75 = sum(mats_3{m}(:,:,3),2)./popbystate*1e5/5;

    v50r = v50(reg_order);
    v25r = v25(reg_order);
    v75r = v75(reg_order);

    bar(v50r, 0.7,'FaceColor','k','EdgeColor',[0.7 0.7 0.7],'LineWidth',0.5)
    errorbar(1:S, v50r, v50r-v25r, v75r-v50r, ...
        'LineStyle','none','Color',[0.3 0.3 0.3],'LineWidth',1.2,'CapSize',3)

    xlim([0.25 S+0.75]); ylim([0 ylims(m)])
    xticks(1:S); xticklabels(codes_reord)
    ylabel(ylabels{m})
    set(gca,'FontName','SansSerif','FontSize',10,'XTickLabelRotation',90)
end

exportgraphics(gcf,'figures/fig_BoD_state_IQR_novaccine.jpg','Resolution',600)

%% ── 7. Scatter: 누적 incidence / riskpop vs DALY IQR ──────────
cases_mat = table2array(readtable('results/cases_age_state_year.csv'));  % [S*T × A5]
cases_mat = cases_mat(:, age_reorder_idx);   % 동일한 연령 순서 버그 보정

data_case_bystate = zeros(T, S);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_case_bystate(:,s) = sum(cases_mat(idx,:), 2);
end

daly_q50 = sum(DALY_age_3(:,:,2),2) ./ popbystate * 1e5 / 5;
daly_q25 = sum(DALY_age_3(:,:,1),2) ./ popbystate * 1e5 / 5;
daly_q75 = sum(DALY_age_3(:,:,3),2) ./ popbystate * 1e5 / 5;

% cum_inc = sum(data_case_bystate(15:24,:), 1)' ./ popbystate * 1e3;  % [S×1]
cum_inc = sum(data_case_bystate(17:26,:), 1)' ./ popbystate * 1e5;  % [S×1]
% cum_inc = sum(data_case_bystate(15:24,:), 1)' ./ popbystate * 1e5;  % [S×1]

region_of = zeros(S,1);
region_of(north_idx)     = 1;
region_of(northeast_idx) = 2;
region_of(centwest_idx)  = 3;
region_of(southeast_idx) = 4;
region_of(south_idx)     = 5;

reg_names = {'North','Northeast','Central-west','Southeast','South'};

figure('Units','centimeters','Position',[0 0 18 7],'PaperPositionMode','auto');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

for panel = 1:2
    nexttile(panel)
    if panel == 1
        xval = cum_inc;
        xl   = 'Past 10-year cumulative incidences (2016-2025)';
    else
        xval = riskpop;
        xl   = 'Proportion ≤1 infection as of 2026';
    end

    p = polyfit(xval, daly_q50, 1);
    xfit = linspace(min(xval)*0.95, max(xval)*1.05, 100);
    plot(xfit, polyval(p,xfit), 'k:', 'LineWidth',1.5); hold on

    h = gobjects(5,1);
    for r = 1:5
        idx = find(region_of == r);
        h(r) = plot(xval(idx), daly_q50(idx), 'o', ...
                    'Color',col_region(r,:), ...
                    'MarkerFaceColor',col_region(r,:), ...
                    'MarkerSize',4,'LineWidth',1);
        errorbar(xval(idx), daly_q50(idx), ...
                 daly_q50(idx)-daly_q25(idx), ...
                 daly_q75(idx)-daly_q50(idx), ...
                 'LineStyle','none','Color',col_region(r,:)*0.7, ...
                 'LineWidth',0.8,'CapSize',3);
    end


    [rho_val, pval] = corr(xval, daly_q50);
    if pval < 0.001
        ptxt = 'p < 0.001';
    elseif pval < 0.01
        ptxt = sprintf('p = %.3f', pval);
    else
        ptxt = sprintf('p = %.2f', pval);
    end
    text(0.05, 0.75, {sprintf('r = %.2f', rho_val); ptxt}, ...
         'Units','normalized','VerticalAlignment','top', ...
         'FontName','SansSerif','FontSize',9)

    xlabel(xl,'FontSize',9)
    ylabel('Annual DALY per 100K, 2026–2030','FontSize',9)
    set(gca,'FontName','SansSerif','FontSize',9)
    box on

    if panel == 1
        legend(h, reg_names,'Location','northeast','FontSize',8); legend('boxoff')
    end
end

exportgraphics(gcf,'figures/fig_DALY_scatter_novaccine.jpg','Resolution',600)

fprintf('\n--- Correlation summary (no vaccine baseline) ---\n')
[r1,p1] = corr(cum_inc, daly_q50);
[r2,p2] = corr(riskpop, daly_q50);
fprintf('Cumul. incidence vs DALY q50:  r = %.3f,  p = %.4f\n', r1, p1)
fprintf('Risk pop (≤1 inf) vs DALY q50: r = %.3f,  p = %.4f\n', r2, p2)

%% ── 8. Fig. 3 inset: 지역 색상 지도 (Mapping Toolbox) ─────────────────────
% 위 섹션들이 실행된 상태여야 함 (state_codes, region_of, col_region 사용)
% 지도를 "최종 크기(cm)"로 그려서 100%로 얹으면 글씨 크기가 그대로 유지됨.
map_w_cm = 3.3;       % 지도 폭 (cm)  ← 패널 F에 들어갈 실제 폭
lab_fs   = 6;         % 주 이름 글씨 크기 (pt) ← 키우려면 7~8 (지도 폭도 같이 키울 것)
mp_xl1   = -74.3;     % 지도 왼쪽 경계 (경도)
mp_yl    = [-34 5.8]; % 지도 아래/위 경계 (위도)

fan_codes = {'RN','PB','PE','AL','SE'};   % 북동부 해안 5개 주: 위아래 부채꼴 배치
% 지시선 라벨 {코드, dx, dy, 정렬}: dx,dy = 주 중심점에서 이동량 (degree)
free_lab = {'DF', -2.6, -2.6, 'center'; ...
            'ES',  4.5, -0.5, 'left'; ...
            'RJ',  0.0, -3.2, 'center'; ...
            'SC',  5.4, -0.9, 'left'};
in_shift = {'GO', -1.2, 1.8};             % 주 안쪽 라벨 위치 보정 {코드, dx, dy}

Shfile  = shaperead('BR_UF_2023.shp');
mp_n    = numel(Shfile);
mp_code = {Shfile.SIGLA_UF};

% 주별 색(지역), 라벨 기준점(가장 큰 덩어리의 무게중심), 그리기 순서용 면적
mp_reg = zeros(mp_n,1);  mp_cx = zeros(mp_n,1);  mp_cy = zeros(mp_n,1);
mp_xmax = zeros(mp_n,1); mp_area = zeros(mp_n,1);
for i = 1:mp_n
    mp_reg(i) = region_of(strcmp(state_codes, mp_code{i}));
    X = Shfile(i).X(:);  Y = Shfile(i).Y(:);
    mp_area(i) = (max(X)-min(X)) * (max(Y)-min(Y));
    brk = [0; find(isnan(X)); numel(X)+1];
    best = 0;
    for k = 1:numel(brk)-1
        xs = X(brk(k)+1:brk(k+1)-1);  ys = Y(brk(k)+1:brk(k+1)-1);
        if numel(xs) < 3, continue, end
        xa = [xs; xs(1)];  ya = [ys; ys(1)];
        cr = xa(1:end-1).*ya(2:end) - xa(2:end).*ya(1:end-1);
        a  = sum(cr)/2;
        if abs(a) > best
            best = abs(a);
            mp_cx(i)   = sum((xa(1:end-1)+xa(2:end)).*cr)/(6*a);
            mp_cy(i)   = sum((ya(1:end-1)+ya(2:end)).*cr)/(6*a);
            mp_xmax(i) = max(xs);
        end
    end
end

% 해안 라벨 열의 x 위치 → 오른쪽 경계를 글씨 크기에 맞춰 자동 계산
[~, fan_i] = ismember(fan_codes, mp_code);
[~, ord]   = sort(mp_cy(fan_i), 'descend');       % 북 → 남
fan_i      = fan_i(ord);
x_lab      = max(mp_xmax(fan_i)) + 1.3;
Lfrac      = (lab_fs/72*2.54) / map_w_cm;          % 글자 높이 / 지도 폭
mp_xr      = (x_lab + 0.3 - 1.8*Lfrac*mp_xl1) / (1 - 1.8*Lfrac);
lab_h      = Lfrac * (mp_xr - mp_xl1);             % 글자 높이 (degree)
map_h_cm   = map_w_cm * diff(mp_yl) / (mp_xr - mp_xl1);

figure('Units','centimeters','Position',[0 0 map_w_cm map_h_cm],'PaperPositionMode','auto');
axes('Position',[0 0 1 1]); hold on

% 주 그리기 (큰 주부터 → 작은 주(DF 등)가 위에 오도록)
[~, draw_ord] = sort(mp_area, 'descend');
for ii = draw_ord(:)'
    X = Shfile(ii).X(:);  Y = Shfile(ii).Y(:);
    brk = [0; find(isnan(X)); numel(X)+1];
    for k = 1:numel(brk)-1
        xs = X(brk(k)+1:brk(k+1)-1);  ys = Y(brk(k)+1:brk(k+1)-1);
        if numel(xs) < 3, continue, end
        patch(xs, ys, col_region(mp_reg(ii),:), 'EdgeColor','w', 'LineWidth',0.4);
    end
end
axis equal off
xlim([mp_xl1 mp_xr]);  ylim(mp_yl)

txt = @(x,y,c,ha) text(x, y, c, 'FontName','SansSerif', 'FontSize',lab_fs, ...
                       'HorizontalAlignment',ha, 'VerticalAlignment','middle');
dot = @(x,y) plot(x, y, 'ko', 'MarkerSize',1.8, 'MarkerFaceColor','k', 'LineWidth',0.1);

% (1) 주 안쪽 라벨
for i = 1:mp_n
    c = mp_code{i};
    if any(strcmp(c, fan_codes)) || any(strcmp(c, free_lab(:,1))), continue, end
    dxy = [0 0];
    j = find(strcmp(in_shift(:,1), c));
    if ~isempty(j), dxy = [in_shift{j,2}, in_shift{j,3}]; end
    txt(mp_cx(i)+dxy(1), mp_cy(i)+dxy(2), c, 'center');
end

% (2) 북동부 해안 5개 주: 중심 높이 기준으로 위아래로 부채꼴 배치
nf    = numel(fan_i);
y_mid = mean(mp_cy(fan_i));
gap   = 1.15 * lab_h;
for k = 1:nf
    i   = fan_i(k);
    y_l = y_mid + ((nf-1)/2 - (k-1)) * gap;
    plot([mp_cx(i) x_lab-0.3], [mp_cy(i) y_l], 'k-', 'LineWidth',0.4)
    dot(mp_cx(i), mp_cy(i))
    txt(x_lab, y_l, mp_code{i}, 'left');
end

% (3) 개별 지시선 라벨 (DF, ES, RJ, SC)
for j = 1:size(free_lab,1)
    i  = find(strcmp(mp_code, free_lab{j,1}));
    lx = mp_cx(i) + free_lab{j,2};   ly = mp_cy(i) + free_lab{j,3};
    if strcmp(free_lab{j,4}, 'left')
        ex = lx - 0.3;  ey = ly;
    else   % 'center': 지시선이 글자 위를 가리지 않도록 글자 상자 가장자리에서 끝냄
        vx = mp_cx(i) - lx;  vy = mp_cy(i) - ly;  vn = hypot(vx, vy);
        ux = vx/vn;  uy = vy/vn;
        tb = 1.15 * min(0.75*lab_h/max(abs(ux),1e-6), 0.45*lab_h/max(abs(uy),1e-6));
        ex = lx + ux*tb;  ey = ly + uy*tb;
    end
    plot([mp_cx(i) ex], [mp_cy(i) ey], 'k-', 'LineWidth',0.4)
    dot(mp_cx(i), mp_cy(i))
    txt(lx, ly, free_lab{j,1}, free_lab{j,4});
end

% exportgraphics(gcf,'figures/fig_map_region.png','Resolution',600,'BackgroundColor','none');
% exportgraphics(gcf,'figures/fig_map_region.pdf','ContentType','vector','BackgroundColor','none');