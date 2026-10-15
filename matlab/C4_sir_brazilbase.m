%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Vaccination BoD averted (R0 IQR 반영)
%
%  Baseline : vacc=0, R0 q25/q50/q75
%  Scenario : 브라질 배분 비율, 전국 연간 500만 명, R0 q25/q50/q75
%  Figure   : bar(q50) + errorbar(q25~q75)
%% ═══════════════════════════════════════════════════════════════
% close all; clear; clc;

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};

ibge_of_state = [12,27,13,16,29,23,53,32,52,21, ...
                 31,50,51,15,25,26,22,41,33,24, ...
                 11,14,43,42,28,35,17];

age_labels = {'0-','5-','10-','15-','20-','25-','30-','35-', ...
              '40-','45-','50-','55-','60-','65-','70-','75-','80+'};
x_age  = 0:5:80;
tspan  = 0:4;
annual_national = 5e6;

%% ── 1. 데이터 로드 ────────────────────────────────────────────
fprintf('Loading data...\n')

load('results/catalytic_ppc_draws.mat', ...
    'op_seroposall_0','op_seroposall_1','op_seroposall_2', ...
    'op_seroposall_3','op_seroposall_4');

load('results/R0predict_new.mat','R0pred_q25','R0pred_q50','R0pred_q75')

load('data_lifeexp.mat','data_lifeexp')
data_lifeexp_age = data_lifeexp(1:5:end)';

load('data_survrate.mat'); load('data_birthrate.mat')
data_birthper1K_since1960(end:end+50) = data_birthper1K_since1960(end);
survrate = survrate(1:5:end)';

load('results/catalytic_lambda_kappa_draws.mat', ...
    'kappa_draws','omega_draws','delta_draws');
% kappa_age 평균 정규화
mean_kappa = mean(kappa_draws, 2);
kappa_adj  = kappa_draws;

% pop 연령 컬럼 순서 버그 보정
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
assert(all(age_reorder_idx > 0), '연령 컬럼 매핑 실패');

pop_mat = table2array(readtable('results/pop_age_state_year.csv'));
pop_mat = pop_mat(:, age_reorder_idx);

vac_tbl      = readtable('doses_by_state_all_remessas_plus_2025_ibge_sorted.csv','NumHeaderLines',1);
ibge_in_file = vac_tbl.Var1;
doses_sum    = vac_tbl.Var12;
doses_alph   = zeros(S,1);
for si = 1:S
    row = find(ibge_in_file == ibge_of_state(si));
    doses_alph(si) = doses_sum(row);
end
doses_annual = annual_national * doses_alph / sum(doses_alph);

%% ── 2. 공통 inputs ────────────────────────────────────────────
inputs = struct();
inputs.agenum    = A5;
inputs.beta0     = median(kappa_adj)';         % [A5×1] 정규화됨, 평균=1
inputs.cageadj   = median(omega_draws)';
inputs.cprop_1   = median(delta_draws(:,1));
inputs.cprop_2   = median(delta_draws(:,2));
inputs.cprop_3   = median(delta_draws(:,3));
inputs.cprop_4   = median(delta_draws(:,4));
inputs.gamma     = 1/(1/52);
inputs.ohm       = 1/1;
inputs.i1=3/4; inputs.i2=2/4; inputs.i3=1/4;
inputs.mu        = 1/(5*1);
inputs.mu_final  = 1/(8.9*1);
inputs.q         = survrate(1:16)';
inputs.birthrate = data_birthper1K_since1960(64);

inputs.e_n_d=0.562;  inputs.e_n_i=0.281;  inputs.e_n_h=0.844; inputs.e_n_m=0.844;
inputs.e_p_d=0.720;  inputs.e_p_i=0.360;  inputs.e_p_h=0.844; inputs.e_p_m=0.844;
inputs.ew_n_d=0.405; inputs.ew_n_i=0.2025; inputs.ew_n_h=0.844; inputs.ew_n_m=0.844;
inputs.ew_p_d=0.531; inputs.ew_p_i=0.2655; inputs.ew_p_h=0.844; inputs.ew_p_m=0.844;
inputs.eww_n_d=0.233; inputs.eww_n_i=0.1165; inputs.eww_n_h=0.844; inputs.eww_n_m=0.844;
inputs.eww_p_d=0.319; inputs.eww_p_i=0.1595; inputs.eww_p_h=0.844; inputs.eww_p_m=0.844;
inputs.ohm_vac=1/(1); inputs.ohm_vac2=1/(4); inputs.ohm_vac3=0;

h_ref=[0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
       0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
f_ref=[0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
       0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';
inputs.h1=h_ref; inputs.h2=h_ref; inputs.h3=h_ref; inputs.h4=h_ref;
inputs.f1=f_ref; inputs.f2=f_ref; inputs.f3=f_ref; inputs.f4=f_ref;

%% ── 3. 헬퍼 함수 ─────────────────────────────────────────────
function y0 = make_y0(A0,A1,A2,A3,A4,R0,A5)
    mmu = 1/76.4;
    sc  = mmu/(mmu+52)*(1-1/R0);
    y0  = zeros(51*A5,1);
    y0(A5*0+(1:A5))  = A0*(1-sc);
    y0(A5*3+(1:A5))  = A1*(1-sc);
    y0(A5*6+(1:A5))  = A2*(1-sc);
    y0(A5*9+(1:A5))  = A3*(1-sc);
    y0(A5*11+(1:A5)) = A4;
    y0(A5*1+(1:A5))  = A0*sc;
    y0(A5*4+(1:A5))  = A1*sc;
    y0(A5*7+(1:A5))  = A2*sc;
    y0(A5*10+(1:A5)) = A3*sc;
end

function [case_a,hosp_a,dead_a,daly_a] = run_sim(si,inputs, ...
        X0,X1,X2,X3,X4,pop_mat,R0pred,data_lifeexp_age,tspan,T,A5)
    pop_2025     = pop_mat((si-1)*T+T,:)';
    A0 = median(squeeze(X0(:,si,T,:)),1)' .* pop_2025;
    A1 = median(squeeze(X1(:,si,T,:)),1)' .* pop_2025;
    A2 = median(squeeze(X2(:,si,T,:)),1)' .* pop_2025;
    A3 = median(squeeze(X3(:,si,T,:)),1)' .* pop_2025;
    A4 = median(squeeze(X4(:,si,T,:)),1)' .* pop_2025;
    R0si         = R0pred(si,T);
    inputs.beta1 = R0pred(si,T+1:T+5)*52;
    y0 = make_y0(A0,A1,A2,A3,A4,R0si,A5);
    [~,y] = ode23(@(t,yy) func_ode_alt_v5(t,yy,inputs), tspan, y0);
    Casen = y(:,A5*48+(1:A5));
    Hosp  = y(:,A5*49+(1:A5));
    Dead  = y(:,A5*50+(1:A5));
    YLD_age = Casen*0.00084 + Hosp*0.0051;
    YLL_age = Dead .* data_lifeexp_age;
    case_a = Casen(end,:);
    hosp_a = Hosp(end,:);
    dead_a = Dead(end,:);
    daly_a = YLD_age(end,:) + YLL_age(end,:);
end

%% ── 4. 메인 루프: R0 q25/q50/q75 × baseline/vaccine ──────────
R0_sets  = {R0pred_q25, R0pred_q50, R0pred_q75};
R0_names = {'q25','q50','q75'};

Case_0=zeros(S,A5,3); Hosp_0=zeros(S,A5,3);
Dead_0=zeros(S,A5,3); DALY_0=zeros(S,A5,3);
Case_V=zeros(S,A5,3); Hosp_V=zeros(S,A5,3);
Dead_V=zeros(S,A5,3); DALY_V=zeros(S,A5,3);

for rc = 1:3
    R0pred = R0_sets{rc};
    fprintf('\nR0 scenario: %s\n', R0_names{rc})

    for si = 1:S
        % Baseline (no vaccine)
        inp0 = inputs; inp0.vacc=0; inp0.vacorder=1;
        [Case_0(si,:,rc),Hosp_0(si,:,rc),Dead_0(si,:,rc),DALY_0(si,:,rc)] = ...
            run_sim(si,inp0, ...
                    op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                    op_seroposall_3,op_seroposall_4, ...
                    pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);

        % Vaccine (브라질 실제 배분 비율, 연간 500만 명)
        inpV = inputs; inpV.vacc=doses_annual(si); inpV.vacorder=1;
        [Case_V(si,:,rc),Hosp_V(si,:,rc),Dead_V(si,:,rc),DALY_V(si,:,rc)] = ...
            run_sim(si,inpV, ...
                    op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                    op_seroposall_3,op_seroposall_4, ...
                    pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);

        fprintf('  %s [%s] done\n', state_codes{si}, R0_names{rc})
    end
end

save('results/SIR_vaccine_brazil_alloc.mat', ...
     'Case_0','Hosp_0','Dead_0','DALY_0', ...
     'Case_V','Hosp_V','Dead_V','DALY_V', ...
     'doses_annual','state_codes')
fprintf('Saved: results/SIR_vaccine_brazil_alloc.mat\n')

%% ── 5. BoD averted: 전국 합산 [1×A5×3] ──────────────────────
avt_case = sum(Case_0-Case_V,1);   % [1×A5×3]
avt_hosp = sum(Hosp_0-Hosp_V,1);
avt_dead = sum(Dead_0-Dead_V,1);
avt_daly = sum(DALY_0-DALY_V,1);

ac_q25=avt_case(:,:,1); ac_q50=avt_case(:,:,2); ac_q75=avt_case(:,:,3);
ah_q25=avt_hosp(:,:,1); ah_q50=avt_hosp(:,:,2); ah_q75=avt_hosp(:,:,3);
ad_q25=avt_dead(:,:,1); ad_q50=avt_dead(:,:,2); ad_q75=avt_dead(:,:,3);
ay_q25=avt_daly(:,:,1); ay_q50=avt_daly(:,:,2); ay_q75=avt_daly(:,:,3);

total_doses_5yr = annual_national*5;
vactot = total_doses_5yr/1000;

fprintf('\n--- BoD averted (q50, 5M/yr × 5yr, Brazil allocation) ---\n')
fprintf('Cases:            %9.0f  (per 1K vac: %.3f)\n', sum(ac_q50), sum(ac_q50)/vactot)
fprintf('Hospitalizations: %9.0f  (per 1K vac: %.4f)\n', sum(ah_q50), sum(ah_q50)/vactot)
fprintf('Deaths:           %9.4f  (per 1K vac: %.6f)\n', sum(ad_q50), sum(ad_q50)/vactot)
fprintf('DALYs:            %9.4f  (per 1K vac: %.6f)\n', sum(ay_q50), sum(ay_q50)/vactot)

%% ── 6. Figure: bar(q50) + errorbar(IQR) ─────────────────────
% dark_green  = [0.0 0.55 0.0];
% light_green = [0.6 1.0  0.6];
% shade_x = [2.5 12.5 12.5 2.5];
% 
% q50s = {ac_q50, ah_q50, ad_q50, ay_q50};
% q25s = {ac_q25, ah_q25, ad_q25, ay_q25};
% q75s = {ac_q75, ah_q75, ad_q75, ay_q75};
% ylabels_avt = {'Cases averted','Hospitalizations averted', ...
%                'Deaths averted','DALYs averted'};
% 
% mkdir('figures')
% figure('Units','centimeters','Position',[0 0 18 12],'PaperPositionMode','auto');
% tiledlayout(2,2,'TileSpacing','compact','Padding','compact')
% 
% for m = 1:4
%     q50 = q50s{m}; q25 = q25s{m}; q75 = q75s{m};
%     ymax = max(q75)*1.3;
% 
%     nexttile(m)
%     fill(shade_x,[0 0 ymax ymax],light_green,'EdgeColor','none','FaceAlpha',0.5)
%     hold on
%     bar(x_age, q50, 0.8,'FaceColor',dark_green,'EdgeColor','none')
%     errorbar(x_age, q50, q50-q25, q75-q50, ...
%         'LineStyle','none','Color',[0.1 0.1 0.1],'LineWidth',1.0,'CapSize',3)
%     xlim([-2.5 82.5]); ylim([0 ymax])
%     xticks(x_age); xticklabels(age_labels)
%     ylabel(ylabels_avt{m})
%     set(gca,'FontName','SansSerif','FontSize',11)
% end
% 
% % sgtitle('BoD averted: Brazil allocation (5M doses/yr × 5yr, vacorder=1)', ...
% %     'FontSize',10,'FontName','SansSerif')
% exportgraphics(gcf,'figures/fig_BoD_averted_brazil_alloc.jpg','Resolution',600)
% fprintf('Saved: figures/fig_BoD_averted_brazil_alloc.jpg\n')

%% ── 6. Figure: bar(q50, per 1000 vaccinated) + errorbar(IQR) + Total 텍스트 ──
dark_green  = [0.0 0.55 0.0];
light_green = [0.6 1.0  0.6];
shade_x = [2.5 12.5 12.5 2.5];

q50s = {ac_q50, ah_q50, ad_q50, ay_q50};
q25s = {ac_q25, ah_q25, ad_q25, ay_q25};
q75s = {ac_q75, ah_q75, ad_q75, ay_q75};
% ylabels_avt = {'Cases averted (per 1000 vaccinated)','Hospitalizations averted (per 1000 vaccinated)', ...
%                'Deaths averted (per 1000 vaccinated)','DALYs averted (per 1000 vaccinated)'};
ylabels_avt = {{'Cases averted',            '(per 1000 vaccinated)'}, ...
               {'Hospitalizations averted', '(per 1000 vaccinated)'}, ...
               {'Deaths averted',           '(per 1000 vaccinated)'}, ...
               {'DALYs averted',            '(per 1000 vaccinated)'}};

mkdir('figures')
figure('Units','centimeters','Position',[0 0 18 12],'PaperPositionMode','auto');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for m = 1:4
    % per 1000 vaccinated로 변환
    q50 = q50s{m}/vactot; q25 = q25s{m}/vactot; q75 = q75s{m}/vactot;
    ymax = max(q75)*1.3;

    nexttile(m)
    fill(shade_x,[0 0 ymax ymax],light_green,'EdgeColor','none','FaceAlpha',0.5)
    hold on
    bar(x_age, q50, 0.8,'FaceColor',dark_green,'EdgeColor','none')
    errorbar(x_age, q50, q50-q25, q75-q50, ...
        'LineStyle','none','Color',[0.1 0.1 0.1],'LineWidth',1.0,'CapSize',3)
    xlim([-2.5 82.5]); ylim([0 ymax])
    xticks(x_age); xticklabels(age_labels)
    set(gca,'FontName','SansSerif','FontSize',11)
    ylabel(ylabels_avt{m}, 'FontSize', 9)

    % 전 연령 합계 (per 1000 vaccinated) — 우측 상단 텍스트
    total_50 = sum(q50); total_25 = sum(q25); total_75 = sum(q75);
    total_lo = min(total_25, total_75); total_hi = max(total_25, total_75);
    % text(0.995, 0.95, sprintf('Total: %.2f [%.2f, %.2f] per 1000 vaccinated', ...
        % total_50, total_lo, total_hi), ...
    % text(0.995, 0.95, sprintf('Total: %.2f [%.2f, %.2f]', ...
    %     total_50, total_lo, total_hi), ...
    if m == 3, fmt_tot = 'Total: %.4f [%.4f, %.4f]';   % Deaths: 소수점 넷째자리
    else,      fmt_tot = 'Total: %.2f [%.2f, %.2f]';
    end
    text(0.995, 0.95, sprintf(fmt_tot, total_50, total_lo, total_hi), ...
        'Units','normalized', 'HorizontalAlignment','right', 'VerticalAlignment','top', ...
        'FontName','SansSerif','FontSize',9)
end

% sgtitle('BoD averted: Brazil allocation (5M doses/yr × 5yr, vacorder=1)', ...
%     'FontSize',10,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_BoD_averted_brazil_alloc.jpg','Resolution',600)
fprintf('Saved: figures/fig_BoD_averted_brazil_alloc.jpg\n')