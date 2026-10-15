%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Mixed vaccine strategy (children-opt × elderly-opt)
%  지역별/전국 averting ratio 변화, R0 설정(q25/q50/q75)별로 3개 figure
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;
annual_national = 5e6;
tspan = 0:4;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};

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

region_idx = {north_idx, northeast_idx, centwest_idx, southeast_idx, south_idx};
region_lbl = {'N','NE','CW','SE','S'};
col_region = [1.000 0.498 0.055; 0.596 0.306 0.639; 0.090 0.745 0.812; ...
              0.173 0.627 0.173; 0.549 0.341 0.271];

w_grid = 0:0.05:1;   % 21개 점, children 전략 비중 (0=순수 elderly, 1=순수 children)

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
load('results/catalytic_lambda_kappa_draws.mat','kappa_draws','omega_draws','delta_draws');
kappa_adj = kappa_draws;

wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
pop_mat = table2array(readtable('results/pop_age_state_year.csv'));
pop_mat = pop_mat(:, age_reorder_idx);

load('results/SIR_opt_result_dan.mat','doses_opt')   % [S x 3(rc) x 2(tc: 1=children,2=elderly)]

R0_sets  = {R0pred_q25, R0pred_q50, R0pred_q75};
R0_names = {'q25','q50','q75'};

%% ── 2. inputs (최적화 때와 동일 — Shen VE 고정값) ─────────────
inputs = struct();
inputs.agenum    = A5;
inputs.beta0     = median(kappa_adj)';
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
inputs.e_n_d=0.8225;  inputs.e_n_i=0.165;   inputs.e_n_h=0.980; inputs.e_n_m=0.980;
inputs.e_p_d=0.9275;  inputs.e_p_i=0.480;   inputs.e_p_h=0.975; inputs.e_p_m=0.975;
inputs.ew_n_d=0.530;  inputs.ew_n_i=0.100;  inputs.ew_n_h=0.7975; inputs.ew_n_m=0.7975;
inputs.ew_p_d=0.7225; inputs.ew_p_i=0.180;  inputs.ew_p_h=0.880; inputs.ew_p_m=0.880;
inputs.eww_n_d=-0.1675; inputs.eww_n_i=0.015;  inputs.eww_n_h=-0.0825; inputs.eww_n_m=-0.0825;
inputs.eww_p_d=0.180;   inputs.eww_p_i=0.0175; inputs.eww_p_h=0.375;  inputs.eww_p_m=0.375;
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

function [case_a,hosp_a,dead_a,daly_a] = run_sim(si,inputs,vacc_targets, ...
        X0,X1,X2,X3,X4,pop_mat,R0pred,data_lifeexp_age,tspan,T,A5)
    pop_2025     = pop_mat((si-1)*T+T,:)';
    A0 = median(squeeze(X0(:,si,T,:)),1)' .* pop_2025;
    A1 = median(squeeze(X1(:,si,T,:)),1)' .* pop_2025;
    A2 = median(squeeze(X2(:,si,T,:)),1)' .* pop_2025;
    A3 = median(squeeze(X3(:,si,T,:)),1)' .* pop_2025;
    A4 = median(squeeze(X4(:,si,T,:)),1)' .* pop_2025;
    R0si         = R0pred(si,T);
    inputs.beta1 = R0pred(si,T+1:T+5)*52;
    inputs.vacc_targets = vacc_targets;
    y0 = make_y0(A0,A1,A2,A3,A4,R0si,A5);
    [~,y] = ode23(@(t,yy) func_ode_alt_v6(t,yy,inputs), tspan, y0);
    Casen = y(:,A5*48+(1:A5));
    Hosp  = y(:,A5*49+(1:A5));
    Dead  = y(:,A5*50+(1:A5));
    YLD_age = Casen*0.00084 + Hosp*0.0051;
    YLL_age = Dead .* data_lifeexp_age;
    case_a = sum(Casen(end,:));
    hosp_a = sum(Hosp(end,:));
    dead_a = sum(Dead(end,:));
    daly_a = sum(YLD_age(end,:) + YLL_age(end,:));
end

%% ── 4. Baseline(무백신) — R0 설정별, 주별 ─────────────────────
Case_base = zeros(S,3); Hosp_base = zeros(S,3);
Dead_base = zeros(S,3); DALY_base = zeros(S,3);

inp0 = inputs;
children_ages = [2,3];      % 5-14세
elderly_ages  = 13:17;      % 60세+
vacc_targets_zero = struct('idx', {children_ages, elderly_ages}, 'dose', {0, 0});

for rc = 1:3
    R0pred = R0_sets{rc};
    fprintf('Baseline (no vaccine), R0=%s...\n', R0_names{rc})
    for si = 1:S
        [Case_base(si,rc),Hosp_base(si,rc),Dead_base(si,rc),DALY_base(si,rc)] = ...
            run_sim(si,inp0,vacc_targets_zero, op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                    op_seroposall_3,op_seroposall_4, ...
                    pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
    end
end

%% ── 5. 혼합 전략 시뮬레이션 (w grid × R0 3개) ────────────────
% children 연령군엔 w*doses_children(si), elderly 연령군엔 (1-w)*doses_elderly(si)
% 총도즈만 넣고, ODE 내부(func_ode_alt_v6)에서 v5와 동일하게 그 시점 감수성
% pool 비례로 동적 분배함 (정적 인구비 근사 불필요 — 두 그룹이 연령상 겹치지 않음)
nW = numel(w_grid);
Case_mix = zeros(S,3,nW); Hosp_mix = zeros(S,3,nW);
Dead_mix = zeros(S,3,nW); DALY_mix = zeros(S,3,nW);

for rc = 1:3
    R0pred = R0_sets{rc};
    doses_children = doses_opt(:,rc,1);   % [S x 1]
    doses_elderly  = doses_opt(:,rc,2);   % [S x 1]

    for wi = 1:nW
        w = w_grid(wi);
        fprintf('R0=%s, w=%.2f (%d/%d)...\n', R0_names{rc}, w, wi, nW)
        for si = 1:S
            vacc_targets = struct( ...
                'idx',  {children_ages,            elderly_ages}, ...
                'dose', {w*doses_children(si),      (1-w)*doses_elderly(si)});

            [Case_mix(si,rc,wi),Hosp_mix(si,rc,wi),Dead_mix(si,rc,wi),DALY_mix(si,rc,wi)] = ...
                run_sim(si,inputs,vacc_targets, op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                        op_seroposall_3,op_seroposall_4, ...
                        pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
        end
    end
end

save('results/SIR_mix_strategy_dan.mat', ...
     'Case_base','Hosp_base','Dead_base','DALY_base', ...
     'Case_mix','Hosp_mix','Dead_mix','DALY_mix', ...
     'w_grid','doses_opt','state_codes','R0_names')
fprintf('Saved: results/SIR_mix_strategy.mat\n')

%% ── 6. 지역별/전국 averting ratio + Figure (R0 설정별 3개) ────
mkdir('figures')
metric_names = {'Case','Hospitalization','Fatality','DALY'};
metric_base  = {Case_base, Hosp_base, Dead_base, DALY_base};
metric_mix   = {Case_mix, Hosp_mix, Dead_mix, DALY_mix};

for rc = 1:3
    figure('Units','centimeters','Position',[0 0 18 12],'PaperPositionMode','auto');
    tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

    for m = 1:4
        nexttile(m)
        hold on
        base_m = metric_base{m}(:,rc);            % [S x 1]
        mix_m  = squeeze(metric_mix{m}(:,rc,:));  % [S x nW]

        hplot = gobjects(6,1);
        x_axis = 1 - w_grid;   % 0=children전략만, 1=elderly전략만 (방향 반대로)
        for r = 1:5
            idx = region_idx{r};
            base_r  = sum(base_m(idx));
            mix_r   = sum(mix_m(idx,:), 1);        % [1 x nW]
            avert_r = (base_r - mix_r) ./ base_r * 100;
            hplot(r) = plot(x_axis, avert_r, '-', 'Color', col_region(r,:), 'LineWidth', 1.5);
        end
        base_nat  = sum(base_m);
        mix_nat   = sum(mix_m, 1);
        avert_nat = (base_nat - mix_nat) ./ base_nat * 100;
        hplot(6) = plot(x_axis, avert_nat, 'k--', 'LineWidth', 1.5);

        ylabel(sprintf('%s averting ratio (%%)', metric_names{m}))
        xlim([0 1]); xticks([0 0.5 1]); xticklabels({})
        set(gca,'FontName','SansSerif','FontSize',10)
        if m == 1
            legend(hplot, [region_lbl, {'Nationwide'}], 'Location','best','FontSize',8); legend('boxoff')
        end
    end
    sgtitle(sprintf('Mixed strategy (children \\leftrightarrow elderly), R0=%s', R0_names{rc}), ...
        'FontSize',11,'FontName','SansSerif')
    exportgraphics(gcf, sprintf('figures/fig_mix_strategy_%s.jpg', R0_names{rc}), 'Resolution', 600)
    fprintf('Saved: figures/fig_mix_strategy_%s.jpg\n', R0_names{rc})
end