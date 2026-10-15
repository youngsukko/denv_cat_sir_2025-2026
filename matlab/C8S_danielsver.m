%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — 대량 백신접종 (5–59세, 인구비례 배분) BoD averted
%
%  Target age : 5–59세 (age index 2:12, vacorder=8 in func_ode_alt_v5)
%  Allocation : 주별 5–59세 인구 비중에 비례하여 전국 도즈를 배분
%  Scenarios  : 전국 연간 1000만 / 2000만 / 3000만 도즈 (5년간 매년 동일 투여)
%  Baseline   : vacc=0, R0 q25/q50/q75
%  Figure     : bar(q50) + errorbar(q25~q75), 시나리오별 비교
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};

age_labels = {'0-','5-','10-','15-','20-','25-','30-','35-', ...
              '40-','45-','50-','55-','60-','65-','70-','75-','80+'};
x_age  = 0:5:80;
tspan  = 0:4;

% 대상 연령: 5–59세 → age index 2:12 (5-9 ... 55-59), func_ode_alt_v5의 vacorder=8과 동일
target_idx = 2:12;
VACORDER_MASS = 8;

% 시나리오: 전국 연간 도즈
scenario_doses = [10e6, 20e6, 30e6];
scenario_names = {'10M','20M','30M'};
nSc = numel(scenario_doses);

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
kappa_adj = kappa_draws;

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

%% ── 2. 공통 inputs ────────────────────────────────────────────
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

%% ── 4. 주별 도즈 배분 (5–59세 인구비례) ───────────────────────
targetpop_state = zeros(S,1);
for si = 1:S
    pop_2025 = pop_mat((si-1)*T+T, :);
    targetpop_state(si) = sum(pop_2025(target_idx));
end
targetpop_share = targetpop_state / sum(targetpop_state);   % [S×1]

% doses_by_state(:,sc) = 시나리오 sc의 전국 도즈를 주별로 배분
doses_by_state = targetpop_share * scenario_doses;          % [S×nSc]

fprintf('\n--- 5-59세 인구비례 배분 (상위 5개 주, 10M 시나리오 기준) ---\n')
[~,ord] = sort(doses_by_state(:,1),'descend');
for k = 1:5
    si = ord(k);
    fprintf('  %s: %.0f doses (%.1f%%)\n', state_codes{si}, ...
        doses_by_state(si,1), 100*targetpop_share(si))
end

%% ── 5. 메인 루프: R0 q25/q50/q75 × baseline/시나리오별 백신 ──
R0_sets  = {R0pred_q25, R0pred_q50, R0pred_q75};
R0_names = {'q25','q50','q75'};

Case_0=zeros(S,A5,3); Hosp_0=zeros(S,A5,3);
Dead_0=zeros(S,A5,3); DALY_0=zeros(S,A5,3);
Case_V=zeros(S,A5,3,nSc); Hosp_V=zeros(S,A5,3,nSc);
Dead_V=zeros(S,A5,3,nSc); DALY_V=zeros(S,A5,3,nSc);

for rc = 1:3
    R0pred = R0_sets{rc};
    fprintf('\nR0 scenario: %s\n', R0_names{rc})

    for si = 1:S
        % Baseline (no vaccine) — 시나리오와 무관하게 1회만 계산
        inp0 = inputs; inp0.vacc=0; inp0.vacorder=VACORDER_MASS;
        [Case_0(si,:,rc),Hosp_0(si,:,rc),Dead_0(si,:,rc),DALY_0(si,:,rc)] = ...
            run_sim(si,inp0, ...
                    op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                    op_seroposall_3,op_seroposall_4, ...
                    pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);

        % Mass vaccination scenarios (5-59세, 인구비례 배분)
        for sc = 1:nSc
            inpV = inputs; inpV.vacc=doses_by_state(si,sc); inpV.vacorder=VACORDER_MASS;
            [Case_V(si,:,rc,sc),Hosp_V(si,:,rc,sc), ...
             Dead_V(si,:,rc,sc),DALY_V(si,:,rc,sc)] = ...
                run_sim(si,inpV, ...
                        op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                        op_seroposall_3,op_seroposall_4, ...
                        pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
        end

        fprintf('  %s [%s] done\n', state_codes{si}, R0_names{rc})
    end
end

mkdir('results')
save('results/SIR_massvax_5to59.mat', ...
     'Case_0','Hosp_0','Dead_0','DALY_0', ...
     'Case_V','Hosp_V','Dead_V','DALY_V', ...
     'doses_by_state','scenario_doses','scenario_names','state_codes')
fprintf('\nSaved: results/SIR_massvax_5to59.mat\n')

%% ── 6. BoD averted: 전국 합산, 시나리오별 [3(R0) × nSc] ──────
nat_case = zeros(3,nSc); nat_hosp = zeros(3,nSc);
nat_dead = zeros(3,nSc); nat_daly = zeros(3,nSc);

for rc = 1:3
    base_case = sum(Case_0(:,:,rc),'all');
    base_hosp = sum(Hosp_0(:,:,rc),'all');
    base_dead = sum(Dead_0(:,:,rc),'all');
    base_daly = sum(DALY_0(:,:,rc),'all');
    for sc = 1:nSc
        nat_case(rc,sc) = base_case - sum(Case_V(:,:,rc,sc),'all');
        nat_hosp(rc,sc) = base_hosp - sum(Hosp_V(:,:,rc,sc),'all');
        nat_dead(rc,sc) = base_dead - sum(Dead_V(:,:,rc,sc),'all');
        nat_daly(rc,sc) = base_daly - sum(DALY_V(:,:,rc,sc),'all');
    end
end
% rc 순서: 1=q25, 2=q50, 3=q75
case_q25=nat_case(1,:); case_q50=nat_case(2,:); case_q75=nat_case(3,:);
hosp_q25=nat_hosp(1,:); hosp_q50=nat_hosp(2,:); hosp_q75=nat_hosp(3,:);
dead_q25=nat_dead(1,:); dead_q50=nat_dead(2,:); dead_q75=nat_dead(3,:);
daly_q25=nat_daly(1,:); daly_q50=nat_daly(2,:); daly_q75=nat_daly(3,:);

fprintf('\n--- BoD averted (q50, 5년 누적, 5-59세 인구비례 배분) ---\n')
for sc = 1:nSc
    vactot_5yr = scenario_doses(sc)*5/1000;
    fprintf('\n[%s doses/yr]\n', scenario_names{sc})
    fprintf('  Cases:            %10.0f  (per 1K vac: %.3f)\n', case_q50(sc), case_q50(sc)/vactot_5yr)
    fprintf('  Hospitalizations: %10.0f  (per 1K vac: %.4f)\n', hosp_q50(sc), hosp_q50(sc)/vactot_5yr)
    fprintf('  Deaths:           %10.4f  (per 1K vac: %.6f)\n', dead_q50(sc), dead_q50(sc)/vactot_5yr)
    fprintf('  DALYs:            %10.4f  (per 1K vac: %.6f)\n', daly_q50(sc), daly_q50(sc)/vactot_5yr)
end

%% ── 7. Figure: 시나리오별 bar(q50) + errorbar(IQR) ───────────
dark_green  = [0.0 0.55 0.0];
light_green = [0.6 1.0  0.6];

q50s = {case_q50, hosp_q50, dead_q50, daly_q50};
q25s = {case_q25, hosp_q25, dead_q25, daly_q25};
q75s = {case_q75, hosp_q75, dead_q75, daly_q75};
ylabels_avt = {'Cases averted','Hospitalizations averted', ...
               'Deaths averted','DALYs averted'};

x_sc = 1:nSc;

mkdir('figures')
figure('Units','centimeters','Position',[0 0 18 12],'PaperPositionMode','auto');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for m = 1:4
    q50 = q50s{m}; q25 = q25s{m}; q75 = q75s{m};
    ymax = max(q75)*1.3;

    nexttile(m)
    bar(x_sc, q50, 0.5,'FaceColor',dark_green,'EdgeColor','none')
    hold on
    errorbar(x_sc, q50, q50-q25, q75-q50, ...
        'LineStyle','none','Color',[0.1 0.1 0.1],'LineWidth',1.0,'CapSize',3)
    xlim([0.4 nSc+0.6]); ylim([0 ymax])
    xticks(x_sc); xticklabels(scenario_names)
    xlabel('National annual doses')
    ylabel(ylabels_avt{m})
    set(gca,'FontName','SansSerif','FontSize',9)
end

sgtitle('BoD averted: Mass vaccination age 5-59 (population-proportional allocation)', ...
    'FontSize',8,'FontName','SansSerif','FontWeight','normal')
exportgraphics(gcf,'figures/fig_BoD_averted_massvax_5to59.jpg','Resolution',600)
fprintf('\nSaved: figures/fig_BoD_averted_massvax_5to59.jpg\n')