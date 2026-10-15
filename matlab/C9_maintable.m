%% ═══════════════════════════════════════════════════════════════
%  Master table: 6개 전략 × case/hosp/dead/daly averted
%  각 셀: "Q50 [Q25,Q75]" 줄바꿈 "per1000vac_Q50 [Q25,Q75]"
%  Daniels/Shen 두 버전 생성 (result_csv/Table_master_Daniels.csv, _Shen.csv)
%
%  전략 6개:
%    1. Brazil distribution      (doses_annual_brazil, vacorder=1, 고정 5M)
%    2. Children-optimized       (doses_opt(:,rc,1), vacorder=1, R0별 배분 다름)
%    3. Elderly-optimized        (doses_opt(:,rc,2), vacorder=6, R0별 배분 다름)
%    4. Mass vaccination 10M     (인구비례, vacorder=8)
%    5. Mass vaccination 20M     (인구비례, vacorder=8)
%    6. Mass vaccination 30M     (인구비례, vacorder=8)
%% ═══════════════════════════════════════════════════════════════
close all; clc;

S = 27; T = 26; A5 = 17;
tspan = 0:4;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};
ibge_of_state = [12,27,13,16,29,23,53,32,52,21, ...
                 31,50,51,15,25,26,22,41,33,24, ...
                 11,14,43,42,28,35,17];

%% ── 1. 데이터 로드 ────────────────────────────────────────────
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

vac_tbl      = readtable('doses_by_state_all_remessas_plus_2025_ibge_sorted.csv','NumHeaderLines',1);
ibge_in_file = vac_tbl.Var1;
doses_sum    = vac_tbl.Var12;
doses_alph   = zeros(S,1);
for si = 1:S
    row = find(ibge_in_file == ibge_of_state(si));
    doses_alph(si) = doses_sum(row);
end
annual_national     = 5e6;
doses_annual_brazil = annual_national * doses_alph / sum(doses_alph);

Sd = load('results/SIR_opt_result_dan.mat', 'doses_opt');   % Daniels 최적화 결과
Ss = load('results/SIR_opt_result.mat', 'doses_opt');        % Shen 최적화 결과

%% ── 2. inputs 기본값 (VE 제외) ───────────────────────────────
inputs_base = struct();
inputs_base.agenum    = A5;
inputs_base.beta0     = median(kappa_adj)';
inputs_base.cageadj   = median(omega_draws)';
inputs_base.cprop_1   = median(delta_draws(:,1));
inputs_base.cprop_2   = median(delta_draws(:,2));
inputs_base.cprop_3   = median(delta_draws(:,3));
inputs_base.cprop_4   = median(delta_draws(:,4));
inputs_base.gamma     = 1/(1/52);
inputs_base.ohm       = 1/1;
inputs_base.i1=3/4; inputs_base.i2=2/4; inputs_base.i3=1/4;
inputs_base.mu        = 1/(5*1);
inputs_base.mu_final  = 1/(8.9*1);
inputs_base.q         = survrate(1:16)';
inputs_base.birthrate = data_birthper1K_since1960(64);
inputs_base.ohm_vac=1/(1); inputs_base.ohm_vac2=1/(4); inputs_base.ohm_vac3=0;

h_ref=[0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
       0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
f_ref=[0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
       0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';
inputs_base.h1=h_ref; inputs_base.h2=h_ref; inputs_base.h3=h_ref; inputs_base.h4=h_ref;
inputs_base.f1=f_ref; inputs_base.f2=f_ref; inputs_base.f3=f_ref; inputs_base.f4=f_ref;

%% ── 3. VE 설정 (Daniels / Shen) ───────────────────────────────
VE_dan.e_n_d=0.8225;  VE_dan.e_n_i=0.165;   VE_dan.e_n_h=0.980; VE_dan.e_n_m=0.980;
VE_dan.e_p_d=0.9275;  VE_dan.e_p_i=0.480;   VE_dan.e_p_h=0.975; VE_dan.e_p_m=0.975;
VE_dan.ew_n_d=0.530;  VE_dan.ew_n_i=0.100;  VE_dan.ew_n_h=0.7975; VE_dan.ew_n_m=0.7975;
VE_dan.ew_p_d=0.7225; VE_dan.ew_p_i=0.180;  VE_dan.ew_p_h=0.880; VE_dan.ew_p_m=0.880;
VE_dan.eww_n_d=-0.1675; VE_dan.eww_n_i=0.015;  VE_dan.eww_n_h=-0.0825; VE_dan.eww_n_m=-0.0825;
VE_dan.eww_p_d=0.180;   VE_dan.eww_p_i=0.0175; VE_dan.eww_p_h=0.375;  VE_dan.eww_p_m=0.375;

VE_shen.e_n_d=0.562;  VE_shen.e_n_i=0.281;  VE_shen.e_n_h=0.844; VE_shen.e_n_m=0.844;
VE_shen.e_p_d=0.720;  VE_shen.e_p_i=0.360;  VE_shen.e_p_h=0.844; VE_shen.e_p_m=0.844;
VE_shen.ew_n_d=0.405; VE_shen.ew_n_i=0.2025; VE_shen.ew_n_h=0.844; VE_shen.ew_n_m=0.844;
VE_shen.ew_p_d=0.531; VE_shen.ew_p_i=0.2655; VE_shen.ew_p_h=0.844; VE_shen.ew_p_m=0.844;
VE_shen.eww_n_d=0.233; VE_shen.eww_n_i=0.1165; VE_shen.eww_n_h=0.844; VE_shen.eww_n_m=0.844;
VE_shen.eww_p_d=0.319; VE_shen.eww_p_i=0.1595; VE_shen.eww_p_h=0.844; VE_shen.eww_p_m=0.844;

function inp = apply_ve(inp, VE)
    fn = fieldnames(VE);
    for i = 1:numel(fn); inp.(fn{i}) = VE.(fn{i}); end
end

function s = comma_num(x)
    % 정수로 반올림 + 천단위 쉼표
    s = sprintf('%.0f', round(x));
    neg = startsWith(s,'-');
    if neg, s = s(2:end); end
    n = length(s);
    if n > 3
        first = mod(n,3); if first==0, first=3; end
        parts = {s(1:first)};
        for i = first+1:3:n
            parts{end+1} = s(i:i+2); %#ok<AGROW>
        end
        s = strjoin(parts, ',');
    end
    if neg, s = ['-' s]; end
end

%% ── 4. 헬퍼 함수 ─────────────────────────────────────────────
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
    case_a = sum(Casen(end,:));
    hosp_a = sum(Hosp(end,:));
    dead_a = sum(Dead(end,:));
    daly_a = sum(YLD_age(end,:) + YLL_age(end,:));
end

%% ── 5. 인구 배분 기준 ──────────────────────────────────────────
% popbystate: 전체 인구(모든 연령/혈청상태) — elderly/children-optimized 무관, 참고용
popbystate = zeros(S,1);
% mass vaccination(vacorder=8, 5-59세)은 반드시 5-59세 인구비로 배분해야 함
% (document22와 동일 방식 — 전체 인구비를 쓰면 배분이 달라짐)
mass_target_idx = 2:12;   % 5-59세
mass_targetpop_state = zeros(S,1);
for si = 1:S
    pop_2025 = pop_mat((si-1)*T+T, :)';
    A0 = median(squeeze(op_seroposall_0(:,si,T,:)),1)' .* pop_2025;
    A1 = median(squeeze(op_seroposall_1(:,si,T,:)),1)' .* pop_2025;
    A2 = median(squeeze(op_seroposall_2(:,si,T,:)),1)' .* pop_2025;
    A3 = median(squeeze(op_seroposall_3(:,si,T,:)),1)' .* pop_2025;
    A4 = median(squeeze(op_seroposall_4(:,si,T,:)),1)' .* pop_2025;
    popbystate(si) = sum([A0;A1;A2;A3;A4],'all');
    mass_targetpop_state(si) = sum(pop_2025(mass_target_idx));
end
mass_pop_frac = mass_targetpop_state / sum(mass_targetpop_state);   % 5-59세 인구비 (mass vax 배분 기준)

%% ── 6. baseline (무백신) — R0 3개, VE 무관 ────────────────────
R0_sets = {R0pred_q25, R0pred_q50, R0pred_q75};
baseline_cache = 'results/master_table_baseline.mat';
if isfile(baseline_cache)
    fprintf('캐시 발견 — baseline 재실행 없이 로드: %s\n', baseline_cache)
    ld = load(baseline_cache);
    Case_0 = ld.Case_0; Hosp_0 = ld.Hosp_0; Dead_0 = ld.Dead_0; DALY_0 = ld.DALY_0;
else
    Case_0 = zeros(S,3); Hosp_0 = zeros(S,3); Dead_0 = zeros(S,3); DALY_0 = zeros(S,3);
    inp0 = inputs_base; inp0.vacc = 0; inp0.vacorder = 1;
    inp0 = apply_ve(inp0, VE_shen);  % VE는 vacc=0이라 실제 안 쓰임, placeholder
    for rc = 1:3
        R0pred = R0_sets{rc};
        fprintf('Baseline R0 scenario %d/3...\n', rc)
        for si = 1:S
            [Case_0(si,rc),Hosp_0(si,rc),Dead_0(si,rc),DALY_0(si,rc)] = ...
                run_sim(si,inp0, op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                        op_seroposall_3,op_seroposall_4, ...
                        pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
        end
    end
    mkdir('results')
    save(baseline_cache, 'Case_0','Hosp_0','Dead_0','DALY_0')
    fprintf('저장(캐시): %s\n', baseline_cache)
end

%% ── 7. 전략 정의 ──────────────────────────────────────────────
% doses_fn(rc) → [S x 1] 그 R0시나리오에서 각 주에 줄 도즈. vacorder는 고정.
strategy_names = {'Brazil distribution','Children-optimized','Elderly-optimized', ...
                   'Mass vaccination 10M','Mass vaccination 20M','Mass vaccination 30M'};
strategy_vacorder = [1, 1, 6, 8, 8, 8];

strategy_doses_dan  = { ...
    @(rc) doses_annual_brazil, ...
    @(rc) Sd.doses_opt(:,rc,1), ...
    @(rc) Sd.doses_opt(:,rc,2), ...
    @(rc) 10e6*mass_pop_frac, @(rc) 20e6*mass_pop_frac, @(rc) 30e6*mass_pop_frac};
strategy_doses_shen = { ...
    @(rc) doses_annual_brazil, ...
    @(rc) Ss.doses_opt(:,rc,1), ...
    @(rc) Ss.doses_opt(:,rc,2), ...
    @(rc) 10e6*mass_pop_frac, @(rc) 20e6*mass_pop_frac, @(rc) 30e6*mass_pop_frac};
strategy_total_doses = [annual_national, annual_national, annual_national, 10e6, 20e6, 30e6];

%% ── 8. 메인 루프: VE(Daniels/Shen) × 전략(6) × R0(3) ──────────

VE_sets_all   = {VE_dan, VE_shen};
strategy_sets = {strategy_doses_dan, strategy_doses_shen};
ve_file_tags  = {'Daniels','Shen'};

for vi = 1:2
    VE = VE_sets_all{vi};
    strategy_doses = strategy_sets{vi};

    cache_file = sprintf('results/master_table_raw_%s.mat', ve_file_tags{vi});
    if isfile(cache_file)
        fprintf('캐시 발견 — 시뮬레이션 재실행 없이 로드: %s\n', cache_file)
        ld = load(cache_file);
        Case_avt = ld.Case_avt; Hosp_avt = ld.Hosp_avt;
        Dead_avt = ld.Dead_avt; DALY_avt = ld.DALY_avt;
    else
        Case_avt = zeros(6,3); Hosp_avt = zeros(6,3);
        Dead_avt = zeros(6,3); DALY_avt = zeros(6,3);

        for st = 1:6
            vacorder = strategy_vacorder(st);
            for rc = 1:3
                R0pred = R0_sets{rc};
                doses_state = strategy_doses{st}(rc);
                fprintf('[%s] %s, R0 scenario %d/3...\n', ve_file_tags{vi}, strategy_names{st}, rc)

                case_v=0; hosp_v=0; dead_v=0; daly_v=0;
                for si = 1:S
                    inpV = apply_ve(inputs_base, VE);
                    inpV.vacc     = doses_state(si);
                    inpV.vacorder = vacorder;
                    [ca,ho,de,da] = run_sim(si,inpV, ...
                        op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                        op_seroposall_3,op_seroposall_4, ...
                        pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
                    case_v = case_v+ca; hosp_v = hosp_v+ho; dead_v = dead_v+de; daly_v = daly_v+da;
                end

                Case_avt(st,rc) = sum(Case_0(:,rc)) - case_v;
                Hosp_avt(st,rc) = sum(Hosp_0(:,rc)) - hosp_v;
                Dead_avt(st,rc) = sum(Dead_0(:,rc)) - dead_v;
                DALY_avt(st,rc) = sum(DALY_0(:,rc)) - daly_v;
            end
        end

        mkdir('results')
        save(cache_file, 'Case_avt','Hosp_avt','Dead_avt','DALY_avt', ...
             'strategy_names','strategy_total_doses')
        fprintf('저장(캐시): %s\n', cache_file)
    end

    %% ── 표 조립: median[Q25,Q75] 줄바꿈 per1000vac[Q25,Q75] ──
    vactot_row = strategy_total_doses * 5 / 1000;   % [1x6], 전략별 5년 누적 도즈/1000

    metric_avt = {Case_avt, Hosp_avt, Dead_avt, DALY_avt};
    col_names  = {'case_averted','hosp_averted','dead_averted','daly_averted'};

    tbl = cell(6,4);
    for mi = 1:4
        avt = metric_avt{mi};   % [6 x 3] (rc=1,2,3 → q25,q50,q75)
        for st = 1:6
            m50 = avt(st,2); m25 = avt(st,1); m75 = avt(st,3);
            m_lo = min(m25,m75); m_hi = max(m25,m75);

            p50 = m50/vactot_row(st); p25 = m25/vactot_row(st); p75 = m75/vactot_row(st);
            p_lo = min(p25,p75); p_hi = max(p25,p75);

            raw_line1 = comma_num(m50);
            raw_line2 = sprintf('[%s, %s]', comma_num(m_lo), comma_num(m_hi));
            per_line1 = sprintf('(%.2f', p50);
            per_line2 = sprintf('[%.2f, %.2f])', p_lo, p_hi);
            tbl{st,mi} = sprintf('%s\n%s\n%s\n%s', raw_line1, raw_line2, per_line1, per_line2);
        end
    end
    T_master = cell2table(tbl, 'VariableNames', col_names, 'RowNames', strategy_names);
    mkdir('result_csv')
    writetable(T_master, sprintf('result_csv/Table_master_%s.csv', ve_file_tags{vi}), 'WriteRowNames', true)
    fprintf('Saved: result_csv/Table_master_%s.csv\n', ve_file_tags{vi})
end