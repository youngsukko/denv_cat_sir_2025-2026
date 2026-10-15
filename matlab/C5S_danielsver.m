%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Optimal vaccine allocation (one case at a time)
%  케이스 번호:
%    1: Children (5-14), R0 q25
%    2: Children (5-14), R0 q50
%    3: Children (5-14), R0 q75
%    4: Elderly  (60+),  R0 q25
%    5: Elderly  (60+),  R0 q50
%    6: Elderly  (60+),  R0 q75
%% ═══════════════════════════════════════════════════════════════
% close all; clear; clc;

% %% ▼▼▼ 여기서 케이스 번호 선택 (1~6) ▼▼▼
case_no = 6;
%% ▲▲▲

%% ── 케이스 해석 ──────────────────────────────────────────────
vacorder_map = [1,1,1, 6,6,6];
r0_map       = [1,2,3, 1,2,3];
target_map   = [1,1,1, 2,2,2];
target_names = {'Children (5-14)', 'Elderly (60+)'};
R0_names     = {'q25','q50','q75'};

vacorder = vacorder_map(case_no);
rc       = r0_map(case_no);
tc       = target_map(case_no);
fprintf('Running: %s | R0=%s | vacorder=%d\n', target_names{tc}, R0_names{rc}, vacorder)

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S = 27; T = 26; A5 = 17;
annual_national = 5e6;
tspan = 0:4;

state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};
ibge_of_state = [12,27,13,16,29,23,53,32,52,21, ...
                 31,50,51,15,25,26,22,41,33,24, ...
                 11,14,43,42,28,35,17];

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
mean_kappa = mean(kappa_draws, 2);
kappa_adj  = kappa_draws;

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
doses_annual_brazil = annual_national * doses_alph / sum(doses_alph);

R0_sets = {R0pred_q25, R0pred_q50, R0pred_q75};
R0pred  = R0_sets{rc};

%% ── 2. inputs ────────────────────────────────────────────────
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

function daly = run_daly(si, inputs, X0,X1,X2,X3,X4, pop_mat, R0pred, ...
        data_lifeexp_age, tspan, T, A5)
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
    Casen = y(end, A5*48+(1:A5));
    Hosp  = y(end, A5*49+(1:A5));
    Dead  = y(end, A5*50+(1:A5));
    YLD = Casen*0.00084 + Hosp*0.0051;
    YLL = Dead .* data_lifeexp_age;
    daly = sum(YLD + YLL);
end

function total_daly = objfun(x, inputs, vacorder, R0pred, annual_national, popbystate, ...
        X0,X1,X2,X3,X4, pop_mat, data_lifeexp_age, tspan, T, A5, S)
    % x → 인구 가중 비율로 doses 변환 (등식/부등식 제약 불필요)
    doses_vec = annual_national * x .* popbystate / sum(x .* popbystate);
    total_daly = 0;
    % if rand < 0.1
        plot(x/mean(x))
        pause(0.01)

    % end
    for si = 1:S
        inp = inputs;
        inp.vacc     = doses_vec(si);
        inp.vacorder = vacorder;
        total_daly = total_daly + run_daly(si, inp, X0,X1,X2,X3,X4, ...
            pop_mat, R0pred, data_lifeexp_age, tspan, T, A5);
    end
end

%% ── 4. popbystate 계산 ────────────────────────────────────────
popbystate = zeros(S,1);
for si = 1:S
    pop_2025 = pop_mat((si-1)*T+T,:)';
    A0tmp = median(squeeze(op_seroposall_0(:,si,T,:)),1)' .* pop_2025;
    A1tmp = median(squeeze(op_seroposall_1(:,si,T,:)),1)' .* pop_2025;
    A2tmp = median(squeeze(op_seroposall_2(:,si,T,:)),1)' .* pop_2025;
    A3tmp = median(squeeze(op_seroposall_3(:,si,T,:)),1)' .* pop_2025;
    A4tmp = median(squeeze(op_seroposall_4(:,si,T,:)),1)' .* pop_2025;
    popbystate(si) = sum([A0tmp;A1tmp;A2tmp;A3tmp;A4tmp],'all');
end

%% ── 5. 초기값: 브라질 배분 → x 역변환 ───────────────────────
% doses = annual_national * x .* pop / sum(x .* pop)
% 역변환: x ∝ doses / pop
% x0 = doses_annual_brazil ./ popbystate;
x0 = popbystate*0 + 1;
x0 = max(x0, 1e-8);   % 0 방지

%% ── 6. 최적화 ────────────────────────────────────────────────
lb = zeros(S,1);   % x > 0

fmincon_opts = optimoptions('fmincon', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations', 300, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-8, ...
    'Algorithm', 'interior-point', ...
    'HessianApproximation', 'bfgs', ...
    'FiniteDifferenceType', 'forward');

obj = @(x) objfun(x, inputs, vacorder, R0pred, annual_national, popbystate, ...
    op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
    op_seroposall_3,op_seroposall_4, ...
    pop_mat, data_lifeexp_age, tspan, T, A5, S);

fprintf('Optimizing...\n')
figure
[x_opt, daly_opt_val] = fmincon(obj, x0, [], [], [], [], lb, [], [], fmincon_opts);

% 최적 doses 복원
doses_opt_si = annual_national * x_opt .* popbystate / sum(x_opt .* popbystate);

daly_brazil_val = obj(x0);

fprintf('\n--- Result: %s | R0=%s ---\n', target_names{tc}, R0_names{rc})
fprintf('  Brazil alloc DALY: %.1f\n', daly_brazil_val)
fprintf('  Optimal      DALY: %.1f\n', daly_opt_val)
fprintf('  Improvement:       %.2f%%\n', (daly_brazil_val-daly_opt_val)/daly_brazil_val*100)

%% ── 7. 기존 결과에 병합해서 저장 ─────────────────────────────
save_path = 'results/SIR_opt_result_dan.mat';
if isfile(save_path)
    prev = load(save_path);
    doses_opt   = prev.doses_opt;
    daly_opt    = prev.daly_opt;
    daly_brazil = prev.daly_brazil;
else
    doses_opt   = zeros(S,3,2);
    daly_opt    = zeros(3,2);
    daly_brazil = zeros(3,2);
end

doses_opt(:, rc, tc) = doses_opt_si;
daly_opt(rc, tc)     = daly_opt_val;
daly_brazil(rc, tc)  = daly_brazil_val;

save(save_path, 'doses_opt','daly_opt','daly_brazil','doses_annual_brazil', ...
    'popbystate','state_codes','target_names','R0_names')
fprintf('Saved (case %d appended): %s\n', case_no, save_path)