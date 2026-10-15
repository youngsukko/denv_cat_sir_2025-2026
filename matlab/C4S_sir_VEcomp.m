%% VE 설정값 추출 (Daniels / Shen / Average)
close all; clear; clc;

dan_sn_symp = [84  55  41  29   16     2  -10  -28  -37 -45 -51;...
              100 100 100 100  100   100  100  100   97  96  93;...
               76  33  18   2  -10   -22  -34  -45  -49 -59 -65;...
               69  24  10  -2  -16   -29  -40  -47  -55 -59 -67];
dan_sp_symp = [95  74  60  50   42    35   29   23   18  14  12;...
              100  91  81  75   68    60   53   45   38  32  25;...
               82  48  32  25   19    15   12   10    7   6   5;...
               94  76  65  57   50    43   36   30   24  20  15];
dan_sn_hosp = [100  87  83  70  62   49   40   23    6  -10  -23;...
              100 100 100 100 100  100  100  100  100  100  100;...
               92  62  49  32  15   -6  -27  -45  -70  -91 -100;...
              100  70  61  49  36   27    2  -10  -19  -32  -40];
dan_sp_hosp = [ 97  91  83  78  72   64   56   49   40   33   29;...
               100  95  94  90  87   84   79   74   66   58   49;...
                95  75  60  51  43   35   30   24   19   16   13;...
                98  91  84  80  75   68   62   56   49   43   36];
dan_sn_inf = [0 0 0 0 0 0 0 0 0 0 0;...
              66 40 32 26 20 17 13 10 8 6 5;...
              0 0 0 0 0 0 0 0 0 0 0;...
              0 0 0 0 0 0 0 0 0 0 0];
dan_sp_inf  = [48  15   8   7   5    3    2    2    2    1    1;...
               71  32  21  16  13   11    7    7    5    5    2;...
               20   5   3   1   1    0    0    0    0    0    0;...
               53  20  12  10   7    5    5    5    2    1    1];

shen_sn_symp = [56.2 40.5 34.8 31.4 29.2 27.5 26.1 25.0 24.1 23.3 22.6];
shen_sp_symp = [72.0 53.1 46.2 42.1 39.2 37.1 35.5 34.1 32.9 31.9 31.1];
shen_sn_inf  = shen_sn_symp / 2;
shen_sp_inf  = shen_sp_symp / 2;
shen_hosp    = 84.4;

idx1 = 1; idx2 = 2; idx3 = 10;

% ── Daniels ──────────────────────────────────────────────────
D.e_n_d   = mean(dan_sn_symp(:,idx1))/100;
D.e_n_i   = mean(dan_sn_inf(:,idx1))/100;
D.e_n_h   = mean(dan_sn_hosp(:,idx1))/100;
D.e_p_d   = mean(dan_sp_symp(:,idx1))/100;
D.e_p_i   = mean(dan_sp_inf(:,idx1))/100;
D.e_p_h   = mean(dan_sp_hosp(:,idx1))/100;

D.ew_n_d  = mean(dan_sn_symp(:,idx2))/100;
D.ew_n_i  = mean(dan_sn_inf(:,idx2))/100;
D.ew_n_h  = mean(dan_sn_hosp(:,idx2))/100;
D.ew_p_d  = mean(dan_sp_symp(:,idx2))/100;
D.ew_p_i  = mean(dan_sp_inf(:,idx2))/100;
D.ew_p_h  = mean(dan_sp_hosp(:,idx2))/100;

D.eww_n_d = mean(dan_sn_symp(:,idx3))/100;
D.eww_n_i = mean(dan_sn_inf(:,idx3))/100;
D.eww_n_h = mean(dan_sn_hosp(:,idx3))/100;
D.eww_p_d = mean(dan_sp_symp(:,idx3))/100;
D.eww_p_i = mean(dan_sp_inf(:,idx3))/100;
D.eww_p_h = mean(dan_sp_hosp(:,idx3))/100;

% ── Shen ─────────────────────────────────────────────────────
H.e_n_d   = shen_sn_symp(idx1)/100;
H.e_n_i   = shen_sn_inf(idx1)/100;
H.e_n_h   = shen_hosp/100;
H.e_p_d   = shen_sp_symp(idx1)/100;
H.e_p_i   = shen_sp_inf(idx1)/100;
H.e_p_h   = shen_hosp/100;

H.ew_n_d  = shen_sn_symp(idx2)/100;
H.ew_n_i  = shen_sn_inf(idx2)/100;
H.ew_n_h  = shen_hosp/100;
H.ew_p_d  = shen_sp_symp(idx2)/100;
H.ew_p_i  = shen_sp_inf(idx2)/100;
H.ew_p_h  = shen_hosp/100;

H.eww_n_d = shen_sn_symp(idx3)/100;
H.eww_n_i = shen_sn_inf(idx3)/100;
H.eww_n_h = shen_hosp/100;
H.eww_p_d = shen_sp_symp(idx3)/100;
H.eww_p_i = shen_sp_inf(idx3)/100;
H.eww_p_h = shen_hosp/100;

% ── Average ───────────────────────────────────────────────────
fields = fieldnames(D);
A = struct();
for fi = 1:numel(fields)
    A.(fields{fi}) = (D.(fields{fi}) + H.(fields{fi})) / 2;
end

VE_Daniels = D;
VE_Shen    = H;
VE_Average = A;

save('results/VE_settings.mat', 'VE_Daniels','VE_Shen','VE_Average')
fprintf('Saved: results/VE_settings.mat\n')

%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Vaccination BoD averted (VE sensitivity: Daniels / Shen / Average)
%% ═══════════════════════════════════════════════════════════════

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

%% ── 1. VE 계산 (doc12 코드로부터) ────────────────────────────
dan_sn_symp = [84  55  41  29   16     2  -10  -28  -37 -45 -51;...
              100 100 100 100  100   100  100  100   97  96  93;...
               76  33  18   2  -10   -22  -34  -45  -49 -59 -65;...
               69  24  10  -2  -16   -29  -40  -47  -55 -59 -67];
dan_sp_symp = [95  74  60  50   42    35   29   23   18  14  12;...
              100  91  81  75   68    60   53   45   38  32  25;...
               82  48  32  25   19    15   12   10    7   6   5;...
               94  76  65  57   50    43   36   30   24  20  15];
dan_sn_hosp = [100  87  83  70  62   49   40   23    6  -10  -23;...
              100 100 100 100 100  100  100  100  100  100  100;...
               92  62  49  32  15   -6  -27  -45  -70  -91 -100;...
              100  70  61  49  36   27    2  -10  -19  -32  -40];
dan_sp_hosp = [ 97  91  83  78  72   64   56   49   40   33   29;...
               100  95  94  90  87   84   79   74   66   58   49;...
                95  75  60  51  43   35   30   24   19   16   13;...
                98  91  84  80  75   68   62   56   49   43   36];
dan_sn_inf  = zeros(4,11);   % Daniels seronegative infection VE = 0
dan_sp_inf  = [48  15   8   7   5    3    2    2    2    1    1;...
               71  32  21  16  13   11    7    7    5    5    2;...
               20   5   3   1   1    0    0    0    0    0    0;...
               53  20  12  10   7    5    5    5    2    1    1];

shen_sn_symp = [56.2 40.5 34.8 31.4 29.2 27.5 26.1 25.0 24.1 23.3 22.6];
shen_sp_symp = [72.0 53.1 46.2 42.1 39.2 37.1 35.5 34.1 32.9 31.9 31.1];
shen_sn_inf  = shen_sn_symp / 2;
shen_sp_inf  = shen_sp_symp / 2;
% Shen hosp: 상수 84.4%
shen_hosp    = 84.4;

% 시점 인덱스 (doc12와 동일)
idx1 = 1;    % year 1  → Sv (full)
idx2 = 2;    % year 2  → Sw (waned once)
idx3 = 10;   % year 10 → Sww (waned twice)

% ── Daniels ──────────────────────────────────────────────────
D.n_d1 = mean(dan_sn_symp(:,idx1))/100;
D.n_d2 = mean(dan_sn_symp(:,idx2))/100;
D.n_d3 = mean(dan_sn_symp(:,idx3))/100;
D.p_d1 = mean(dan_sp_symp(:,idx1))/100;
D.p_d2 = mean(dan_sp_symp(:,idx2))/100;
D.p_d3 = mean(dan_sp_symp(:,idx3))/100;
D.n_i1 = mean(dan_sn_inf(:,idx1))/100;
D.n_i2 = mean(dan_sn_inf(:,idx2))/100;
D.n_i3 = mean(dan_sn_inf(:,idx3))/100;
D.p_i1 = mean(dan_sp_inf(:,idx1))/100;
D.p_i2 = mean(dan_sp_inf(:,idx2))/100;
D.p_i3 = mean(dan_sp_inf(:,idx3))/100;
D.n_h1 = mean(dan_sn_hosp(:,idx1))/100;
D.n_h2 = mean(dan_sn_hosp(:,idx2))/100;
D.n_h3 = mean(dan_sn_hosp(:,idx3))/100;
D.p_h1 = mean(dan_sp_hosp(:,idx1))/100;
D.p_h2 = mean(dan_sp_hosp(:,idx2))/100;
D.p_h3 = mean(dan_sp_hosp(:,idx3))/100;

% ── Shen ─────────────────────────────────────────────────────
H.n_d1 = shen_sn_symp(idx1)/100;
H.n_d2 = shen_sn_symp(idx2)/100;
H.n_d3 = shen_sn_symp(idx3)/100;
H.p_d1 = shen_sp_symp(idx1)/100;
H.p_d2 = shen_sp_symp(idx2)/100;
H.p_d3 = shen_sp_symp(idx3)/100;
H.n_i1 = shen_sn_inf(idx1)/100;
H.n_i2 = shen_sn_inf(idx2)/100;
H.n_i3 = shen_sn_inf(idx3)/100;
H.p_i1 = shen_sp_inf(idx1)/100;
H.p_i2 = shen_sp_inf(idx2)/100;
H.p_i3 = shen_sp_inf(idx3)/100;
H.n_h1 = shen_hosp/100;
H.n_h2 = shen_hosp/100;
H.n_h3 = shen_hosp/100;
H.p_h1 = shen_hosp/100;
H.p_h2 = shen_hosp/100;
H.p_h3 = shen_hosp/100;

% ── Average ───────────────────────────────────────────────────
fields = fieldnames(D);
A = struct();
for fi = 1:numel(fields)
    A.(fields{fi}) = (D.(fields{fi}) + H.(fields{fi})) / 2;
end

VE_sets   = {D, H};
VE_names  = {'Daniels','Shen'};

fprintf('VE settings (n_d/p_d/n_i/p_i/n_h/p_h @ year1/2/10):\n')
for v = 1:2
    VE = VE_sets{v};
    fprintf('%s:\n', VE_names{v})
    fprintf('  SN: d=[%.3f %.3f %.3f]  i=[%.3f %.3f %.3f]  h=[%.3f %.3f %.3f]\n', ...
        VE.n_d1,VE.n_d2,VE.n_d3, VE.n_i1,VE.n_i2,VE.n_i3, VE.n_h1,VE.n_h2,VE.n_h3)
    fprintf('  SP: d=[%.3f %.3f %.3f]  i=[%.3f %.3f %.3f]  h=[%.3f %.3f %.3f]\n', ...
        VE.p_d1,VE.p_d2,VE.p_d3, VE.p_i1,VE.p_i2,VE.p_i3, VE.p_h1,VE.p_h2,VE.p_h3)
end

%% ── 2. 데이터 로드 ────────────────────────────────────────────
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
% mean_kappa = mean(kappa_draws, 2);
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
doses_annual = annual_national * doses_alph / sum(doses_alph);

%% ── 3. 공통 inputs 기본값 ─────────────────────────────────────
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
inputs_base.ohm_vac   = 1/(1);    % Sv → Sw
inputs_base.ohm_vac2  = 1/(4);    % Sw → Sww
inputs_base.ohm_vac3  = 0;   % Sww → S (fully waned)

h_ref=[0.0482 0.0428 0.0344 0.0256 0.0229 0.0237 0.0248 0.0257 ...
       0.0263 0.0274 0.0294 0.0323 0.0369 0.0441 0.0559 0.0721 0.1109]';
f_ref=[0.00137 0.000660 0.000516 0.000582 0.000538 0.000611 0.000707 0.000858 ...
       0.001023 0.001190 0.001453 0.001825 0.002244 0.003219 0.004579 0.007225 0.016830]';
inputs_base.h1=h_ref; inputs_base.h2=h_ref; inputs_base.h3=h_ref; inputs_base.h4=h_ref;
inputs_base.f1=f_ref; inputs_base.f2=f_ref; inputs_base.f3=f_ref; inputs_base.f4=f_ref;

%% ── 4. VE 셋팅별 inputs 생성 헬퍼 ───────────────────────────
function inp = apply_ve(inp, VE)
    % Sv (full, year 1): e_n_*, e_p_*
    inp.e_n_d = VE.n_d1; inp.e_n_i = VE.n_i1; inp.e_n_h = VE.n_h1; inp.e_n_m = VE.n_h1;
    inp.e_p_d = VE.p_d1; inp.e_p_i = VE.p_i1; inp.e_p_h = VE.p_h1; inp.e_p_m = VE.p_h1;
    % Sw (waned once, year 2): ew_n_*, ew_p_*
    inp.ew_n_d = VE.n_d2; inp.ew_n_i = VE.n_i2; inp.ew_n_h = VE.n_h2; inp.ew_n_m = VE.n_h2;
    inp.ew_p_d = VE.p_d2; inp.ew_p_i = VE.p_i2; inp.ew_p_h = VE.p_h2; inp.ew_p_m = VE.p_h2;
    % Sww (waned twice, year 10): eww_n_*, eww_p_*
    inp.eww_n_d = VE.n_d3; inp.eww_n_i = VE.n_i3; inp.eww_n_h = VE.n_h3; inp.eww_n_m = VE.n_h3;
    inp.eww_p_d = VE.p_d3; inp.eww_p_i = VE.p_i3; inp.eww_p_h = VE.p_h3; inp.eww_p_m = VE.p_h3;
end

%% ── 5. 헬퍼 함수 ─────────────────────────────────────────────
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

%% ── 6. 메인 루프: VE셋팅 × R0 q50 (baseline은 공통) ──────────
% Baseline은 VE와 무관하므로 q50 한 번만 계산
R0_sets  = {R0pred_q25, R0pred_q50, R0pred_q75};
R0_names = {'q25','q50','q75'};

% 결과: [S × A5 × 3(R0)] × 3(VE) × 2(baseline/vaccine)
Case_0 = zeros(S,A5,3);  % baseline은 VE 무관 → R0 3개만
Hosp_0 = zeros(S,A5,3);
Dead_0 = zeros(S,A5,3);
DALY_0 = zeros(S,A5,3);

% [S × A5 × 3(R0) × 3(VE)]
Case_V = zeros(S,A5,3,2);
Hosp_V = zeros(S,A5,3,2);
Dead_V = zeros(S,A5,3,2);
DALY_V = zeros(S,A5,3,2);

popbystate = zeros(S,1);

inp0 = inputs_base; inp0.vacc=0; inp0.vacorder=1;
% VE 없이 baseline 먼저 (VE 무관하므로 VE=D로 apply해도 어차피 vacc=0)
inp0 = apply_ve(inp0, D);

for rc = 1:3
    R0pred = R0_sets{rc};
    fprintf('\nBaseline R0=%s...\n', R0_names{rc})
    for si = 1:S
        [Case_0(si,:,rc),Hosp_0(si,:,rc),Dead_0(si,:,rc),DALY_0(si,:,rc)] = ...
            run_sim(si,inp0, ...
                    op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                    op_seroposall_3,op_seroposall_4, ...
                    pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
        if rc == 2
            pop_2025 = pop_mat((si-1)*T+T,:)';
            A0tmp = median(squeeze(op_seroposall_0(:,si,T,:)),1)' .* pop_2025;
            A1tmp = median(squeeze(op_seroposall_1(:,si,T,:)),1)' .* pop_2025;
            A2tmp = median(squeeze(op_seroposall_2(:,si,T,:)),1)' .* pop_2025;
            A3tmp = median(squeeze(op_seroposall_3(:,si,T,:)),1)' .* pop_2025;
            A4tmp = median(squeeze(op_seroposall_4(:,si,T,:)),1)' .* pop_2025;
            popbystate(si) = sum([A0tmp;A1tmp;A2tmp;A3tmp;A4tmp],'all');
        end
    end
end

for vc = 1:2
    VE = VE_sets{vc};
    inpV = apply_ve(inputs_base, VE);
    for rc = 1:3
        R0pred = R0_sets{rc};
        fprintf('\nVE=%s R0=%s...\n', VE_names{vc}, R0_names{rc})
        for si = 1:S
            inpV.vacc     = doses_annual(si);
            inpV.vacorder = 1;
            [Case_V(si,:,rc,vc),Hosp_V(si,:,rc,vc),Dead_V(si,:,rc,vc),DALY_V(si,:,rc,vc)] = ...
                run_sim(si,inpV, ...
                        op_seroposall_0,op_seroposall_1,op_seroposall_2, ...
                        op_seroposall_3,op_seroposall_4, ...
                        pop_mat,R0pred,data_lifeexp_age,tspan,T,A5);
        end
    end
end

save('results/SIR_vaccine_VE_sensitivity.mat', ...
     'Case_0','Hosp_0','Dead_0','DALY_0', ...
     'Case_V','Hosp_V','Dead_V','DALY_V', ...
     'doses_annual','popbystate','state_codes','VE_names')
fprintf('Saved: results/SIR_vaccine_VE_sensitivity.mat\n')

%% ── 7. BoD averted: 전국 합산 + 콘솔 출력 ───────────────────
total_doses_5yr = annual_national*5;
vactot = total_doses_5yr/1000;

fprintf('\n━━━ BoD averted summary (5M doses/yr × 5yr, q50, Brazil allocation) ━━━\n')
for vc = 1:2
    avt_case = sum(Case_0(:,:,2) - Case_V(:,:,2,vc), 1);  % q50
    avt_hosp = sum(Hosp_0(:,:,2) - Hosp_V(:,:,2,vc), 1);
    avt_dead = sum(Dead_0(:,:,2) - Dead_V(:,:,2,vc), 1);
    avt_daly = sum(DALY_0(:,:,2) - DALY_V(:,:,2,vc), 1);
    fprintf('\n[%s]\n', VE_names{vc})
    fprintf('  Cases:            %9.0f  (per 1K vac: %.3f)\n', sum(avt_case), sum(avt_case)/vactot)
    fprintf('  Hospitalizations: %9.0f  (per 1K vac: %.4f)\n', sum(avt_hosp), sum(avt_hosp)/vactot)
    fprintf('  Deaths:           %9.4f  (per 1K vac: %.6f)\n', sum(avt_dead), sum(avt_dead)/vactot)
    fprintf('  DALYs:            %9.4f  (per 1K vac: %.6f)\n', sum(avt_daly), sum(avt_daly)/vactot)
end

%% ── 8. Figure: VE 3가지 비교 (q50, 연령별) ─────────────────
mkdir('figures')

dark_green  = [0.0 0.55 0.0];
light_green = [0.6 1.0  0.6];
shade_x = [5 15 15 5];

col_dan  = [0, 62, 116]/255;
col_shen = [225, 36, 42]/255;
col_ve   = [col_dan; col_shen];

ylabels_avt = {'Cases averted','Hospitalizations averted', ...
               'Deaths averted','DALYs averted'};
outcome_keys = {'case','hosp','dead','daly'};

figure('Units','centimeters','Position',[0 0 18 12],'PaperPositionMode','auto');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

for m = 1:4
    nexttile(m)
    hold on
    fill(shade_x,[0 0 1 1]*Inf,light_green,'EdgeColor','none','FaceAlpha',0.4)

    for vc = 1:2
        switch m
            case 1; q50 = sum(Case_0(:,:,2)-Case_V(:,:,2,vc),1);
                    q25 = sum(Case_0(:,:,1)-Case_V(:,:,1,vc),1);
                    q75 = sum(Case_0(:,:,3)-Case_V(:,:,3,vc),1);
            case 2; q50 = sum(Hosp_0(:,:,2)-Hosp_V(:,:,2,vc),1);
                    q25 = sum(Hosp_0(:,:,1)-Hosp_V(:,:,1,vc),1);
                    q75 = sum(Hosp_0(:,:,3)-Hosp_V(:,:,3,vc),1);
            case 3; q50 = sum(Dead_0(:,:,2)-Dead_V(:,:,2,vc),1);
                    q25 = sum(Dead_0(:,:,1)-Dead_V(:,:,1,vc),1);
                    q75 = sum(Dead_0(:,:,3)-Dead_V(:,:,3,vc),1);
            case 4; q50 = sum(DALY_0(:,:,2)-DALY_V(:,:,2,vc),1);
                    q25 = sum(DALY_0(:,:,1)-DALY_V(:,:,1,vc),1);
                    q75 = sum(DALY_0(:,:,3)-DALY_V(:,:,3,vc),1);
        end
        offset = (vc-2)*1.5;  % 살짝 x 오프셋으로 겹침 방지
        XXX(vc) = errorbar(x_age+offset, q50, q50-q25, q75-q50, 'o', ...
            'Color',col_ve(vc,:),'MarkerFaceColor',col_ve(vc,:), ...
            'MarkerSize',3,'LineWidth',1.0,'CapSize',3);
    end

    yline(0,'k:','LineWidth',0.8)
    xlim([-2.5 82.5]); xticks(x_age); xticklabels(age_labels)
    ylabel(ylabels_avt{m})
    set(gca,'FontName','SansSerif','FontSize',9)
    if m == 1
        lge = legend([XXX], VE_names,'Location','northeast','FontSize',8); legend('boxoff');
    end
end
title(lge, 'VE setting','FontWeight','normal')
% sgtitle('BoD averted: VE sensitivity (5M doses/yr × 5yr, vacorder=1, Brazil allocation)', ...
%     'FontSize',10,'FontName','SansSerif')
exportgraphics(gcf,'figures/fig_BoD_averted_VE_sensitivity.jpg','Resolution',600)
fprintf('Saved: figures/fig_BoD_averted_VE_sensitivity.jpg\n')

%% ── 8-new. Figure: VE 비교 (구 레이아웃 + point/interval 스타일) ──
%% ── 8-new. Figure: VE 비교 (구 레이아웃 + point/interval 스타일) ──
mkdir('figures')

light_green = [0.6 1.0 0.6];
shade_x = [2.5 12.5 12.5 2.5];   % 소아 타겟 연령대 강조 (기존과 동일)

col_dan  = [0, 62, 116]/255;
col_shen = [225, 36, 42]/255;
col_ve   = [col_dan; col_shen];

function plot_avt_panel(metric_0, metric_V, vactot, col_ve, VE_names, ...
                         age_labels, x_age, shade_x, ylab, fname)
    nVE = size(metric_V,4);
light_green = [0.6 1.0 0.6];
    q50 = zeros(nVE, numel(x_age)); q25 = q50; q75 = q50;
    tot50 = zeros(nVE,1); tot25 = tot50; tot75 = tot50;
    for vc = 1:nVE
        q50(vc,:) = sum(metric_0(:,:,2) - metric_V(:,:,2,vc), 1) / vactot;
        q25(vc,:) = sum(metric_0(:,:,1) - metric_V(:,:,1,vc), 1) / vactot;
        q75(vc,:) = sum(metric_0(:,:,3) - metric_V(:,:,3,vc), 1) / vactot;
        tot50(vc) = sum(q50(vc,:));
        tot25(vc) = sum(q25(vc,:));
        tot75(vc) = sum(q75(vc,:));
    end
    ylimval = max(q75(:)) * 1.15;

    figure('Units','centimeters','Position',[0 0 9 5],'PaperPositionMode','auto');
    tiledlayout(1,4,'TileSpacing','compact')

    nexttile([1,3])
    fill(shade_x, [0 0 ylimval ylimval], light_green, 'EdgeColor','none','FaceAlpha',0.4)
    hold on
    bb = bar(x_age, q50', 'EdgeColor','none');
    for vc = 1:nVE
        bb(vc).FaceColor = col_ve(vc,:);
    end
    for vc = 1:nVE
        offset = bb(vc).XEndPoints - x_age;   % 실제 bar 중심 x좌표 오프셋
        errorbar(x_age+offset, q50(vc,:), q50(vc,:)-q25(vc,:), q75(vc,:)-q50(vc,:), ...
            '.', 'Color', 'k', 'LineWidth', 1.0, 'CapSize', 3, 'LineStyle','none')
    end
    yline(0,'k:','LineWidth',0.8)
    xticks(x_age); xticklabels(age_labels); xlim([-2.5 82.5]); ylim([0 ylimval])
    ylabel(ylab)
    set(gca,'FontName','SansSerif','FontSize',10)
    legend(bb, VE_names, 'Location','northeast','FontSize',8); legend('boxoff')

    nexttile(4)
    hold on
    bbt = bar(1, tot50(:)', 'EdgeColor','none');
    for vc = 1:nVE
        bbt(vc).FaceColor = col_ve(vc,:);
        offset = bbt(vc).XEndPoints(1) - 1;
        errorbar(1+offset, tot50(vc), tot50(vc)-tot25(vc), tot75(vc)-tot50(vc), ...
            '.', 'Color', 'k', 'LineWidth', 1.0, 'CapSize', 3, 'LineStyle','none')
    end
    xlim([0.5 1.5]); xticks(1); xticklabels({'Total'})
    ylim([0 max(tot75)*1.4])
    set(gca,'FontName','SansSerif','FontSize',10,'XTickLabelRotation',0)

    exportgraphics(gcf, fname, 'Resolution', 600)
    fprintf('Saved: %s\n', fname)
end

plot_avt_panel(Case_0, Case_V, vactot, col_ve, VE_names, ...
    age_labels, x_age, shade_x, 'Cases averted', ...
    'figures/fig_vac_byage_case.jpg')

plot_avt_panel(Hosp_0, Hosp_V, vactot, col_ve, VE_names, ...
    age_labels, x_age, shade_x, 'Hospitalizations averted', ...
    'figures/fig_vac_byage_hosp.jpg')

plot_avt_panel(Dead_0, Dead_V, vactot, col_ve, VE_names, ...
    age_labels, x_age, shade_x, 'Fatality averted', ...
    'figures/fig_vac_byage_dead.jpg')

plot_avt_panel(DALY_0, DALY_V, vactot, col_ve, VE_names, ...
    age_labels, x_age, shade_x, 'DALY averted', ...
    'figures/fig_vac_byage_daly.jpg')

%% ═══════════════════════════════════════════════════════════════
%  CSV 출력: VE 셋팅(Daniels/Shen)별 averted case/hosp/dead/daly
%  — 연령별 + 전국합계, raw count와 per-1000-vaccinated 둘 다
%  median [IQR] (R0 q25/q50/q75 시나리오 기반)
%% ═══════════════════════════════════════════════════════════════
mkdir('result_csv')

fmt_iqr = @(m,lo,hi) sprintf('%.2f [%.2f, %.2f]', m, lo, hi);

metric_names = {'case','hosp','dead','daly'};
metric_0_set = {Case_0, Hosp_0, Dead_0, DALY_0};
metric_V_set = {Case_V, Hosp_V, Dead_V, DALY_V};

row_labels = [age_labels, {'Total'}];   % 17 연령군 + 전국합계

for vc = 1:2
    col_names = {};
    col_data  = {};

    for mi = 1:4
        Metric_0 = metric_0_set{mi};
        Metric_V = metric_V_set{mi};

        % [3(rc) x A5] — 연령별, R0 시나리오별 averted (raw count, 주 전체 합산)
        avt_age = zeros(3, A5);
        for rc = 1:3
            avt_age(rc,:) = sum(Metric_0(:,:,rc) - Metric_V(:,:,rc,vc), 1);
        end
        avt_total = sum(avt_age, 2);   % [3 x 1], 전국+전연령 합계

        raw_col   = cell(A5+1,1);
        per1k_col = cell(A5+1,1);
        for a = 1:A5
            raw_col{a}   = fmt_iqr(avt_age(2,a), avt_age(1,a), avt_age(3,a));
            per1k_col{a} = fmt_iqr(avt_age(2,a)/vactot, avt_age(1,a)/vactot, avt_age(3,a)/vactot);
        end
        raw_col{A5+1}   = fmt_iqr(avt_total(2), avt_total(1), avt_total(3));
        per1k_col{A5+1} = fmt_iqr(avt_total(2)/vactot, avt_total(1)/vactot, avt_total(3)/vactot);

        col_names = [col_names, {sprintf('%s_averted', metric_names{mi})}, ...
                                 {sprintf('%s_averted_per1000vac', metric_names{mi})}];
        col_data  = [col_data, {raw_col}, {per1k_col}];
    end

    tbl = [col_data{:}];
    T_ve = cell2table(tbl, 'VariableNames', col_names, 'RowNames', row_labels);
    writetable(T_ve, sprintf('result_csv/VE_averted_%s.csv', VE_names{vc}), 'WriteRowNames', true)
end

fprintf('Saved: result_csv/VE_averted_Daniels.csv, result_csv/VE_averted_Shen.csv\n')