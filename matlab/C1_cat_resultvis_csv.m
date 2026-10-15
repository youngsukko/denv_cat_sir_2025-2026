%% ═══════════════════════════════════════════════════════════════
%  논문 부록용 표 일괄 CSV 출력 (result_csv/)
%  포맷: "median [lo, hi]" — FOI/R0는 IQR(25-75%), 나머지는 95% CrI(2.5-97.5%)
%  ※ 메인 스크립트(위)에서 계산된 변수들을 그대로 사용함. 이 블록은
%    lambda_adj, lambda_pre_adj, kappa_adj, omega_draws, delta_draws,
%    rho_draws, op_seroposall_pre_0..4, op_seroposall_0..4,
%    ratio_med/ratio_lo/ratio_hi, state_codes, S, T, A5, agelabels,
%    pct025/pct975/pct25/pct75 가 이미 워크스페이스에 있다는 전제로 돌아감.
%% ═══════════════════════════════════════════════════════════════
mkdir('result_csv')

fmt_pct  = @(m,lo,hi) sprintf('%.2f [%.0f, %.0f]', m*100, lo*100, hi*100);  % (미사용, 참고용)
fmt_pct_nl = @(m,lo,hi) sprintf('%.2f\n[%.0f, %.0f]', m*100, lo*100, hi*100); % FOI, rho: median 다음 줄바꿈, ×100
fmt_iqr  = @(m,lo,hi) sprintf('%.2f [%.2f, %.2f]', m, lo, hi);              % 나머지: 전부 2dp, 스케일 그대로

L_exp = 76.4;   % R0 = FOI(비×4) * L_exp + 1, section 19와 동일 컨벤션

% ── 지역순 재배열 (North→Northeast→Central-west→Southeast→South, 지역 내 주번호순) ──
north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];
reg_order   = [north_order, northeast_order, centwest_order, southeast_order, south_order];
codes_reord = state_codes(reg_order);

period_labels = [{'pre-00'}, arrayfun(@(y) sprintf('%02d', mod(y,100)), 2000:2025, 'UniformOutput', false)];
year_labels   = arrayfun(@(y) sprintf('%02d', mod(y,100)), 2000:2025, 'UniformOutput', false);

%% ── 1. FOI table (27 states × 27 periods, IQR, ×100) ────────────
FOI_tbl = cell(S, numel(period_labels));
for k = 1:S
    si = reg_order(k);
    vals_pre = 4*lambda_pre_adj(:,si);
    FOI_tbl{k,1} = fmt_pct_nl(pmed(vals_pre), pct25(vals_pre), pct75(vals_pre));
    for ti = 1:T
        vals = 4*squeeze(lambda_adj(:,si,ti));
        FOI_tbl{k,ti+1} = fmt_pct_nl(pmed(vals), pct25(vals), pct75(vals));
    end
end
T_foi = cell2table(FOI_tbl, 'VariableNames', period_labels, 'RowNames', codes_reord);
writetable(T_foi, 'result_csv/FOI_table.csv', 'WriteRowNames', true)

%% ── 2. R0 table (27 states × 27 periods, IQR) ───────────────────
% pre-00은 원 스크립트에 없던 값이라 동일 공식(L_exp, ×4 없이)으로 새로 계산한 것.
R0_tbl = cell(S, numel(period_labels));
for k = 1:S
    si = reg_order(k);
    vals_pre = lambda_pre_adj(:,si)*L_exp + 1;
    R0_tbl{k,1} = fmt_iqr(pmed(vals_pre), pct25(vals_pre), pct75(vals_pre));
    for ti = 1:T
        vals = squeeze(lambda_adj(:,si,ti))*L_exp + 1;
        R0_tbl{k,ti+1} = fmt_iqr(pmed(vals), pct25(vals), pct75(vals));
    end
end
T_r0 = cell2table(R0_tbl, 'VariableNames', period_labels, 'RowNames', codes_reord);
writetable(T_r0, 'result_csv/R0_table.csv', 'WriteRowNames', true)

%% ── 3. Report rate (rho) table (27 states × 26 년, 95% CrI) ────
rho_tbl = cell(S, T);
for k = 1:S
    si = reg_order(k);
    for ti = 1:T
        vals = squeeze(rho_draws(:,si,ti));
        rho_tbl{k,ti} = fmt_pct_nl(pmed(vals), pct025(vals), pct975(vals));
    end
end
T_rho = cell2table(rho_tbl, 'VariableNames', year_labels, 'RowNames', codes_reord);
writetable(T_rho, 'result_csv/rho_table.csv', 'WriteRowNames', true)

%% ── 4/5. 면역 프로파일 (연령군 인구비 정규화 seropositivity), 2000/2025 ──
imm_2000 = cell(S, A5);
imm_2025 = cell(S, A5);
for k = 1:S
    si = reg_order(k);
    X0b = squeeze(op_seroposall_pre_0(:,si,:));
    X1b = squeeze(op_seroposall_pre_1(:,si,:));
    X2b = squeeze(op_seroposall_pre_2(:,si,:));
    X3b = squeeze(op_seroposall_pre_3(:,si,:));
    X4b = squeeze(op_seroposall_pre_4(:,si,:));
    denom = X0b+X1b+X2b+X3b+X4b;
    seropos = 1 - X0b./denom;   % [n_draw x A5]
    for a = 1:A5
        vals = seropos(:,a);
        imm_2000{k,a} = fmt_iqr(pmed(vals), pct025(vals), pct975(vals));
    end

    X0b=squeeze(op_seroposall_0(:,si,T,:));
    X1b=squeeze(op_seroposall_1(:,si,T,:));
    X2b=squeeze(op_seroposall_2(:,si,T,:));
    X3b=squeeze(op_seroposall_3(:,si,T,:));
    X4b=squeeze(op_seroposall_4(:,si,T,:));
    denom = X0b+X1b+X2b+X3b+X4b;
    seropos = 1 - X0b./denom;
    for a = 1:A5
        vals = seropos(:,a);
        imm_2025{k,a} = fmt_iqr(pmed(vals), pct025(vals), pct975(vals));
    end
end
T_imm2000 = cell2table(imm_2000, 'VariableNames', agelabels, 'RowNames', codes_reord);
T_imm2025 = cell2table(imm_2025, 'VariableNames', agelabels, 'RowNames', codes_reord);
writetable(T_imm2000, 'result_csv/immune_profile_2000.csv', 'WriteRowNames', true)
writetable(T_imm2025, 'result_csv/immune_profile_2025.csv', 'WriteRowNames', true)

%% ── 4b/5b. 주별 인구가중 면역 프로파일 (0~4회 감염, 2000/2025, ×100) ──
fmt_iqr_pct_nl = @(m,lo,hi) sprintf('%.2f\n[%.2f, %.2f]', m*100, lo*100, hi*100);

exposure_labels = {'0','1','2','3','4'};
imm_state_cols = [strcat('2000_', exposure_labels), strcat('2025_', exposure_labels)];  % 10열

imm_state_tbl = cell(S, 10);
for k = 1:S
    si = reg_order(k);
    pop_2000_si = pop_mat((si-1)*T+1, :);   % [1 x A5]
    pop_2025_si = pop_mat((si-1)*T+T, :);

    Xpre = {squeeze(op_seroposall_pre_0(:,si,:)), squeeze(op_seroposall_pre_1(:,si,:)), ...
            squeeze(op_seroposall_pre_2(:,si,:)), squeeze(op_seroposall_pre_3(:,si,:)), ...
            squeeze(op_seroposall_pre_4(:,si,:))};
    denom_pre = Xpre{1}+Xpre{2}+Xpre{3}+Xpre{4}+Xpre{5};
    for kx = 1:5
        prop_k = Xpre{kx}./denom_pre;                       % [n_draw x A5]
        vals = (prop_k * pop_2000_si') / sum(pop_2000_si);  % [n_draw x1], 주 내 연령가중
        imm_state_tbl{k,kx} = fmt_iqr_pct_nl(pmed(vals), pct025(vals), pct975(vals));
    end

    Xcur = {squeeze(op_seroposall_0(:,si,T,:)), squeeze(op_seroposall_1(:,si,T,:)), ...
            squeeze(op_seroposall_2(:,si,T,:)), squeeze(op_seroposall_3(:,si,T,:)), ...
            squeeze(op_seroposall_4(:,si,T,:))};
    denom_cur = Xcur{1}+Xcur{2}+Xcur{3}+Xcur{4}+Xcur{5};
    for kx = 1:5
        prop_k = Xcur{kx}./denom_cur;
        vals = (prop_k * pop_2025_si') / sum(pop_2025_si);
        imm_state_tbl{k,5+kx} = fmt_iqr_pct_nl(pmed(vals), pct025(vals), pct975(vals));
    end
end
T_imm_state = cell2table(imm_state_tbl, 'VariableNames', imm_state_cols, 'RowNames', codes_reord);
writetable(T_imm_state, 'result_csv/immune_profile_bystate.csv', 'WriteRowNames', true)

%% ── 6. kappa_age (17행, 95% CrI) ────────────────────────────────
kappa_tbl = cell(A5,1);
for a = 1:A5
    vals = kappa_adj(:,a);
    kappa_tbl{a} = fmt_iqr(pmed(vals), pct025(vals), pct975(vals));
end
T_kappa = cell2table(kappa_tbl, 'VariableNames', {'kappa_age'}, 'RowNames', agelabels);
writetable(T_kappa, 'result_csv/kappa_age.csv', 'WriteRowNames', true)

%% ── 7. omega_age (17행, 95% CrI) ────────────────────────────────
omega_tbl = cell(A5,1);
for a = 1:A5
    vals = omega_draws(:,a);
    omega_tbl{a} = fmt_iqr(pmed(vals), pct025(vals), pct975(vals));
end
T_omega = cell2table(omega_tbl, 'VariableNames', {'omega_age'}, 'RowNames', agelabels);
writetable(T_omega, 'result_csv/omega_age.csv', 'WriteRowNames', true)

%% ── 8. delta (4행, 95% CrI) ──────────────────────────────────────
delta_labels = {'Primary','Secondary','Tertiary','Quaternary'};
delta_tbl = cell(4,1);
for k = 1:4
    vals = delta_draws(:,k);
    delta_tbl{k} = fmt_iqr(pmed(vals), pct025(vals), pct975(vals));
end
T_delta = cell2table(delta_tbl, 'VariableNames', {'delta'}, 'RowNames', delta_labels);
writetable(T_delta, 'result_csv/delta.csv', 'WriteRowNames', true)

%% ── 9. Case/Infection ratio (27행, 95% CrI) — 이미 계산된 값 사용 ──
ratio_tbl = cell(S,1);
for k = 1:S
    si = reg_order(k);
    ratio_tbl{k} = fmt_iqr(ratio_med(si), ratio_lo(si), ratio_hi(si));
end
T_ratio = cell2table(ratio_tbl, 'VariableNames', {'case_infection_ratio'}, 'RowNames', codes_reord);
writetable(T_ratio, 'result_csv/case_infection_ratio.csv', 'WriteRowNames', true)

fprintf('Saved 9 CSV files to result_csv/\n')

%% ── 10/11. 전국(인구가중) 연도별 FOI·보고율 — median 다음 줄바꿈 후 구간 ──
n_draw_here = size(lambda_adj,1);
nat_FOI = cell(1, numel(period_labels));
nat_rho = cell(1, numel(year_labels));

% pre-00: 인구가중은 2000년 인구 사용
pop_2000_bystate = zeros(S,1);
pop_2025_bystate = zeros(S,1);
for si = 1:S
    pop_2000_bystate(si) = sum(pop_mat((si-1)*T+1, :));
    pop_2025_bystate(si) = sum(pop_mat((si-1)*T+T, :));
end

vals_pre = zeros(n_draw_here,1);
for si = 1:S
    vals_pre = vals_pre + 4*lambda_pre_adj(:,si) * pop_2000_bystate(si);
end
vals_pre = vals_pre / sum(pop_2000_bystate);
nat_FOI{1,1} = fmt_pct_nl(pmed(vals_pre), pct25(vals_pre), pct75(vals_pre));

for ti = 1:T
    pop_ti_bystate = zeros(S,1);
    for si = 1:S
        pop_ti_bystate(si) = sum(pop_mat((si-1)*T+ti, :));
    end
    vals_foi = zeros(n_draw_here,1);
    vals_rho = zeros(n_draw_here,1);
    for si = 1:S
        vals_foi = vals_foi + 4*squeeze(lambda_adj(:,si,ti)) * pop_ti_bystate(si);
        vals_rho = vals_rho + squeeze(rho_draws(:,si,ti)) * pop_ti_bystate(si);
    end
    vals_foi = vals_foi / sum(pop_ti_bystate);
    vals_rho = vals_rho / sum(pop_ti_bystate);
    nat_FOI{1,ti+1} = fmt_pct_nl(pmed(vals_foi), pct25(vals_foi), pct75(vals_foi));
    nat_rho{1,ti}   = fmt_pct_nl(pmed(vals_rho), pct025(vals_rho), pct975(vals_rho));
end

T_natFOI = cell2table(nat_FOI, 'VariableNames', period_labels, 'RowNames', {'Brazil'});
T_natRho = cell2table(nat_rho, 'VariableNames', year_labels,   'RowNames', {'Brazil'});
writetable(T_natFOI, 'result_csv/FOI_national_byyear.csv', 'WriteRowNames', true)
writetable(T_natRho, 'result_csv/rho_national_byyear.csv', 'WriteRowNames', true)

fprintf('Saved 2 more CSV files (national FOI/rho by year) to result_csv/\n')