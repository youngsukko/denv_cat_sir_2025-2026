%% ═══════════════════════════════════════════════════════════════
%  CSV 출력: 주별 백신분배량, 주별 연간 BoD per 100K (median [IQR])
%  ※ reg_order/codes_reord는 위 섹션 5에서 이미 정의됨을 전제로 함
%% ═══════════════════════════════════════════════════════════════
mkdir('result_csv')

fmt_iqr = @(m,lo,hi) sprintf('%.2f [%.2f, %.2f]', m, lo, hi);

%% ── 1. 주별 백신분배량 (27×2): 총 도즈, 어린이 타겟인구 1000명당 ──
opt = load('results/SIR_opt_result.mat', 'doses_annual_brazil');
doses_annual_brazil = opt.doses_annual_brazil;   % [S x 1]

children_ages = [2,3];   % 5-14세
target_pop_children = zeros(S,1);
for si = 1:S
    pop_2025_si = pop_mat((si-1)*T+T, :);
    target_pop_children(si) = sum(pop_2025_si(children_ages));
end

doses_per1000_children = doses_annual_brazil ./ target_pop_children * 1000;

vac_tbl = cell(S,2);
for k = 1:S
    si = reg_order(k);
    vac_tbl{k,1} = sprintf('%.0f', doses_annual_brazil(si));
    vac_tbl{k,2} = sprintf('%.2f', doses_per1000_children(si));
end
T_vac = cell2table(vac_tbl, 'VariableNames', {'doses_total','doses_per1000_children'}, ...
    'RowNames', codes_reord);
writetable(T_vac, 'result_csv/vaccine_allocation_bystate.csv', 'WriteRowNames', true)

%% ── 2. 주별 연간 BoD per 100K (27×4): median [IQR] ──────────────
mats_3   = {Case_age_3, Hosp_age_3, Dead_age_3, DALY_age_3};
bod_cols = {'case_per100k','hosp_per100k','fatal_per100k','daly_per100k'};

bod_tbl = cell(S, 4);
for m = 1:4
    v50 = sum(mats_3{m}(:,:,2),2)./popbystate*1e5/5;
    v25 = sum(mats_3{m}(:,:,1),2)./popbystate*1e5/5;
    v75 = sum(mats_3{m}(:,:,3),2)./popbystate*1e5/5;
    for k = 1:S
        si = reg_order(k);
        bod_tbl{k,m} = fmt_iqr(v50(si), v25(si), v75(si));
    end
end
T_bod = cell2table(bod_tbl, 'VariableNames', bod_cols, 'RowNames', codes_reord);
writetable(T_bod, 'result_csv/BoD_per100k_bystate.csv', 'WriteRowNames', true)

fprintf('Saved: result_csv/vaccine_allocation_bystate.csv, result_csv/BoD_per100k_bystate.csv\n')