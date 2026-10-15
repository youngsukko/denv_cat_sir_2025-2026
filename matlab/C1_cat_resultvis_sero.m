%% ── helper: validation 플롯 (원본 scpaux_serovaliplot_* 5종 통합 재현) ──
% obs_style: 'step'(원본 scpaux_serovaliplot / pre2000m)
%            'range'(원본 scpaux_serovaliplot_pre2000_range — 전체구간 단일밴드)
%            'multibin'(원본 scpaux_serovaliplot_age_multi / 2000_agerange — 구간별 밴드)
% ti = 0  → pre-2000 baseline(seropos_pre_0) 사용, 그 외 → post-2000 t시점(op_seroposall_0) 사용
function draw_vali(seropos_pre_0, op_seroposall_0, si, ti, agerangeidx, ...
                    tmp_dat_ages, tmp_prop, minidx, maxidx, ...
                    col_model, col_model_light, obs_color, obs_style, title_str)

    if ti == 0
        p_band = squeeze(1 - seropos_pre_0(:,si,:));
    else
        X0b = squeeze(op_seroposall_0(:,si,ti,:));
        p_band = 1 - X0b;
    end
    p_s = sort(p_band);
    fill([agerangeidx, fliplr(agerangeidx)], [p_s(maxidx,:), fliplr(p_s(minidx,:))], ...
         col_model_light,'EdgeColor','none','FaceAlpha',0.7); hold on
    plot(agerangeidx, median(p_s), '-','color',col_model,'LineWidth',1.5)

    switch obs_style
        case 'step'
            % 원본 scpaux_serovaliplot / pre2000m: 구간 경계 그대로, 선형인덱싱 tmp_prop(i)
            for i = 1:length(tmp_dat_ages)-1
                plot([tmp_dat_ages(i) tmp_dat_ages(i+1)], [tmp_prop(i) tmp_prop(i)], ...
                     '-', 'Color', obs_color, 'LineWidth', 1.5)
            end
        case 'multibin'
            % 원본 scpaux_serovaliplot_age_multi / 2000_agerange: 구간별 min-max 밴드 + 점선 중앙값
            for i = 1:length(tmp_dat_ages)-1
                col_i = tmp_prop(:,i);
                fill([tmp_dat_ages(i) tmp_dat_ages(i+1) tmp_dat_ages(i+1) tmp_dat_ages(i)], ...
                     [min(col_i) min(col_i) max(col_i) max(col_i)], obs_color, ...
                     'EdgeColor','none','FaceAlpha',0.2)
                plot([tmp_dat_ages(i) tmp_dat_ages(i+1)], [median(col_i) median(col_i)], ...
                     ':', 'Color', obs_color, 'LineWidth', 1.5)
            end
        case 'range'
            % 원본 scpaux_serovaliplot_pre2000_range: tmp_dat_ages 무시, 전체구간 단일 min-max 밴드
            fill([agerangeidx(1) agerangeidx(end) agerangeidx(end) agerangeidx(1)], ...
                 [min(tmp_prop) min(tmp_prop) max(tmp_prop) max(tmp_prop)], obs_color, ...
                 'EdgeColor','none','FaceAlpha',0.2)
            plot([agerangeidx(1) agerangeidx(end)], [median(tmp_prop) median(tmp_prop)], ...
                 ':', 'Color', obs_color, 'LineWidth', 1.5)
    end

    % xlabel('Age'); ylabel('Seroprevalence')
    xlim([0 80]); ylim([0 1])
    title(title_str,'FontSize',9,'FontWeight','normal')
set(gca, "XTick",[0:10:80],"FontName",'sansserif','FontSize',12)
end

%% ── 18. Pre-2000/Post-2000 검증 (validation) — 원본 코드 완전 재현 ──
col_g = [0 1 0];   % 원본 'g'
col_m = [0 1 0];   % 원본 'g'

%% -- North: 2 panels --
fprintf('Fig: validation North...\n');
figure('Units','centimeters','Position',[0 0 5 8],'PaperPositionMode','auto');
tiledlayout(2,1,'TileSpacing','tight')
nexttile(1)  % PA (1997) — pre2000m
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(5), 0, agerangeidx, ...
    [0 5 10 15 25 35 45 55 82], [26.7 37.5 43.3 41.3 45.6 43.3 37.4 52.5]/100, ...
    minidx, maxidx, col_model, col_model_light, col_m, 'step', 'PA (1997)')
set(gca, "XTick",[0:10:80],"FontName",'sansserif','FontSize',12)
nexttile(2)  % AM (2014) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(3), 14, agerangeidx, ...
    [0 17 48 80], [41.9 92.4 99.5]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'AM (2014)')
set(gca, "XTick",[0:10:80],"FontName",'sansserif','FontSize',12)
exportgraphics(gcf,'figures/fig_vali_north.jpg','Resolution',600)

%% -- Northeast: 14 panels, 4x4 grid (원본과 동일) --
fprintf('Fig: validation Northeast...\n');
figure('Units','centimeters','Position',[0 0 15 20],'PaperPositionMode','auto');
tiledlayout(5,3,'TileSpacing','tight')

nexttile(1)  % CE (1994) — pre2000_range
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(10), 0, agerangeidx, ...
    [0 80], [71 50 21 25 31 55 23 35 58]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'range', 'CE (1994)')

nexttile(2)  % MA (1996) — pre2000m
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(8), 0, agerangeidx, ...
    [0 15 25 35 45 55 80], 1-[76.4 60 59 54.4 57.3 44.8]/100, minidx, maxidx, ...
    col_model, col_model_light, col_m, 'step', 'MA (1996)')
set(gca, "XTick",[0:10:80],"FontName",'sansserif','FontSize',12)

nexttile(3)  % BA (1998) main — pre2000_range
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(16), 0, agerangeidx, ...
    [0 80], [62.8 69.2 78.4]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'range', 'BA (1998)')
set(gca, "XTick",[0:10:80],"FontName",'sansserif','FontSize',12)

nexttile(4)  % BA (1998) infant — pre2000m
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(16), 0, agerangeidx, ...
    [0 3], [26.6]/100, minidx, maxidx, ...
    col_model, col_model_light, col_m, 'step', 'BA (1998)')

nexttile(5)  % PB (2005) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(12), 5, agerangeidx, ...
    [5 15 65], [84 96.8; 78.8 94.1; 69.9 77.7]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'PB (2005)')

nexttile(6)  % BA (2010) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(16), 10, agerangeidx, ...
    [5 10 15 80], [43.9 76.9 95.9]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'BA (2010)')

nexttile(7)  % CE (2012) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(10), 12, agerangeidx, ...
    [5 13], [46.6 46.6; 61.4 61.4]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'CE (2012)')

nexttile(8)  % PE (2015) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(13), 15, agerangeidx, ...
    [21.5-3, 21.5+3], [87.1 87.1]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'PE (2015)')

nexttile(9)  % BA (2015, Rodriguez-Barraquer 2019) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(16), 15, agerangeidx, ...
    [0 5 10 15 20 30 40 50], [57 67 77 85 90 92 95]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'BA (2015)')

nexttile(10)  % RN (2016) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(11), 16, agerangeidx, ...
    [0 17 48 80], [65.7 99 100]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'RN (2016)')

nexttile(11)  % BA (2017) — 원본이 age_multi 아니라 serovaliplot(step)로 호출됨,
              % tmp_prop(i) 선형인덱싱이라 87.44%만 실제로 사용됨 (원본 그대로 재현)
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(16), 17, agerangeidx, ...
    [21 33], [87.44 87.44]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'BA (2017)')

nexttile(12)  % PE (2017, Braga 2023) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(13), 17, agerangeidx, ...
    [5 15 25 35 45 55 66], [52.9 89.3 92.9 96.0 96.4 94.8]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'PE (2017)')

nexttile(13)  % CE (2018) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(10), 18, agerangeidx, ...
    [4 20 46 60 82], [33/44 161/177 68/69 111/114], minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'CE (2018)')

nexttile(14)  % CE (2019) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(10), 19, agerangeidx, ...
    [2 5 12], [3.4 32.0; 10.9 43.5]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'CE (2019)')

exportgraphics(gcf,'figures/fig_vali_northeast.jpg','Resolution',600)

%% -- Central-west: 4 panels, 세로 (원본과 동일: 폭 5cm) --
fprintf('Fig: validation Central-west...\n');
figure('Units','centimeters','Position',[0 0 5 16],'PaperPositionMode','auto');
tiledlayout(4,1,'TileSpacing','tight')

nexttile(1)  % GO (2002) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(25), 2, agerangeidx, ...
    [15 50], [37.3]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'GO (2002)')

nexttile(2)  % GO (2009) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(25), 9, agerangeidx, ...
    [25.8-6.4 25.8+6.4], [49.4 49.4; 58.3 58.3]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'GO (2009)')

nexttile(3)  % DF (2020) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(26), 20, agerangeidx, ...
    [7 18 41 61 87], ...
    [40.3 57.9 59.9 71.1; 76.2 65.6 68.2 83.4]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'DF (2020)')

nexttile(4)  % MS (2024) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(23), 24, agerangeidx, ...
    [4 9 13 16], [13.4 22.6 27.7]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'MS (2024)')

exportgraphics(gcf,'figures/fig_vali_centralwest.jpg','Resolution',600)

%% -- Southeast: 6 panels, 2x3 grid --
fprintf('Fig: validation Southeast...\n');
figure('Units','centimeters','Position',[0 0 15 8],'PaperPositionMode','auto');
tiledlayout(2,3,'TileSpacing','tight')

nexttile(1)  % SP (1998) — pre2000m
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(19), 0, agerangeidx, ...
    [0 5 10 15 20 25 30 40 50 80], ...
    [19.63 15.48 11.61 15.93 17.01 15.30 14.62 11.23 15.59 14.80]/100, ...
    minidx, maxidx, col_model, col_model_light, col_m, 'step', 'SP (1998)')

nexttile(2)  % MG (2000) — 2000_agerange (pre-2000 baseline 모델 + multibin 관측치)
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(17), 0, agerangeidx, ...
    [0 16 60 80], [26 18 34.3; 17.02 26.7 18.5]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'MG (2000)')

nexttile(3)  % RJ (2006) — age_multi
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(18), 6, agerangeidx, ...
    [1 5 10 20], ...
    [50 71.4 83.3 20; 13 38.1 75 71; 38.9 41.5 61.8 82.1]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'multibin', 'RJ (2006)')

nexttile(4)  % RJ (2014) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(18), 14, agerangeidx, ...
    [0 17 48 80], [45 96.1 98.2]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'RJ (2014)')

nexttile(5)  % SP (2015) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(19), 15, agerangeidx, ...
    [10 20 40 60 92], [63.1 73.7 74.7 85.9]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'SP (2015)')

nexttile(6)  % SP (2017) — serovaliplot
draw_vali(op_seroposall_pre_0, op_seroposall_0, old_to_new(19), 17, agerangeidx, ...
    [2 5 10 14 17], [33.1 33.1 33.1 33.1]/100, minidx, maxidx, ...
    col_model, col_model_light, col_g, 'step', 'SP (2017)')

exportgraphics(gcf,'figures/fig_vali_southeast.jpg','Resolution',600)