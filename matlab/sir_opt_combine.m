%% ═══════════════════════════════════════════════════════════════
%  Brazil DENV — Optimal allocation results: combined figure
%  6개 mat 파일 (opt_vac1_R0{1,2,3}.mat, opt_vac6_R0{1,2,3}.mat) 로드
%% ═══════════════════════════════════════════════════════════════
close all; clear; clc;

%% ── 0. 설정 ───────────────────────────────────────────────────
S = 27;
state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
               'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
               'RO','RR','RS','SC','SE','SP','TO'};

north_order     = [21,1,3,22,14,4,27];
northeast_order = [10,17,6,20,15,16,2,25,5];
centwest_order  = [12,13,9,7];
southeast_order = [11,8,19,26];
south_order     = [18,24,23];
reg_order  = [north_order, northeast_order, centwest_order, southeast_order, south_order];
codes_reord = state_codes(reg_order);
reg_bound  = [0, cumsum([7,9,4,4,3])];
reg_labels = {'North','Northeast','Central-west','Southeast','South'};

col_region = [1.000 0.498 0.055;
              0.596 0.306 0.639;
              0.090 0.745 0.812;
              0.173 0.627 0.173;
              0.549 0.341 0.271];

col_child   = [0.122 0.467 0.706];   % blue  (5-14yr)
col_elderly = [0.839 0.153 0.157];   % red   (60+yr)
col_brazil  = [0.3 0.3 0.3];         % gray  (Brazil actual)
col_equal   = [0.6 0.6 0.6];         % light gray (equal)

r0_labels  = {'q25','q50','q75'};
vac_labels = {'5-14yr','60+yr'};
vac_orders = [1, 6];
fav = 0.4;

mkdir('figures')

%% ── 1. 결과 로드 ─────────────────────────────────────────────
% results struct: [2 vac × 3 R0]
res = cell(2,3);
for vi = 1:2
    for rc = 1:3
        fname = sprintf('results/opt_vac%d_R0%d.mat', vac_orders(vi), rc);
        if exist(fname,'file')
            res{vi,rc} = load(fname);
            fprintf('Loaded: %s\n', fname)
        else
            warning('Missing: %s', fname)
        end
    end
end

%% ── 2. Fig A: 최적 배분 비율 패널맵 (2×3) ───────────────────
figure('Units','centimeters','Position',[0 0 22 12],'PaperPositionMode','auto');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact')

for vi = 1:2
    for rc = 1:3
        nexttile((vi-1)*3+rc)
        r = res{vi,rc};
        if isempty(r); continue; end

        opt_r = r.xopt_ratio(reg_order)*100;
        bra_r = r.brazil_ratio(reg_order)*100;
        ymax  = max([opt_r; bra_r])*1.35;

        % 지역 배경
        for rr = 1:5
            fill([reg_bound(rr) reg_bound(rr+1) reg_bound(rr+1) reg_bound(rr)]+0.5, ...
                 [0 0 ymax ymax], col_region(rr,:),'EdgeColor','none','FaceAlpha',fav)
            hold on
        end

        col = col_child*(vi==1) + col_elderly*(vi==2);
        bar(opt_r, 0.65,'FaceColor',col,'EdgeColor','none','DisplayName','Optimal')
        hold on
        plot(bra_r,'o-','Color',col_brazil,'MarkerSize',3,'LineWidth',1, ...
            'MarkerFaceColor',col_brazil,'DisplayName','Brazil actual')

        xlim([0.5 S+0.5]); ylim([0 ymax])
        xticks(1:S); xticklabels(codes_reord)
        set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',6.5)
        ylabel('Allocation (%)','FontSize',7)
        title(sprintf('%s / R0-%s  (−%.1f%%)', ...
            vac_labels{vi}, r0_labels{rc}, ...
            (r.daly_brazil-r.daly_opt)/r.daly_brazil*100), ...
            'FontSize',8,'FontWeight','bold')

        if vi==1 && rc==1
            legend('Location','northeast','FontSize',7,'Box','off')
        end
    end
end

exportgraphics(gcf,'figures/fig_opt_allocation_all.jpg','Resolution',600)
fprintf('Saved: fig_opt_allocation_all.jpg\n')

%% ── 3. Fig B: q50 기준 최적 vs 브라질 vs 균등 비교 (2패널) ──
figure('Units','centimeters','Position',[0 0 16 6],'PaperPositionMode','auto');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

for vi = 1:2
    nexttile(vi)
    r = res{vi,2};   % q50
    if isempty(r); continue; end

    opt_r = r.xopt_ratio(reg_order)*100;
    bra_r = r.brazil_ratio(reg_order)*100;
    equ_r = ones(S,1)/S*100;
    ymax  = max([opt_r; bra_r])*1.3;

    for rr = 1:5
        fill([reg_bound(rr) reg_bound(rr+1) reg_bound(rr+1) reg_bound(rr)]+0.5, ...
             [0 0 ymax ymax], col_region(rr,:),'EdgeColor','none','FaceAlpha',fav)
        hold on
    end

    col = col_child*(vi==1) + col_elderly*(vi==2);
    bar(opt_r, 0.65,'FaceColor',col,'EdgeColor','none','DisplayName','Optimal')
    hold on
    plot(bra_r,'o-','Color',col_brazil,'MarkerSize',3,'LineWidth',1.2, ...
        'MarkerFaceColor',col_brazil,'DisplayName','Brazil actual')
    plot(equ_r,'s--','Color',col_equal,'MarkerSize',2,'LineWidth',1, ...
        'DisplayName','Equal')

    xlim([0.5 S+0.5]); ylim([0 ymax])
    xticks(1:S); xticklabels(codes_reord)
    set(gca,'XTickLabelRotation',90,'FontName','SansSerif','FontSize',7)
    ylabel('Allocation (%)','FontSize',8)
    title(sprintf('%s (q50)  Optimal −%.1f%% vs Brazil', ...
        vac_labels{vi}, (r.daly_brazil-r.daly_opt)/r.daly_brazil*100), ...
        'FontSize',9,'FontWeight','bold')
    legend('Location','northeast','FontSize',8,'Box','off')
end

exportgraphics(gcf,'figures/fig_opt_allocation_q50.jpg','Resolution',600)
fprintf('Saved: fig_opt_allocation_q50.jpg\n')

%% ── 4. Fig C: DALY 절감률 (Brazil 대비) ─────────────────────
% [2 vac × 3 R0] 행렬
reduction_mat = zeros(2,3);
for vi = 1:2
    for rc = 1:3
        r = res{vi,rc};
        if ~isempty(r)
            reduction_mat(vi,rc) = (r.daly_brazil-r.daly_opt)/r.daly_brazil*100;
        end
    end
end

figure('Units','centimeters','Position',[0 0 10 6],'PaperPositionMode','auto');
x_pos = [1 2 3];
b1 = bar(x_pos-0.2, reduction_mat(1,:), 0.35, ...
    'FaceColor',col_child,'EdgeColor','none'); hold on
b2 = bar(x_pos+0.2, reduction_mat(2,:), 0.35, ...
    'FaceColor',col_elderly,'EdgeColor','none');

for rc = 1:3
    text(rc-0.2, reduction_mat(1,rc)+0.05, sprintf('%.1f%%',reduction_mat(1,rc)), ...
        'HorizontalAlignment','center','FontSize',7,'FontName','SansSerif')
    text(rc+0.2, reduction_mat(2,rc)+0.05, sprintf('%.1f%%',reduction_mat(2,rc)), ...
        'HorizontalAlignment','center','FontSize',7,'FontName','SansSerif')
end

xticks(1:3); xticklabels({'R0-q25','R0-q50','R0-q75'})
ylabel('DALY reduction vs Brazil allocation (%)')
ylim([0 max(reduction_mat(:))*1.3])
legend([b1,b2], vac_labels,'Location','northeast','FontSize',9,'Box','off')
set(gca,'FontName','SansSerif','FontSize',9)
exportgraphics(gcf,'figures/fig_opt_daly_reduction.jpg','Resolution',600)
fprintf('Saved: fig_opt_daly_reduction.jpg\n')

%% ── 5. Fig D: 상관관계 — 브라질 배분 vs 최적 배분 (q50) ─────
figure('Units','centimeters','Position',[0 0 12 6],'PaperPositionMode','auto');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact')

for vi = 1:2
    nexttile(vi)
    r = res{vi,2};
    if isempty(r); continue; end

    col = col_child*(vi==1) + col_elderly*(vi==2);
    bra = r.brazil_ratio*100;
    opt = r.xopt_ratio*100;

    % 지역별 색으로 scatter
    region_of = zeros(S,1);
    for rr=1:5
        idx_list = {north_order,northeast_order,centwest_order,southeast_order,south_order};
        region_of(idx_list{rr}) = rr;
    end
    for rr = 1:5
        idx = region_of==rr;
        scatter(bra(idx), opt(idx), 30, col_region(rr,:), 'filled','DisplayName',reg_labels{rr})
        hold on
    end

    % 대각선 (브라질=최적 기준선)
    lim_val = max([bra;opt])*1.05;
    plot([0 lim_val],[0 lim_val],'k--','LineWidth',1,'HandleVisibility','off')

    % 회귀선
    p = polyfit(bra, opt, 1);
    xfit = linspace(0, lim_val, 100);
    plot(xfit, polyval(p,xfit), '-','Color',col,'LineWidth',1.5,'HandleVisibility','off')

    [rho,pval] = corr(bra,opt);
    text(0.05,0.95,sprintf('r = %.2f\np = %.3f',rho,pval), ...
        'Units','normalized','VerticalAlignment','top', ...
        'FontName','SansSerif','FontSize',8)

    xlabel('Brazil actual allocation (%)','FontSize',8)
    ylabel('Optimal allocation (%)','FontSize',8)
    title(sprintf('%s (q50)',vac_labels{vi}),'FontSize',9,'FontWeight','bold')
    xlim([0 lim_val]); ylim([0 lim_val]); axis square
    if vi==1
        legend('Location','southeast','FontSize',7,'Box','off')
    end
    set(gca,'FontName','SansSerif','FontSize',8)
end

exportgraphics(gcf,'figures/fig_opt_scatter.jpg','Resolution',600)
fprintf('Saved: fig_opt_scatter.jpg\n')

fprintf('\nAll figures saved to figures/\n')