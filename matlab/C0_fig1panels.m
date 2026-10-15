close all; clear; clc;

%% ── 0. 기본 설정 ──────────────────────────────────────────────
S    = 27; T = 26; A5 = 17;
years = 2000:2025;

state_codes_alph = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
                    'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
                    'RO','RR','RS','SC','SE','SP','TO'};
state_names_alph = {'Acre','Alagoas','Amazonas','Amapá','Bahia','Ceará', ...
                    'Distrito Federal','Espírito Santo','Goiás','Maranhão', ...
                    'Minas Gerais','Mato Grosso do Sul','Mato Grosso','Pará', ...
                    'Paraíba','Pernambuco','Piauí','Paraná','Rio de Janeiro', ...
                    'Rio Grande do Norte','Rondônia','Roraima','Rio Grande do Sul', ...
                    'Santa Catarina','Sergipe','São Paulo','Tocantins'};

ibge_codes = [12,27,13,16,29,23,53,32,52,21, ...
              31,50,51,15,25,26,22,41,33,24, ...
              11,14,43,42,28,35,17];

age_labels = {'0-','5-','10-','15-','20-','25-','30-','35-', ...
              '40-','45-','50-','55-','60-','65-','70-','75-','80+'};

% 지역 인덱스 (알파벳순 기준)
north_idx     = find(ismember(state_codes_alph,{'AC','AM','AP','PA','RO','RR','TO'}));
northeast_idx = find(ismember(state_codes_alph,{'AL','BA','CE','MA','PB','PE','PI','RN','SE'}));
centwest_idx  = find(ismember(state_codes_alph,{'DF','GO','MS','MT'}));
southeast_idx = find(ismember(state_codes_alph,{'ES','MG','RJ','SP'}));
south_idx     = find(ismember(state_codes_alph,{'PR','RS','SC'}));

col_region = [1.000 0.498 0.055;   % North     - orange
              0.596 0.306 0.639;   % Northeast - purple
              0.090 0.745 0.812;   % Central-west - teal
              0.173 0.627 0.173;   % Southeast - green
              0.549 0.341 0.271];  % South     - brown

region_of = zeros(S,1);
region_of(north_idx)=1; region_of(northeast_idx)=2; region_of(centwest_idx)=3;
region_of(southeast_idx)=4; region_of(south_idx)=5;

%% ── 1. CSV 로드 + 연령 컬럼 순서 버그 보정 ────────────────────
wrong_order = {'0-4','10-14','15-19','20-24','25-29','30-34','35-39', ...
               '40-44','45-49','5-9','50-54','55-59','60-64', ...
               '65-69','70-74','75-79','80+'};
correct_order_full = {'0-4','5-9','10-14','15-19','20-24','25-29','30-34', ...
                       '35-39','40-44','45-49','50-54','55-59','60-64', ...
                       '65-69','70-74','75-79','80+'};
[~, age_reorder_idx] = ismember(correct_order_full, wrong_order);
assert(all(age_reorder_idx > 0), '연령 컬럼 매핑 실패');

cases_mat = table2array(readtable('results/cases_age_state_year.csv'));
pop_mat   = table2array(readtable('results/pop_age_state_year.csv'));

% 연령 순서 보정
cases_mat = cases_mat(:, age_reorder_idx);
pop_mat   = pop_mat(:, age_reorder_idx);

% 주별 연도별 총 케이스/인구 [T × S]
data_case_bystate = zeros(T, S);
data_pop_bystate  = zeros(T, S);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_case_bystate(:,s) = sum(cases_mat(idx,:), 2);
    data_pop_bystate(:,s)  = sum(pop_mat(idx,:),   2);
end

% 전국 연도×연령 케이스 [T × A5]
data_case_nationwide = zeros(T, A5);
for s = 1:S
    idx = (s-1)*T + (1:T);
    data_case_nationwide = data_case_nationwide + cases_mat(idx,:);
end

%% ── 2. 인구 대비 케이스 (주별 총합) ──────────────────────────
data_pop_ref  = sum(data_pop_bystate, 1)';   % [S × 1]
cases_per_pop = sum(data_case_bystate)' ./ data_pop_ref;

figure('Units','centimeters','Position',[0 0 8 5],'PaperPositionMode','auto');
bar(cases_per_pop * 100)
xticks(1:S); xticklabels(state_codes_alph)
ylabel('Cases per 100 population (cumulative)')
set(gca,'FontName','SansSerif','FontSize',8,'XTickLabelRotation',90)
%% ── 3. 브라질 choropleth 지도 (누적 인시던스 per 100) ──────────
Shfile = shaperead('BR_UF_2023.shp');
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    valid = ~isnan(x) & ~isnan(y);
    pg = polyshape(x(valid), y(valid));
    [cx,cy] = centroid(pg);
    cxvec(i)=cx; cyvec(i)=cy;
end

inc_vals = cases_per_pop * 100;   % [S × 1]
inc_norm = inc_vals / max(inc_vals);

figure('Units','centimeters','Position',[0 0 10 12],'PaperPositionMode','auto');
hold on
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    cd_uf = str2double(Shfile(i).CD_UF);
    si = find(ibge_codes == cd_uf);
    isnan_idx = isnan(x);
    parts = find(isnan_idx);
    start_idx = [1, parts+1]; end_idx = [parts-1, length(x)];
    for j = 1:length(start_idx)
        idx = start_idx(j):end_idx(j);
        if ~isempty(si)
            patch(x(idx), y(idx), 'r', 'EdgeColor','k', 'LineWidth',0.8, ...
                  'FaceAlpha', inc_norm(si));
        else
            patch(x(idx), y(idx), [0.7 0.7 0.7], 'EdgeColor','k', 'LineWidth',0.8);
        end
    end
end
axis equal off
set(gca,'FontName','SansSerif','FontSize',11)
exportgraphics(gcf,'figures/fig_choropleth_incidence.jpg','Resolution',600)


%% ── 3. 브라질 지역 색상 지도 (첨부 이미지 스타일) ──────────────
Shfile = shaperead('BR_UF_2023.shp');

% 각 주의 centroid 계산
cxvec = zeros(1, length(Shfile));
cyvec = zeros(1, length(Shfile));
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    valid = ~isnan(x) & ~isnan(y);
    pg = polyshape(x(valid), y(valid));
    [cx, cy] = centroid(pg);
    cxvec(i) = cx; cyvec(i) = cy;
end

figure('Units','centimeters','Position',[0 0 10 12],'PaperPositionMode','auto');
hold on
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    cd_uf = str2double(Shfile(i).CD_UF);
    si = find(ibge_codes == cd_uf);

    if ~isempty(si)
        fc = col_region(region_of(si), :);
    else
        fc = [0.7 0.7 0.7];
    end

    % NaN으로 분리된 멀티폴리곤 처리
    nan_idx = find(isnan(x));
    starts  = [1, nan_idx+1];
    ends    = [nan_idx-1, length(x)];
    for j = 1:length(starts)
        seg = starts(j):ends(j);
        if length(seg) >= 3
            patch(x(seg), y(seg), fc, 'EdgeColor','k', 'LineWidth',0.8, 'FaceAlpha',1)
        end
    end

    % % 주 코드 레이블
    % if ~isempty(si)
    %     text(cxvec(i), cyvec(i), state_codes_alph{si}, ...
    %         'HorizontalAlignment','center','FontSize',7.5, ...
    %         'FontName','SansSerif','FontWeight','bold','Color','k')
    % end

label_offset = containers.Map( ...
    {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
     'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
     'RO','RR','RS','SC','SE','SP','TO'}, ...
    {[0,0.5],[0,0],[0,0],[0,0],[0,0],[0,0],[0,1.2],[0,0],[0,0],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,-1.2],[0,0],[0,0],[0.3,-1.5],[0,0], ...
     [0,0],[0,0],[0,0],[1.4,0],[0,0],[0,0],[0,0]});

% 레이블 출력 시 보정 적용
if ~isempty(si)
    off = label_offset(state_codes_alph{si});
    text(cxvec(i)+off(1), cyvec(i)+off(2), state_codes_alph{si}, ...
        'HorizontalAlignment','center','FontSize',7.5, ...
        'FontName','SansSerif','FontWeight','normal','Color','k')
end

end
axis equal off
exportgraphics(gcf,'figures/fig_map_region.jpg','Resolution',600)

%% ── 4. 연도별 연령 분포 (색 = 연도) ──────────────────────────
figure('Units','centimeters','Position',[0 0 9 6],'PaperPositionMode','auto');
for t = 1:T
    prop = data_case_nationwide(t,:) / sum(data_case_nationwide(t,:));
    plot(prop, 'Color', [(t-1)/(T-1), 1-(t-1)/(T-1), 1], 'LineWidth', 1)
    hold on
end
xticks(1:A5); xticklabels(age_labels)
ylabel('Proportion by age group')
xlabel('Age group')
xlim([1 A5])
pp1 = plot(nan, nan, '-', 'Color',[0 1 1], 'LineWidth',1.5);
pp2 = plot(nan, nan, '-', 'Color',[1 0 1], 'LineWidth',1.5);
% legend([pp1, pp2],{'2000','2025'},'Location','northeast','FontSize',8); legend('boxoff')
set(gca,'FontName','SansSerif','FontSize',11)
exportgraphics(gcf,'figures/fig_age_dist_byyear.jpg','Resolution',600)