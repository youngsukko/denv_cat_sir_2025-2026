%% Brazil region color map
close all; clc;

col_region = [1.000 0.498 0.055;   % North
              0.596 0.306 0.639;   % Northeast
              0.090 0.745 0.812;   % Central-west
              0.173 0.627 0.173;   % Southeast
              0.549 0.341 0.271];  % South

ibge_north     = [11,12,13,14,15,16,17];
ibge_northeast = [21,22,23,24,25,26,27,28,29];
ibge_centwest  = [50,51,52,53];
ibge_southeast = [31,32,33,35];
ibge_south     = [41,42,43];

label_offset = containers.Map( ...
    {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
     'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
     'RO','RR','RS','SC','SE','SP','TO'}, ...
    {[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0.8],[0,0],[0,0],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,-0.8],[0,0],[0,0],[0,-0.8],[0,0], ...
     [0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]});

Shfile = shaperead('BR_UF_2023.shp');
cxvec = zeros(1,length(Shfile));
cyvec = zeros(1,length(Shfile));
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    valid = ~isnan(x) & ~isnan(y);
    pg = polyshape(x(valid),y(valid));
    [cx,cy] = centroid(pg);
    cxvec(i) = cx; cyvec(i) = cy;
end
cxvec(26) = cxvec(26) - 1;

figure('Units','centimeters','Position',[0 0 10 12],'PaperPositionMode','auto');
hold on
for i = 1:length(Shfile)
    x = Shfile(i).X; y = Shfile(i).Y;
    cd_uf = str2double(Shfile(i).CD_UF);
    if     ismember(cd_uf, ibge_north);     fc = col_region(1,:);
    elseif ismember(cd_uf, ibge_northeast); fc = col_region(2,:);
    elseif ismember(cd_uf, ibge_centwest);  fc = col_region(3,:);
    elseif ismember(cd_uf, ibge_southeast); fc = col_region(4,:);
    elseif ismember(cd_uf, ibge_south);     fc = col_region(5,:);
    else;                                    fc = [0.8 0.8 0.8];
    end
    nan_idx = find(isnan(x));
    starts = [1, nan_idx+1]; ends = [nan_idx-1, length(x)];
    for j = 1:length(starts)
        seg = starts(j):ends(j);
        if length(seg) >= 3
            patch(x(seg),y(seg),fc,'EdgeColor','w','LineWidth',0.8,'FaceAlpha',1)
        end
    end
end

for i = 1:length(Shfile)
    off = label_offset(Shfile(i).SIGLA_UF);
    text(cxvec(i)+off(1), cyvec(i)+off(2), Shfile(i).SIGLA_UF, ...
        'HorizontalAlignment','center','FontName','SansSerif')
end
axis equal off
set(gca,'FontName','SansSerif','FontSize',14)
exportgraphics(gcf,'figures/fig_map_region.jpg','Resolution',600)