close all

dan_sn_symp = [84  55  41  29   16     2  -10  -28  -37 -45 -51;...
100 100 100 100 100 100  100  100  97 96  93;...
76  33  18  2  -10  -22  -34  -45 -49 -59 -65;...
69  24  10  -2  -16  -29 -40  -47 -55  -59 -67];
dan_sp_symp = [95  74  60   50  42   35  29  23   18  14   12;...
100 91  81   75  68   60  53  45   38  32   25;...
82  48  32   25  19   15   12 10   7   6    5;...
94  76  65   57  50   43   36  30  24  20   15];
dan_sn_hosp = [100  87  83  70  62   49  40  23    6  -10  -23;...
100 100 100 100 100 100   100 100   100 100 100;...
92   62  49  32  15   -6  -27 -45  -70  -91  -100;...
100  70  61  49  36  27   2  -10    -19  -32  -40];
dan_sp_hosp = [97   91  83  78   72   64 56   49   40   33   29;...
100  95  94  90  87   84  79  74    66   58   49;...
95   75  60  51  43   35  30   24  19   16   13;...
98  91  84   80  75   68  62   56  49   43   36];
dan_sn_inf = [0 0 0 0 0 0 0 0 0 0 0;...
66 40 32 26 20 17 13 10 8 6 5;...
0 0 0 0 0 0 0 0 0 0 0;...
0 0 0 0 0 0 0 0 0 0 0];
dan_sp_inf = [48 15 8 7 5 3 2 2 2 1 1;...
71 32 21 16 13 11 7 7 5 5 2;...
20 5 3 1 1 0 0 0 0 0 0;...
53 20 12 10 7 5 5 5 2 1 1];

shen_sn_symp = [56.2, 40.5, 34.8, 31.4, 29.2, 27.5, 26.1, 25.0, 24.1, 23.3, 22.6, 22.0, 21.5, 21.0, 20.6, 20.2, 19.8, 19.4, 19.1, 18.8, 18.5, 18.3, 18.0, 17.8, 17.6, 17.4, 17.2, 17.0, 16.8, 16.6];
shen_sp_symp = [72.0, 53.1, 46.2, 42.1, 39.2, 37.1, 35.5, 34.1, 32.9, 31.9, 31.1, 30.3, 29.6, 29.0, 28.4, 27.9, 27.4, 27.0, 26.6, 26.2, 25.8, 25.5, 25.2, 24.9, 24.6, 24.3, 24.1, 23.8, 23.6, 23.4];
shen_sn_inf = shen_sn_symp/2; shen_sp_inf = shen_sp_symp/2;
shen_hosp = shen_sn_symp*0 + 84.4;

col_dan = 	[0, 62, 116]/255;
col_shen = [225 36 42]/255;

% figure('Units','centimeters',...
%     'Position',[0 0 26 11],...
%     'PaperPositionMode','auto');
% 
% tiledlayout(2,3,"TileSpacing","compact")
% 
% nexttile(1)
% plot(1:11,sum(dan_sn_inf)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_sn_inf,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seronegative')
% legend('Daniels','Shen'); legend('boxoff')
% set(gca,'FontName','times','FontSize',14)
% 
% nexttile(2)
% plot(1:11,sum(dan_sn_symp)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_sn_symp,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seronegative')
% set(gca,'FontName','times','FontSize',14)
% 
% nexttile(3)
% plot(1:11,sum(dan_sn_hosp)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_hosp,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seronegative')
% set(gca,'FontName','times','FontSize',14)
% 
% nexttile(4)
% plot(1:11,sum(dan_sp_inf)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_sp_inf,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seropositive')
% set(gca,'FontName','times','FontSize',14)
% 
% nexttile(5)
% plot(1:11,sum(dan_sp_symp)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_sp_symp,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seropositive')
% set(gca,'FontName','times','FontSize',14)
% 
% nexttile(6)
% plot(1:11,sum(dan_sp_hosp)/4,'^:','LineWidth',1.5,'color',col_dan)
% hold on
% plot(1:30,shen_hosp,'s:','LineWidth',1.5,'color',col_shen)
% xlim([0.5 10.5]); ylim([-20 100])
% yline(0, 'k:')
% % title('Seronegative')
% set(gca,'FontName','times','FontSize',14)


global ohm1 ohm2
ohm1 = 1/1;
ohm2 = 1/4;

t = 1:0.1:10;
mididx = 2;
endidx = 10;

dan_ve_i_sn1 = mean(dan_sn_inf(:,1));   dan_ve_i_sn2 = mean(dan_sn_inf(:,mididx));   dan_ve_i_sn3 = mean(dan_sn_inf(:,endidx));
dan_ve_d_sn1 = mean(dan_sn_symp(:,1));  dan_ve_d_sn2 = mean(dan_sn_symp(:,mididx));   dan_ve_d_sn3 = mean(dan_sn_symp(:,endidx));
dan_ve_h_sn1 = mean(dan_sn_hosp(:,1));  dan_ve_h_sn2 = mean(dan_sn_hosp(:,mididx));   dan_ve_h_sn3 = mean(dan_sn_hosp(:,endidx));
dan_ve_i_sp1 = mean(dan_sp_inf(:,1));  dan_ve_i_sp2 = mean(dan_sp_inf(:,mididx));   dan_ve_i_sp3 = mean(dan_sp_inf(:,endidx));
dan_ve_d_sp1 = mean(dan_sp_symp(:,1));  dan_ve_d_sp2 = mean(dan_sp_symp(:,mididx));   dan_ve_d_sp3 = mean(dan_sp_symp(:,endidx));
dan_ve_h_sp1 = mean(dan_sp_hosp(:,1));  dan_ve_h_sp2 = mean(dan_sp_hosp(:,mididx));   dan_ve_h_sp3 = mean(dan_sp_hosp(:,endidx));
shen_ve_i_sn1 = shen_sn_inf(1); shen_ve_i_sn2 = shen_sn_inf(mididx); shen_ve_i_sn3 = shen_sn_inf(endidx);
shen_ve_d_sn1 = shen_sn_symp(1);  shen_ve_d_sn2 = shen_sn_symp(mididx);   shen_ve_d_sn3 = shen_sn_symp(endidx);  
shen_ve_h_sn1 = shen_hosp(1);  shen_ve_h_sn2 = shen_hosp(mididx);   shen_ve_h_sn3 = shen_hosp(endidx);  
shen_ve_i_sp1 = shen_sp_inf(1); shen_ve_i_sp2 = shen_sp_inf(mididx); shen_ve_i_sp3 = shen_sp_inf(endidx);
shen_ve_d_sp1 = shen_sp_symp(1);  shen_ve_d_sp2 = shen_sp_symp(mididx);   shen_ve_d_sp3 = shen_sp_symp(endidx);  
shen_ve_h_sp1 = shen_hosp(1);  shen_ve_h_sp2 = shen_hosp(mididx);   shen_ve_h_sp3 = shen_hosp(endidx);  

[dan_ve_i_sn1 dan_ve_i_sn2 dan_ve_i_sn3; dan_ve_i_sp1 dan_ve_i_sp2 dan_ve_i_sp3;...
    dan_ve_d_sn1 dan_ve_d_sn2 dan_ve_d_sn3; dan_ve_d_sp1 dan_ve_d_sp2 dan_ve_d_sp3;...
    dan_ve_h_sn1 dan_ve_h_sn2 dan_ve_h_sn3; dan_ve_h_sp1 dan_ve_h_sp2 dan_ve_h_sp3]/100
[shen_ve_i_sn1 shen_ve_i_sn2 shen_ve_i_sn3; shen_ve_i_sp1 shen_ve_i_sp2 shen_ve_i_sp3;...
    shen_ve_d_sn1 shen_ve_d_sn2 shen_ve_d_sn3; shen_ve_d_sp1 shen_ve_d_sp2 shen_ve_d_sp3;...
    shen_ve_h_sn1 shen_ve_h_sn2 shen_ve_h_sn3; shen_ve_h_sp1 shen_ve_h_sp2 shen_ve_h_sp3]/100


[T,Y]=ode45(@vacwan, t, [1, 0 0]);
y1=Y(:,1);
y2=Y(:,2);
y3=Y(:,3);
% veval_i_sn = ve_i_sn1*y1 + ve_i_sn2*y2 + ve_i_sn3*y3 ;
% veval_i_sp = ve_i_sp1*y1 + ve_i_sp2*y2 + ve_i_sp3*y3 ;
% veval_d_sn = ve_d_sn1*y1 + ve_d_sn2*y2 + ve_d_sn3*y3 ;
% veval_d_sp = ve_d_sp1*y1 + ve_d_sp2*y2 + ve_d_sp3*y3 ;
% veval_h_sn = ve_h_sn1*y1 + ve_h_sn2*y2 + ve_h_sn3*y3 ;
% veval_h_sp = ve_h_sp1*y1 + ve_h_sp2*y2 + ve_h_sp3*y3 ;
dan_veval_i_sn = dan_ve_i_sn1*y1 + dan_ve_i_sn2*y2 + dan_ve_i_sn3*y3;
dan_veval_i_sp = dan_ve_i_sp1*y1 + dan_ve_i_sp2*y2 + dan_ve_i_sp3*y3;
dan_veval_d_sn = dan_ve_d_sn1*y1 + dan_ve_d_sn2*y2 + dan_ve_d_sn3*y3;
dan_veval_d_sp = dan_ve_d_sp1*y1 + dan_ve_d_sp2*y2 + dan_ve_d_sp3*y3;
dan_veval_h_sn = dan_ve_h_sn1*y1 + dan_ve_h_sn2*y2 + dan_ve_h_sn3*y3;
dan_veval_h_sp = dan_ve_h_sp1*y1 + dan_ve_h_sp2*y2 + dan_ve_h_sp3*y3;

shen_veval_i_sn = shen_ve_i_sn1*y1 + shen_ve_i_sn2*y2 + shen_ve_i_sn3*y3;
shen_veval_i_sp = shen_ve_i_sp1*y1 + shen_ve_i_sp2*y2 + shen_ve_i_sp3*y3;
shen_veval_d_sn = shen_ve_d_sn1*y1 + shen_ve_d_sn2*y2 + shen_ve_d_sn3*y3;
shen_veval_d_sp = shen_ve_d_sp1*y1 + shen_ve_d_sp2*y2 + shen_ve_d_sp3*y3;
shen_veval_h_sn = shen_ve_h_sn1*y1 + shen_ve_h_sn2*y2 + shen_ve_h_sn3*y3;
shen_veval_h_sp = shen_ve_h_sp1*y1 + shen_ve_h_sp2*y2 + shen_ve_h_sp3*y3;


figure('Units','centimeters',...
    'Position',[0 0 26 13],...
    'PaperPositionMode','auto');
tiledlayout(2,3,"TileSpacing","compact")
nexttile(1)
ppp1 = plot(mean(dan_sn_inf),'^','LineWidth',1.5,'color',col_dan);
hold on
plot(t, dan_veval_i_sn,'-','LineWidth',1.5,'color',col_dan)
ppp2 = plot(shen_sn_inf,'s','LineWidth',1.5,'color',col_shen);
plot(t, shen_veval_i_sn,'-','LineWidth',1.5,'color',col_shen)
ppp3 = plot(nan, nan,'-','LineWidth',1.5,'color','k')
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
legend([ppp1, ppp2, ppp3],'Daniels','Shen','Model'); legend('boxoff')
set(gca,'FontName','SansSerif','FontSize',13)

nexttile(4)
plot(mean(dan_sp_inf),'^','LineWidth',1.5,'color',col_dan)
hold on
plot(t, dan_veval_i_sp,'-','LineWidth',1.5,'color',col_dan)
plot(shen_sp_inf,'s','LineWidth',1.5,'color',col_shen)
plot(t, shen_veval_i_sp,'-','LineWidth',1.5,'color',col_shen)
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
set(gca,'FontName','SansSerif','FontSize',13)

nexttile(2)
plot(mean(dan_sn_symp),'^','LineWidth',1.5,'color',col_dan)
hold on
plot(t, dan_veval_d_sn,'-','LineWidth',1.5,'color',col_dan)
plot(shen_sn_symp,'s','LineWidth',1.5,'color',col_shen)
plot(t, shen_veval_d_sn,'-','LineWidth',1.5,'color',col_shen)
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
set(gca,'FontName','SansSerif','FontSize',13)

nexttile(5)
plot(mean(dan_sp_symp),'^','LineWidth',1.5,'color',col_dan)
hold on
plot(t, dan_veval_d_sp,'-','LineWidth',1.5,'color',col_dan)
plot(shen_sp_symp,'s','LineWidth',1.5,'color',col_shen)
plot(t, shen_veval_d_sp,'-','LineWidth',1.5,'color',col_shen)
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
set(gca,'FontName','SansSerif','FontSize',13)


nexttile(3)
plot(mean(dan_sn_hosp),'^','LineWidth',1.5,'color',col_dan)
hold on
plot(t, dan_veval_h_sn,'-','LineWidth',1.5,'color',col_dan)
plot(shen_hosp,'s','LineWidth',1.5,'color',col_shen)
plot(t, shen_veval_h_sn,'-','LineWidth',1.5,'color',col_shen)
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
set(gca,'FontName','SansSerif','FontSize',13)

nexttile(6)
plot(mean(dan_sp_hosp),'^','LineWidth',1.5,'color',col_dan)
hold on
plot(t, dan_veval_h_sp,'-','LineWidth',1.5,'color',col_dan)
plot(shen_hosp,'s','LineWidth',1.5,'color',col_shen)
plot(t, shen_veval_h_sp,'-','LineWidth',1.5,'color',col_shen)
xlim([0.5 10.5]); ylim([-20 100])
yline(0, 'k:')
set(gca,'FontName','SansSerif','FontSize',13)

% 
% figure
% plot(mean(dan_sn_symp),'rs')
% hold on
% plot(t, veval_d_sn,'k-')
% 
% figure
% plot(mean(dan_sn_symp),'rs')
% hold on
% plot(t, veval_d_sn,'k-')

% figure
% plot(mean(dan_sp_inf),'rs')
% hold on
% plot(t, veval_i_sp,'k-')
% 
% figure
% plot(mean(dan_sp_symp),'rs')
% hold on
% plot(t, veval_d_sp,'k-')
% 
% figure
% plot(mean(dan_sp_symp),'rs')
% hold on
% plot(t, veval_d_sp,'k-')

function odeout = vacwan(t,y)
global ohm1 ohm2
dy1 = -ohm1*y(1);
dy2 = ohm1*y(1) - ohm2*y(2);
dy3 = ohm2*y(2);
odeout = [dy1; dy2; dy3];
end

