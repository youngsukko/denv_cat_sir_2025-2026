%% Fig. S11 and Table S9: vaccine effectiveness from Daniels and Shen and the
%  waning curve implied by the three vaccine-derived immunity states.
%
%  Output: output/figures/figS11_vaccine_effectiveness.png
%          output/tables/tableS9_vaccine_effectiveness.csv

clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));
cfg = project_config();
src = ve_source_data();

col_dan  = [0, 62, 116] / 255;
col_shen = [225, 36, 42] / 255;

%% Occupancy of the <1 y, 1-5 y and >5 y states after vaccination at t = 0
ohm1 = 1; ohm2 = 1/4;   % as in the SIR-type model
t = (1:0.1:10)';
[~, Y] = ode45(@(~, y) [-ohm1*y(1); ohm1*y(1) - ohm2*y(2); ohm2*y(2)], t, [1 0 0]);

ve_dan  = ve_settings('daniels');
ve_shen = ve_settings('shen');

% {outcome field in ve_settings, Daniels table (neg, pos), Shen series (neg, pos)}
outcomes = {'inf', src.daniels.neg_inf,  src.daniels.pos_inf,  src.shen.neg_inf,  src.shen.pos_inf
            'dis', src.daniels.neg_symp, src.daniels.pos_symp, src.shen.neg_symp, src.shen.pos_symp
            'hsp', src.daniels.neg_hosp, src.daniels.pos_hosp, src.shen.hosp,     src.shen.hosp};

%% Fig. S11: columns infection / disease / hospitalization, rows seronegative / seropositive
figure('Units', 'centimeters', 'Position', [0 0 26 13], 'PaperPositionMode', 'auto');
tiledlayout(2, 3, 'TileSpacing', 'compact')
for j = 1:3
    for sero = 1:2
        nexttile((sero - 1) * 3 + j)
        dan_markers  = mean(outcomes{j, 1 + sero}, 1);
        shen_markers = outcomes{j, 3 + sero};
        dan_curve    = Y * ve_dan.(outcomes{j, 1})(:, sero) * 100;
        shen_curve   = Y * ve_shen.(outcomes{j, 1})(:, sero) * 100;

        p1 = plot(dan_markers, '^', 'LineWidth', 1.5, 'Color', col_dan); hold on
        plot(t, dan_curve, '-', 'LineWidth', 1.5, 'Color', col_dan)
        p2 = plot(shen_markers, 's', 'LineWidth', 1.5, 'Color', col_shen);
        plot(t, shen_curve, '-', 'LineWidth', 1.5, 'Color', col_shen)
        yline(0, 'k:')
        xlim([0.5 10.5]); ylim([-20 100])
        set(gca, 'FontName', 'SansSerif', 'FontSize', 13)
        if j == 1 && sero == 1
            p3 = plot(nan, nan, '-', 'LineWidth', 1.5, 'Color', 'k');
            legend([p1 p2 p3], {'Daniels', 'Shen', 'Model'}); legend('boxoff')
        end
    end
end
save_figure(cfg, 'figS11_vaccine_effectiveness')

%% Table S9: VE by immunity state and serostatus
rows  = {'Infection', 'Disease', 'Hospitalization'};
state = {'lt1y', 'lt5y', 'gt5y'};
sero  = {'neg', 'pos'};
vars  = {};
for st = 1:3
    for sg = 1:2
        vars{end+1} = sprintf('%s_%s', state{st}, sero{sg}); %#ok<SAGROW>
    end
end
ves  = {ve_shen, ve_dan};
tbls = cell(1, 2);
for k = 1:2
    ve   = ves{k};
    vals = [reshape(ve.inf', 1, []); reshape(ve.dis', 1, []); reshape(ve.hsp', 1, [])];
    tbls{k} = array2table(round(vals, 2), 'VariableNames', strcat(ve.label, '_', vars), 'RowNames', rows);
end
out = [tbls{:}];
writetable(out, fullfile(cfg.table_dir, 'tableS9_vaccine_effectiveness.csv'), 'WriteRowNames', true)
disp(out)
