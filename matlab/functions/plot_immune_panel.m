function plot_immune_panel(Xk, ages, cmap, panel_title)
%PLOT_IMMUNE_PANEL Stacked bars of the proportion with 0..4 prior infections.
%   Xk: 1x5 cell of [n_draw x A5] draws; the black curve marks the
%   proportion with zero or one prior infection.

denom = Xk{1} + Xk{2} + Xk{3} + Xk{4} + Xk{5};
prop  = zeros(5, numel(ages));
for k = 1:5
    prop(k, :) = mean(Xk{k} ./ denom, 1);
end
bb = bar(ages, prop, 1.1, 'stacked', 'EdgeColor', 'none');
for k = 1:5; bb(k).FaceColor = cmap(k, :); end
hold on
plot(ages, prop(1, :) + prop(2, :), 'k-', 'LineWidth', 1.5)
xlim([0 80]); ylim([0 1]); xticks([0 20 40 60 80])
set(gca, 'FontName', 'SansSerif', 'FontSize', 8)
title(panel_title, 'FontSize', 11, 'FontWeight', 'normal')
end
