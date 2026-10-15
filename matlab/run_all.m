function run_all(run_optimization)
%RUN_ALL Run the MATLAB analyses in order.
%   run_all       reuses the optimized allocations in results/
%   run_all(true) also reruns the optimization (M09, slow)
%
%   M01-M05 need results_big/ (posterior draws); M06-M12 only need the
%   small files in results/ created by M05 and M09.

if nargin < 1; run_optimization = false; end
here = fileparts(mfilename('fullpath'));
steps = {'M01_fig1_surveillance', 'M02_catalytic_figures', 'M03_catalytic_tables', ...
         'M04_serosurvey_validation', 'M05_r0_projection', 'M06_sir_no_vaccine', ...
         'M07_vaccine_effectiveness', 'M08_baseline_vaccination', 'M09_optimize_allocation', ...
         'M10_optimal_allocation_figures', 'M11_mixed_strategy', 'M12_table1'};
for k = 1:numel(steps)
    if strcmp(steps{k}, 'M09_optimize_allocation') && ~run_optimization
        fprintf('\nSkipping %s (call run_all(true) to rerun it)\n', steps{k});
        continue
    end
    fprintf('\n===== %s =====\n', steps{k});
    run_step(fullfile(here, [steps{k} '.m']));
end
end

function run_step(script)
% separate workspace, because every script starts with clear
run(script);
end
