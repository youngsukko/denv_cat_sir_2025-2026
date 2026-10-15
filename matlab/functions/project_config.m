function cfg = project_config()
%PROJECT_CONFIG Paths, dimensions, state metadata and scenario settings.
%   All MATLAB scripts call this first. Inputs produced by the R pipeline
%   are read from the project root (data/processed, results_big), so the
%   MATLAB folder never holds copies of them.

fn_dir         = fileparts(mfilename('fullpath'));
cfg.matlab_dir = fileparts(fn_dir);
cfg.root_dir   = fileparts(cfg.matlab_dir);

%% Inputs from the R pipeline
cfg.file_cases = fullfile(cfg.root_dir, 'data', 'processed', 'dengue_cases_year_age5_state_2000_2025.csv');
cfg.file_pop   = fullfile(cfg.root_dir, 'data', 'processed', 'pop_state_age5_2000_2025.csv');
cfg.file_ppc   = fullfile(cfg.root_dir, 'results_big', 'catalytic_ppc_draws.mat');
cfg.file_param = fullfile(cfg.root_dir, 'results_big', 'catalytic_param_draws.mat');

%% Small results shared between MATLAB scripts (tracked by git)
cfg.results_dir     = fullfile(cfg.root_dir, 'results');
cfg.file_sir_inputs = fullfile(cfg.results_dir, 'sir_inputs.mat');
cfg.file_opt        = @(ve_name) fullfile(cfg.results_dir, sprintf('opt_allocation_%s.mat', ve_name));

%% MATLAB input data
cfg.data_dir   = fullfile(cfg.matlab_dir, 'data');
cfg.file_shape = fullfile(cfg.data_dir, 'shapefile', 'BR_UF_2023.shp');
cfg.file_doses = fullfile(cfg.data_dir, 'vaccine_doses_by_state_2024_2025.csv');

%% Outputs (not tracked by git)
cfg.fig_dir   = fullfile(cfg.matlab_dir, 'output', 'figures');
cfg.table_dir = fullfile(cfg.matlab_dir, 'output', 'tables');
cfg.sim_dir   = fullfile(cfg.matlab_dir, 'output', 'sim');
out_dirs = {cfg.results_dir, cfg.fig_dir, cfg.table_dir, cfg.sim_dir};
for k = 1:numel(out_dirs)
    if ~isfolder(out_dirs{k}); mkdir(out_dirs{k}); end
end

%% Dimensions
cfg.years     = 2000:2025;
cfg.S         = 27;
cfg.T         = numel(cfg.years);
cfg.A5        = 17;
cfg.age_lower = 0:5:80;
cfg.age_labels      = [arrayfun(@(a) sprintf('%d-', a), 0:5:75, 'UniformOutput', false), {'80+'}];
cfg.age_labels_full = [arrayfun(@(a) sprintf('%d-%d', a, a+4), 0:5:75, 'UniformOutput', false), {'80+'}];

%% States, alphabetical by abbreviation (the order used by the R pipeline)
cfg.state_codes = {'AC','AL','AM','AP','BA','CE','DF','ES','GO','MA', ...
                   'MG','MS','MT','PA','PB','PE','PI','PR','RJ','RN', ...
                   'RO','RR','RS','SC','SE','SP','TO'};
cfg.state_names = {'Acre','Alagoas','Amazonas','Amapá','Bahia','Ceará', ...
                   'Distrito Federal','Espírito Santo','Goiás','Maranhão', ...
                   'Minas Gerais','Mato Grosso do Sul','Mato Grosso','Pará', ...
                   'Paraíba','Pernambuco','Piauí','Paraná','Rio de Janeiro', ...
                   'Rio Grande do Norte','Rondônia','Roraima','Rio Grande do Sul', ...
                   'Santa Catarina','Sergipe','São Paulo','Tocantins'};
cfg.ibge_codes  = [12,27,13,16,29,23,53,32,52,21, ...
                   31,50,51,15,25,26,22,41,33,24, ...
                   11,14,43,42,28,35,17];

%% Geographic regions
cfg.region_names  = {'North','Northeast','Central-west','Southeast','South'};
cfg.region_short  = {'N','NE','CW','SE','S'};
cfg.region_colors = [1.000 0.498 0.055;   % North
                     0.596 0.306 0.639;   % Northeast
                     0.090 0.745 0.812;   % Central-west
                     0.173 0.627 0.173;   % Southeast
                     0.549 0.341 0.271];  % South
region_members = {{'AC','AM','AP','PA','RO','RR','TO'}, ...
                  {'AL','BA','CE','MA','PB','PE','PI','RN','SE'}, ...
                  {'DF','GO','MS','MT'}, ...
                  {'ES','MG','RJ','SP'}, ...
                  {'PR','RS','SC'}};
cfg.region_of = zeros(cfg.S, 1);
for r = 1:5
    cfg.region_of(ismember(cfg.state_codes, region_members{r})) = r;
end

% State order on the x-axis of state-level figures (grouped by region)
order_codes = {'RO','AC','AM','RR','PA','AP','TO', ...
               'MA','PI','CE','RN','PB','PE','AL','SE','BA', ...
               'MS','MT','GO','DF', ...
               'MG','ES','RJ','SP', ...
               'PR','SC','RS'};
cfg.reg_order = cellfun(@(c) find(strcmp(cfg.state_codes, c)), order_codes);
cfg.reg_bound = [0, cumsum(arrayfun(@(r) sum(cfg.region_of == r), 1:5))];

%% Tile layout of the state panel maps (approximate geography of Brazil)
cfg.tile_rows   = 8;
cfg.tile_cols   = 6;
cfg.legend_tile = 5;
tile_def = {'RR',1,2; 'AP',1,3; 'AM',2,2; 'PA',2,3; 'MA',2,4; 'CE',2,5;
            'AC',3,1; 'RO',3,2; 'TO',3,3; 'PI',3,4; 'PB',3,5; 'RN',3,6;
            'MT',4,2; 'GO',4,3; 'BA',4,4; 'PE',4,5; 'AL',4,6;
            'MS',5,2; 'DF',5,3; 'MG',5,4; 'SE',5,5;
            'SP',6,3; 'RJ',6,4; 'ES',6,5;
            'PR',7,3; 'SC',7,4; 'RS',8,3};
cfg.state_tile = zeros(cfg.S, 1);
for k = 1:size(tile_def, 1)
    si = strcmp(cfg.state_codes, tile_def{k,1});
    cfg.state_tile(si) = (tile_def{k,2} - 1) * cfg.tile_cols + tile_def{k,3};
end

%% Catalytic-to-SIR conversion
cfg.life_expectancy = 76.4;   % years; mu* = 1 / life expectancy (endemic equilibrium)
cfg.R0_quantiles    = [0.25 0.50 0.75];
cfg.R0_names        = {'q25','q50','q75'};
cfg.n_proj_years    = 10;     % years projected with Gaussian process regression

%% Vaccination scenarios
cfg.n_sim_years   = 5;        % simulation horizon, 2026-2030
cfg.annual_doses  = 5e6;      % limited supply, individuals vaccinated per year
cfg.mass_doses    = [10e6 20e6 30e6];
cfg.target_ages.children = [2 3];     % 5-14 years
cfg.target_ages.elderly  = 13:17;     % 60+ years
cfg.target_ages.mass     = 2:12;      % 5-59 years
cfg.vacorder.children    = 1;         % the same age groups as selected in func_ode_alt_v5
cfg.vacorder.elderly     = 6;
cfg.vacorder.mass        = 8;
cfg.ve_names      = {'shen','daniels'};
cfg.ve_labels     = {'Shen','Daniels'};
cfg.target_names  = {'Children (5-14)','Elderly (60+)'};
end
