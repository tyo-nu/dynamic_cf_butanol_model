% RUN_DYNAMIC_EM Set up and run Ensemble Modeling/Optimization for dynamic
% kinetic model

% Revision history:
%{
2020-09-30: jpm
    Script created
2022-02-24
    Publication version started

%}

function output_struct = run_dynamic_EM(...
    load_file, varargin)

% Set up input parser
isName = @(x) isstring(x) || ischar(x) || iscell(x);
isBool = @(x) isnumeric(x) || islogical(x);
isTwoColumnNum = @(x) isnumeric(x) && size(x, 2) == 2;

p = inputParser;

% Set defaults
% These are absolute value ranges
default_k_range = [1e0, 1e3]; % Rate constants
default_unknown_Ki_range = [1e-1, 1e1]; % Inhibition constants
default_unknown_sat_constant_range = [0.01, 1]; % Set range for saturation constants (will be used for unknowns if sat_constant_bounds are also given)
default_min_abs_kcat = 1e-2; 

default_event_fxn_timeout = 5; % seconds

% Get inputs
addRequired(p, 'load_file', isName)
addParameter(p, 'best_ks_loadfile', '', isName)
addParameter(p, 'n_init_ksets', [], @isnumeric)
addParameter(p, 'UseParallel', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'k_range', default_k_range, @(x) isnumeric(x) && ~isempty(x));
addParameter(p, 'savename', [], isName);
addParameter(p, 'best_n_ksets_as_init', [], @isnumeric)
addParameter(p, 'final_timepoint_weight', [], @(x) (0 <= x) && (x <= 1))
addParameter(p, 'use_ED_for_init_concs', true, isBool);
addParameter(p, 'default_init_conc', 0, @isnumeric);
addParameter(p, 'exp_conditions', [], @isnumeric)
addParameter(p, 'spec_init_met_concs', [], @(x) isnumeric(x) || iscell(x))
addParameter(p, 'fitness_fxn', [], @(x) isName(x) || iscell(x))
addParameter(p, 'sat_constant_bounds', [], @isnumeric) % Fold-change bounds from lit values within which to sampling
addParameter(p, 'rate_constant_bounds', [], @isnumeric) % Fold-change bounds from lit values within which to sampling
addParameter(p, 'Options_struct_loadfile', [], isName)
addParameter(p, 'unknown_sat_constant_range', default_unknown_sat_constant_range, ... % Set range for saturation constants (will be used for unknowns if sat_constant_bounds are also given)
    @(x) isnumeric(x) && ~isempty(x))
addParameter(p, 'opt_alg', 'EM', isName)
addParameter(p, 'xover_ratio', [], @(x)isnumeric(x) && x >= 1)
addParameter(p, 'rate_const_distribution', 'log', isName)
addParameter(p, 'sat_const_distribution', 'log', isName)
addParameter(p, 'use_event_fxn', true, isBool)
addParameter(p, 'event_fxn_timeout', [], @isnumeric)
addParameter(p, 'best_n_ksets', [], @isnumeric) % Screen in EM to best n
addParameter(p, 'fitness_cutoff', [], @isnumeric)
addParameter(p, 'manual_k_values', [], @isnumeric) % Set individual k's to a specific value (no range)
addParameter(p, 'manual_k_centers', [], @isnumeric) % Set new centers for sampling (use existing fc change "ranges")
addParameter(p, 'manual_k_ranges', [], @isnumeric) % Set absolute ranges for specific parameters
addParameter(p, 'dataset_for_sampling', 1, @isnumeric) % TODO - clean this up
addParameter(p, 'rng_setting', 'default', isName)
addParameter(p, 'ode_nonneg', false, isBool)
addParameter(p, 'const_metabs', [], @iscell) % TODO - figure this out and do setup fxns
addParameter(p, 'write_mass_bal_elast_fxns', false, isBool) % TODO - figure this out and do setup fxns
addParameter(p, 'use_spec_rng_seeds', [], @isnumeric)
addParameter(p, 'use_reg_rxns', true, isBool) % ONLY USE IN SETUP
addParameter(p, 'unknown_Ki_range', default_unknown_Ki_range, @isnumeric)
addParameter(p, 'min_abs_kcat', default_min_abs_kcat, @isnumeric);
addParameter(p, 'metabs_for_fitness_fxn', [], isName)
addParameter(p, 'save_EM_population', true, isBool); % Save entire initial ensemble
addParameter(p, 'optim_rate_const_bounds', [], @isnumeric) % Set bounds for optimization ONLY
addParameter(p, 'optim_sat_const_bounds', [], @isnumeric) % Set bounds for optimization ONLY
addParameter(p, 'optim_inhib_const_bounds', [], @isnumeric) % Set bounds for optimization ONLY
addParameter(p, 'make_bounds_include_best_ksets', false, isBool) % Regardless of final_xxx_const_bounds, make sure that bounds include the min/max values from best n ksets, if that was used
addParameter(p, 'metab_fitness_multipliers', {}, @iscell)
addParameter(p, 'MaxTime', Inf, @isnumeric) % Timeout for optimization ONLY
addParameter(p, 'MaxFunEvals', [], @isnumeric) % Max evals for optimization ONLY (for patternsearch, will be called many times per iteration)
addParameter(p, 'init_rng_seed', 0, @isnumeric) % Increase all rng seeds for multiple parallel runs
addParameter(p, 'keep_all_fits', true, isBool) % Whether to save fitness values for whole ensemble or only X best
addParameter(p, 'MaxIterations', [], @isnumeric) % Max iterations in patternsearch
addParameter(p, 'use_mex', true, isBool); % Option to use compiled mex mass balance function - gives minimal gain over pre-written function
addParameter(p, 'PenaltyFactor', [], @(x)isnumeric(x) && x > 1) % For GA
addParameter(p, 'InitialPenalty', [], @(x)isnumeric(x) && x >= 1)  % For GA
addParameter(p, 'NonlinConAlgorithm', 'auglag', isName) % Algorithm to apply nonlinear (thermo) constraints to GA
addParameter(p, 'keep_only_fits', false, isBool) % Option to save ONLY fitness values (no results from integration) - useful when memory is very limited
addParameter(p, 'use_nonlincon', true, isBool) % Option to use thermodynamic constraints
addParameter(p, 'use_log_params', false, isBool) % log transform param values
addParameter(p, 'best_init_kset_inds', [], @isnumeric)
addParameter(p, 'optim_param_bounds', [], @isnumeric) % Set new bounds on params during optimization - WILL ALLOW EXCLUDING INITIAL POINTS
addParameter(p, 'optim_init_params', [], @isnumeric) % Set specific intial parameters to new values during optimization
addParameter(p, 'initial_concentration_perturbations', {}, @iscell)
addParameter(p, 'ignore_unused_params', false, isBool)
addParameter(p, 'param_inds_for_optimization', [], @isnumeric) % Option to only otpimize certain parameters
addParameter(p, 'condition_fitness_weights', [], isTwoColumnNum)
addParameter(p, 'AbsTol', [], @isnumeric)
addParameter(p, 'RelTol', [], @isnumeric)
addParameter(p, 'run_fake_optimization', false, isBool) % For testing - bypass actual patternsearch to test rest of workflow

p.KeepUnmatched = false;
p.CaseSensitive = false;

% Allow inputing all varargin in a single cell (you probably shouldn't do this)
if length(varargin) == 1
    varargin = varargin{1};
end

parse(p, load_file, varargin{:});

best_ks_loadfile = p.Results.best_ks_loadfile;
n_init_ksets = p.Results.n_init_ksets;
UseParallel = p.Results.UseParallel;
base_k_range = p.Results.k_range;
savename = p.Results.savename;
best_n_ksets_as_init = p.Results.best_n_ksets_as_init;
final_timepoints_weight = p.Results.final_timepoint_weight;
fitness_fxn = p.Results.fitness_fxn;
default_init_conc = p.Results.default_init_conc;
exp_conditions = p.Results.exp_conditions;
spec_init_met_concs = p.Results.spec_init_met_concs;
sat_constant_bounds = p.Results.sat_constant_bounds;
rate_constant_bounds = p.Results.rate_constant_bounds; 
Options_struct_loadfile = p.Results.Options_struct_loadfile;
unknown_sat_constant_range = p.Results.unknown_sat_constant_range; 
opt_alg = p.Results.opt_alg;
xover_ratio = p.Results.xover_ratio;
rate_const_distribution = p.Results.rate_const_distribution;
sat_const_distribution = p.Results.sat_const_distribution;
event_fxn_timeout = p.Results.event_fxn_timeout;
use_event_fxn = p.Results.use_event_fxn;
best_n_ksets = p.Results.best_n_ksets;
fitness_cutoff = p.Results.fitness_cutoff;
dataset_for_sampling = p.Results.dataset_for_sampling;
rng_setting = p.Results.rng_setting;
unknown_Ki_range = p.Results.unknown_Ki_range;
save_EM_population = p.Results.save_EM_population;
keep_all_fits = p.Results.keep_all_fits;
keep_only_fits = p.Results.keep_only_fits;
best_init_kset_inds = p.Results.best_init_kset_inds;
run_fake_optimization = p.Results.run_fake_optimization;

% Load model and get needed fields
loadstruct = load(load_file);
model = loadstruct.model;
if ~isempty(Options_struct_loadfile)
    Options_struct = load(Options_struct_loadfile);
    Options = Options_struct.Options;
elseif isfield(loadstruct, 'Options')
    Options = loadstruct.Options;
else
    Options = loadstruct.default_Options;
end
Experimental_Data = loadstruct.Experimental_Data;
if isfield(loadstruct, 'Initial_Ensemble')
    Initial_Ensemble = loadstruct.Initial_Ensemble;
end
if isfield(loadstruct, 'Perturbation_Results')
    Perturbation_Results = loadstruct.Perturbation_Results;
end
if isfield(loadstruct, 'EM_population_trimmed')
    EM_population_trimmed = loadstruct.EM_population_trimmed;
elseif isfield(loadstruct, 'EM_population')
    EM_population = loadstruct.EM_population;    
end
if isfield(loadstruct, 'EM_PR_combined')
    EM_PR_combined = loadstruct.EM_PR_combined;  
end

if isempty(fitness_fxn)
    Options.fitness_fxn = 'RMSE';
else
    error("TODO\n")
end

% Set rng
rng(rng_setting)

Options.rate_const_distribution = rate_const_distribution;
Options.sat_const_distribution = sat_const_distribution;
Options.use_event_fxn = use_event_fxn;
Options.manual_k_values = p.Results.manual_k_values;
Options.manual_k_centers = p.Results.manual_k_centers;
Options.manual_k_ranges = p.Results.manual_k_ranges;
Options.ode_nonneg = p.Results.ode_nonneg;
Options.const_metabs = p.Results.const_metabs;
Options.write_mass_bal_elast_fxns = p.Results.write_mass_bal_elast_fxns;
Options.use_spec_rng_seeds = p.Results.use_spec_rng_seeds;
Options.use_reg_rxns = p.Results.use_reg_rxns;
Options.min_abs_kcat = p.Results.min_abs_kcat;
Options.metabs_for_fitness_fxn = p.Results.metabs_for_fitness_fxn;
Options.metab_fitness_multipliers = p.Results.metab_fitness_multipliers;
Options.MaxTime = p.Results.MaxTime;
Options.MaxFunEvals = p.Results.MaxFunEvals;
Options.init_rng_seed = p.Results.init_rng_seed;
Options.MaxIterations = p.Results.MaxIterations;
Options.use_mex = p.Results.use_mex;
Options.PenaltyFactor = p.Results.PenaltyFactor;
Options.InitialPenalty = p.Results.InitialPenalty;
Options.NonlinConAlgorithm = p.Results.NonlinConAlgorithm;
Options.use_nonlincon = p.Results.use_nonlincon;
Options.use_log_params = p.Results.use_log_params;
Options.optim_rate_const_bounds = p.Results.optim_rate_const_bounds;
Options.optim_sat_const_bounds = p.Results.optim_sat_const_bounds;
Options.optim_inhib_const_bounds = p.Results.optim_inhib_const_bounds;
Options.optim_param_bounds = p.Results.optim_param_bounds;
Options.optim_init_params = p.Results.optim_init_params;
Options.initial_concentration_perturbations = p.Results.initial_concentration_perturbations;
Options.use_ED_for_init_concs = p.Results.use_ED_for_init_concs;
Options.ignore_unused_params = p.Results.ignore_unused_params;
Options.param_inds_for_optimization = p.Results.param_inds_for_optimization;
Options.condition_fitness_weights = p.Results.condition_fitness_weights;
Options.AbsTol = p.Results.AbsTol;
Options.RelTol = p.Results.RelTol;

%% Set Options if varargin are given
% If specific EF timeout given, turn on EF - doing it here so that giving
% an EF time in varargin automatically turns on EF
if ~isempty(event_fxn_timeout)
    Options.event_fxn_timeout = event_fxn_timeout;
    Options.use_event_fxn = true;
else
    % Otherwise use default value but don't auto turn on
    Options.event_fxn_timeout = default_event_fxn_timeout;
end

% Get model properties
[n_metabs, n_rxns] = size(model.S);
[n_metabs_and_fracs, n_ks] = size(model.S_f_b);

% If optioned, trim down to only some experimental conditions
if ~isempty(exp_conditions)
    [Experimental_Data, Options] = ...
        select_experimental_data_conditions(Experimental_Data, Options,...
        exp_conditions);
end

if ~isempty(Options.initial_concentration_perturbations)
    Experimental_Data = add_perturbations_to_ED(...
        model, Experimental_Data, Options);
end

% Set new options from input parser into default_Options
if ~isempty(final_timepoints_weight)
    Options.fitness_fxn = 'wgt_timecourse_rmse';
    Options.final_timepoint_weight = final_timepoints_weight;
elseif ~isempty(fitness_fxn)
    Options.fitness_fxn = fitness_fxn;
end

% Get reference condition from default_Options
if isfield(Options, 'ref_from_ED_experimental_condition')
    ref_condition = Options.ref_from_ED_experimental_condition;
else
    ref_condition = 1;
end

metab_inds_and_concs = ...
    Experimental_Data.initial_metab_concs{ref_condition};
init_ED_inds = metab_inds_and_concs(:, 1);
if iscell(metab_inds_and_concs)
    init_ED_inds = cellfun(@(x)find(ismember(model.mets, x)), ...
        init_ED_inds);
end
init_ED_concs = metab_inds_and_concs(:, 2);
if iscell(init_ED_concs)
    init_ED_concs = cell2mat(init_ED_concs);
end
% Get from timecourse if available
timecourse = Experimental_Data.metab_timecourse{ref_condition};
if ~isempty(timecourse) && timecourse{1,2}(1,1) == 0
    for m = size(timecourse, 1)
        met_ind = find(strcmpi(model.mets, timecourse{m,1}));
        init_conc = timecourse{m,2}(2,1);
        
        prev_met_ind = find(init_ED_inds == met_ind);
        if isempty(prev_met_ind)
            init_ED_inds = [init_ED_inds; met_ind];
            init_ED_concs = [init_ED_concs; init_conc];
        else
            init_ED_concs(prev_met_ind) = init_conc;
        end
    end
end

init_metab_bounds = ones(n_metabs, 2) .* default_init_conc;
init_metab_bounds(init_ED_inds, :) = repmat(init_ED_concs, 1, 2);

% If optioned, overwrite any metabs given in Experimental_Data
if Options.use_ED_for_init_concs && ~isempty(timecourse)
    init_t_ind = find(timecourse{1, 2}(1, :) == 0);
    assert(~isempty(init_t_ind), ...
        "Experimental_Data.metab_timecourse does not have t=0 measurement\n")
        
    % Just need a dummy value to pass into this function
    second_timepoint = Experimental_Data.metab_timecourse{ref_condition}...
        {1, 2}(1, 2);
    
    [~, init_tc_concs, init_tc_inds] = get_metab_dxdt_from_timecourse_data(...
        model, Experimental_Data, Options,...
        'ref_condition', ref_condition,...
        'ref_timeframe', [0, second_timepoint]);
    
    init_metab_bounds(init_tc_inds, :) = repmat(init_tc_concs, 1, 2);
    
end

% If specific concs were give, use those as final option
if ~isempty(spec_init_met_concs)
    assert(size(spec_init_met_concs, 2) == 2, ...
        "Input `spec_init_met_concs` must have two columns\n")
    init_met_inds = spec_init_met_concs(:, 1);
    init_met_concs = spec_init_met_concs(:,2);
    if iscell(spec_init_met_concs)
        [~, init_met_inds] = ismember(init_met_inds, model.mets);
        init_met_concs = cell2mat(init_met_concs);
    end
    init_metab_bounds(init_met_inds, :) = repmat(init_met_concs, 1, 2);
    
end

% Loop through lit datasets - if only one, make into cells so that we
% can loop anyway
if ~iscell(model.lit_values.rate_constants)
    model.lit_values.rate_constants = {model.lit_values.rate_constants};
    model.lit_values.sat_constants = {model.lit_values.sat_constants};
else
    assert(length(model.lit_values.sat_constants) == ...
        length(model.lit_values.rate_constants), "Fix this")
end

% Loop over each lit dataset
n_lit_datasets = length(model.lit_values.rate_constants);
% Make sure `dataset_for_sampling` is in range
assert(dataset_for_sampling <= n_lit_datasets,...
    "Field `dataset_for_sampling` is greater than # of datasets")

ds = dataset_for_sampling;
        
lit_kcats = model.lit_values.rate_constants{ds};

% If new ranges from varargin given to override specific rate
% constants from model.lit_values field, use those
if ~isempty(Options.manual_k_centers)
    assert(size(Options.manual_k_centers, 2) == 2, ...
        ["Field `manual_k_centers` must be nx2 numeric,",...
        "first column of indices, second of values"])

    lit_kcats(Options.manual_k_centers(:, 1)) = ...
        Options.manual_k_centers(:, 2);

end

% If optioned, set Ki (inhibition constant) ranges
reg_rxns = transpose(find(ismember(model.rxn_type, [6,7])));
reg_elem_rxns = [];
for rr = reg_rxns
    reg_elem_rxns = [reg_elem_rxns, ...
        model.rxn_indices(rr,1) : model.rxn_indices(rr,2) ];
end
if Options.use_reg_rxns
    Options.reg_elem_rxns = reg_elem_rxns;
    geomean_Ki = geomean(unknown_Ki_range);
    % Only replace for rxns with no existing value
    lit_kcats(reg_elem_rxns(isnan(lit_kcats(reg_elem_rxns)))) = ...
        0;
    lit_kcats(reg_elem_rxns(lit_kcats(reg_elem_rxns) == 0)) = ...
        geomean_Ki;

else
    lit_kcats(reg_elem_rxns) = 0;
end

% If general (multiplier) bounds are given, use those around the lit
% values
if ~isempty(rate_constant_bounds)
    all_kcat_ranges = lit_kcats .* rate_constant_bounds;
else
    % Otherwise use general ranges
    all_kcat_ranges = repmat(base_k_range, n_ks, 1);
end

% Get info for sat constants
% Get indices of elem rxns that are type 5 or 50 (multiplicative or
% additive saturation)
any_type_sat_net_rxns = find(ismember(model.rxn_type, [5, 50]));
any_type_sat_net_rxns = any_type_sat_net_rxns(:)';
any_type_sat_elem_rxns = [];
for sat_net = any_type_sat_net_rxns
    any_type_sat_elem_rxns = [any_type_sat_elem_rxns, ...
        model.rxn_indices(sat_net, 1) : model.rxn_indices(sat_net, 2)];
end

if any(model.rxn_type == 4)
    error("TODO: implement this")
end

%Get as a vector
sat_const_matrix = model.lit_values.sat_constants{ds};
all_sat_const_ranges = [];
% Save matrix locations
sat_constants_matrix_locs = [];
for sat_elem_rxn = any_type_sat_elem_rxns

    % See if cell array or not
    if isnumeric(sat_const_matrix)
        subs_inds = find(~isnan(sat_const_matrix(sat_elem_rxn,:)));
    else
        subs_inds = find(cellfun(@(x)~isempty(x),sat_const_matrix(sat_elem_rxn,:)));
    end
    % Get substrates (or subs & products) in stoich
    net_rxn_ind = intersect(find(model.rxn_indices(:,1) <= sat_elem_rxn), ...
        find(model.rxn_indices(:,2) >= sat_elem_rxn) );
    if model.rxn_type(net_rxn_ind) == 50
        alt_subs_inds = find(model.S_f_b(:, sat_elem_rxn) < 0);
    else
        alt_subs_inds = find(model.S_f_b(:, sat_elem_rxn) < 0);
    end

    if ~isequal(sort(subs_inds(:)), sort(alt_subs_inds(:)))
        fprintf("Warning: lit_values.sat_constants has %s as substrates,\nand S_f_b has %s\n",...
            strjoin(model.mets(subs_inds), ", "),...
            strjoin(model.mets(alt_subs_inds),", "))
    end

    excluded_species = ["h_c","h_e"];
    excluded_inds = find(ismember(model.mets, excluded_species));
    subs_inds = setdiff(subs_inds, excluded_inds);

    for s = subs_inds
        % If numeric, just get value, apply range if needed
        if isnumeric(sat_const_matrix)
            elem_rxn_const = sat_const_matrix(sat_elem_rxn, s);
            % Add min and max to make ranges within the loop
            elem_rxn_sat_min = elem_rxn_const;
            elem_rxn_sat_max = elem_rxn_const;
        else
            % If cell, see if multiple values - apply ranges to min and
            % max
            elem_rxn_cell_val = sat_const_matrix{sat_elem_rxn, s};
            if isnumeric(elem_rxn_cell_val)
                % If single value, just keep
                cell_vals = elem_rxn_cell_val;
            else
                % If multiple values, will be saved as a cell of chars
                cell_vals = cellfun(@str2num, split(elem_rxn_cell_val) );
            end

            % Get min and max for range
            elem_rxn_sat_min = min(cell_vals);
            elem_rxn_sat_max = max(cell_vals);

        end

        if isempty(sat_constant_bounds)
            elem_rxn_sat_const_bounds = [elem_rxn_sat_min, elem_rxn_sat_max];
        else
            elem_rxn_sat_const_bounds = ...
                [sat_constant_bounds(1)*elem_rxn_sat_min, ...
                 sat_constant_bounds(2)*elem_rxn_sat_max];
        end
        all_sat_const_ranges = [all_sat_const_ranges; elem_rxn_sat_const_bounds];

        sat_constants_matrix_locs = [sat_constants_matrix_locs; ...
            [sat_elem_rxn, s] ];
    end
end

% If any sat constants are zero (unknown), replace ranges with that
%  used for unknown constants
% A row *should* have zeros in both columns or not at all
min_zeros = find(all_sat_const_ranges(:,1) == 0);
max_zeros = find(all_sat_const_ranges(:,2) == 0);
missing_cols = setxor(min_zeros, max_zeros);
if ~isempty(missing_cols)
    error("Fix this")
end
missing_rows = union(min_zeros, max_zeros);

all_sat_const_ranges(missing_rows, :) = repmat(unknown_sat_constant_range,...
    length(missing_rows), 1);

n_sat_consts = size(sat_constants_matrix_locs, 1);

% If no fc bounds were given, overwrite with single default range (for unknowns)
if isempty(sat_constant_bounds)
    fprintf("Warning: overriding lit values for saturation constants and using default value\n")
    all_sat_const_ranges = repmat(unknown_sat_constant_range, n_sat_consts, 1);
end

% Save info into model struct 
model.lit_values.sat_constants_matrix_locs = sat_constants_matrix_locs;

% Adjust specific ranges
if ~isempty(Options.manual_k_ranges)
    assert(size(Options.manual_k_ranges, 2) == 3, ...
        ["Field `manual_k_ranges` must be nx3 numeric,",...
        "first column of indices, then min/max values"])

    n_rate_consts = size(all_kcat_ranges, 1);

    for ind = 1:size(Options.manual_k_ranges, 1)

        if Options.manual_k_ranges(ind, 1) <= n_rate_consts

            all_kcat_ranges(Options.manual_k_ranges(ind, 1), :) = ...
                Options.manual_k_ranges(ind, 2:3);
        else
            all_sat_const_ranges(Options.manual_k_ranges(ind, 1) - n_rate_consts, :) = ...
                Options.manual_k_ranges(ind, 2:3);
        end
    end
end

% Add on saturation constants
all_k_ranges = [all_kcat_ranges; all_sat_const_ranges];
all_bounds = all_k_ranges;
% Save initial concentrations in model struct
model.init_concs = geomean(init_metab_bounds, 2);

lb = all_bounds(:, 1);
ub = all_bounds(:, 2);
n_vars = length(lb);

model.kcat_range = base_k_range;
        
% Add specific ranges for each kcat/Km
model.all_kcat_ranges = all_kcat_ranges;
model.Km_range = unknown_sat_constant_range;
model.all_Km_ranges = all_sat_const_ranges;

% Pull out important info
n_rate_consts = size(all_kcat_ranges, 1);
rate_const_inds = 1:n_rate_consts;
n_sat_consts = size(all_sat_const_ranges, 1);
sat_const_inds = (length(lit_kcats) + 1) : (length(lit_kcats) + n_sat_consts);

% If optioned, get k's and initial fractions of best n ksets to use as
% subset of initial swarm population
if ~isempty(best_n_ksets_as_init) || ~isempty(best_ks_loadfile)

    % If separate file given for best n ksets, load that
    if ~isempty(best_ks_loadfile)
        best_n_ksets_loadstruct = load(best_ks_loadfile);
    
        if isfield(best_n_ksets_loadstruct,'EM_PR_combined')
            PR = best_n_ksets_loadstruct.EM_PR_combined;
        elseif isfield(best_n_ksets_loadstruct,'Perturbation_Results')
            PR = best_n_ksets_loadstruct.Perturbation_Results;
        else
            error("Load file needs Perturbation_Results struct to set initial params")
        end

        % Prioritize `EM_population` field
        if isfield(best_n_ksets_loadstruct,'EM_population') && ...
                ~isempty(best_n_ksets_loadstruct.EM_population)
            all_ks = best_n_ksets_loadstruct.EM_population;
        elseif isfield(best_n_ksets_loadstruct,'EM_population_trimmed')
            all_ks = best_n_ksets_loadstruct.EM_population_trimmed;
        elseif isfield(best_n_ksets_loadstruct,'Initial_Ensemble')
            all_ks = best_n_ksets_loadstruct.Initial_Ensemble.K_sets;
        else
            error("Load file needs Initial_Ensemble/EM_populationto set initial params")
        end
    else
        % Get PR and IE, whatever theyre called
        if exist('EM_PR_combined','var')
            PR = EM_PR_combined;
        elseif exist('Perturbation_Results','var')
            PR = Perturbation_Results;
        else
            error("Load file needs Perturbation_Results struct to set initial params")
        end

        if exist('EM_population_trimmed','var')
            all_ks = EM_population_trimmed;
        elseif exist('EM_population','var')
            all_ks = EM_population;
        elseif exist('Initial_Ensemble','var')
            all_ks = Initial_Ensemble.K_sets;
        else
            error("Load file needs Initial_Ensemble/EM_populationto set initial params")
        end
    end
    
    % If only `best_ks_loadfile` was given and not `best_n_ksets_as_init`,
    % use all
    if isempty(best_n_ksets_as_init)
        % Get init_ksets_inds so that this varargin can also be used with
        % `best_init_kset_inds`
        init_ksets_inds = 1 : size(all_ks, 2);
    else
        % Otherwise, take best n ksets
        [~, ~, init_ksets_fits, init_ksets_inds] = find_best_Ksets(PR, best_n_ksets_as_init);
    end
    
    if ~isempty(best_init_kset_inds)
        init_ksets_fits = init_ksets_fits(best_init_kset_inds);
        init_ksets_inds = init_ksets_inds(best_init_kset_inds);
    end

    % Get reference condition
    if isfield(loadstruct, 'Options') &&...
            isfield(loadstruct.Options, 'ref_from_ED_experimental_condition')
        ref_cond = loadstruct.Options.ref_from_ED_experimental_condition;
    else
        ref_cond = Options.ref_from_ED_experimental_condition;
    end

    init_ks = all_ks(:, init_ksets_inds);
    get_cells_init_concs = @(x) x(1, :)';
    init_all_concs_cells = cellfun(get_cells_init_concs, ...
        PR.last_metab_concs(ref_cond, init_ksets_inds),...
        'UniformOutput', false);
    init_all_concs = cell2mat(init_all_concs_cells);
    
    init_all_vars = init_ks;
    
    % For now, overwrite anything else
    InitialSwarmMatrix = init_all_vars;

end
        
rng default

Options.adjust_only_ED_enz_levels = false; % TESTING - should be false

%% Set up options used across any solver
% Add constraints on minimum rate constants
% Set lb to be non-zero
min_rate_const = 1e-3;
lb(isnan(lb)) = 0;
ub(isnan(ub)) = 0;
lb(lb < min_rate_const) = min_rate_const;
ub(ub < lb) = lb(ub < lb) .* 1.1;

% Save things into Options/model
Options.rate_const_inds = rate_const_inds;
Options.sat_const_inds = sat_const_inds;
Options.reg_elem_rxns = reg_elem_rxns;
Options.lb = lb;
Options.ub = ub;

% Set up constraints
if ~strcmpi(opt_alg, 'EM')
    
    Options.opt_alg = opt_alg;
    optimization_struct = setup_optimization(model, Experimental_Data,...
        Options, InitialSwarmMatrix);
    % Pull stuff out
    Options = optimization_struct.Options;
    InitialSwarmMatrix = optimization_struct.params;
    Options.lb = optimization_struct.lb;
    Options.ub = optimization_struct.ub;
    A_eq = optimization_struct.A;
    b_eq = optimization_struct.b;
    nonlcon = optimization_struct.nonlcon;
    n_vars = size(InitialSwarmMatrix, 1);
    
    % Set up fields needed for optimization - output function
    optim_output_history.x = [];
    optim_output_history.fval = [];
    optim_output_history.meshsize = [];
    optim_output_history.searchdir = [];
end

% Set objective function AFTER - needs to have updated Options
obj_fxn = @(vars) calc_rates_and_fitness(model, Experimental_Data, vars, Options);

%% Which solver to use
switch opt_alg
    
% Ensemble modeling, either for initial ensemble generation/screening 
%  or for testing perturbatinos
case 'EM' 
    
    if isempty(n_init_ksets)
        error("Must give a population size")
    end
    n_ksets = n_init_ksets;
    % Make sure it agrees if other args were given
    if ~isempty(best_init_kset_inds)
        n_ksets = length(best_init_kset_inds);
    elseif ~isempty(best_n_ksets_as_init)
        n_ksets = best_n_ksets_as_init;
    end
        
    % Start overall timer
    starttime = tic;
    
    % Set up options to generate population in parallel
    EM_creation_fxn_options.PopulationSize = 1;
    
    if ~isempty(fitness_cutoff)
        passing_population = cell(1, n_ksets);
    end
    
    % TODO - make this general to any opt method
    if Options.write_mass_bal_elast_fxns
        
        single_creation_options.PopulationSize = 1;

        temp_params = create_param_set_population(n_vars, [], single_creation_options, ...
            model, Options, []);
   
        write_compiled_mass_bal_elast_fxns(model, Experimental_Data, Options,...
            temp_params );

    end
    
    % Always save fitnesses   
    EM_fits = NaN(n_ksets, 1);
    % Only save others if not trimming to best n ksets
    if isempty(best_n_ksets)
        EM_PR_all = cell(n_ksets, 1);
        Population = zeros(n_vars, n_ksets);
    else
        % Placeholder to not enter loop later when consolidating into 1 struct
        EM_PR_all = {};
    end
    
    % Allocate any field used within parfor, regardless of options - parfor
    % doesn't know what's going to be used
    kset_runtimes = zeros(n_ksets, 1);
    all_start_time = tic();
    % Allocate fields that don't exist, even if not used with any set of
    % options
    if ~exist('InitialSwarmMatrix', 'var')
        InitialSwarmMatrix = zeros(0,n_ksets);
    end
    
%     for ind = 1:n_ksets % Switch commenting to debug within parfor loop
    parfor ind = 1:n_ksets
        
        % Check if using pre-made params
        if ~isempty(best_n_ksets_as_init)
            ith_param_set = InitialSwarmMatrix(:, ind);
        else
        	% Otherwise generate params in parallel
            ith_param_set = create_param_set_population(n_vars, [], EM_creation_fxn_options, ...
                model, Options, ind);
            ith_param_set = ith_param_set(:);
        end
        
        % Start timer
        kset_start = tic();
        
        % Get kset run
        [EM_fits(ind), ith_PR] = calc_rates_and_fitness(model, Experimental_Data, ...
            ith_param_set, Options);
        
        % Stop timer
        kset_time = toc(kset_start);
        if ~keep_only_fits
            kset_runtimes(ind) = kset_time;
        end
        
        % If saving all (not saving best n ksets), add to cell array
        if isempty(best_n_ksets) 
            EM_PR_all{ind} = ith_PR;
            Population(:, ind) = ith_param_set;
        end  
       
        if mod(ind, floor(n_ksets / 100)) == 0
            fprintf("Finished %i ksets - total time: %0.0f seconds\n",...
                ind, toc(all_start_time) );
        end
        
    end
    
    % Get single best kset's rng seed
    [~, best_kset_ind] = min(EM_fits);
    if ~isempty(Options.use_spec_rng_seeds)
        best_kset_rng_seed = Options.use_spec_rng_seeds(best_kset_ind);
    else
        best_kset_rng_seed = Options.init_rng_seed + best_kset_ind;
    end
    % And then get optimum param set
    opt_params = create_param_set_population(n_vars, [], EM_creation_fxn_options, ...
            model, Options, best_kset_rng_seed);
        
    % If trimming to best n ksets, do that
    if ~isempty(best_n_ksets)
        [~, sorted_ksets] = sort(EM_fits, 'ascend');
        best_kset_inds = sorted_ksets(1 : best_n_ksets);
        
        % Save old fits
        EM_fits_all = EM_fits;
        
        % Make new fields to save just best n ksets
        EM_PR_all = cell(best_n_ksets, 1);
        EM_fits = zeros(best_n_ksets, 1);
        Population = zeros(n_vars, best_n_ksets);
        parfor k = 1:best_n_ksets
            % Make kset with that rng
            kth_best_param_set = create_param_set_population(n_vars, [], EM_creation_fxn_options, ...
                model, Options, best_kset_inds(k));
            % Run kth best kset
            [EM_fits(k), EM_PR_all{k}] = calc_rates_and_fitness(model, Experimental_Data, ...
                kth_best_param_set, Options);
            Population(:, k) = kth_best_param_set;
            
            if mod(k, floor(n_ksets / 100)) == 0
                fprintf("Finished best %i ksets - total time: %0.0f seconds\n",...
                    k, toc(all_start_time) );
            end
        end

    end
    
    exitflag = [];
    output = [];
    options_struct = [];
    
    % Combine into single PR struct
    passed_inds = find(cellfun(@(x) ~isempty(x), (EM_PR_all)));
    EM_PR_combined = struct();
    for ind = 1:length(passed_inds)
        EM_PR_combined.last_metab_concs(:, ind) = ...
            EM_PR_all{passed_inds(ind)}.last_metab_concs;
        
        EM_PR_combined.fitness_values(ind, :) = ...
            EM_PR_all{passed_inds(ind)}.fitness_values;
        
        EM_PR_combined.avg_fitness_values(ind, 1) = ...
            EM_PR_all{passed_inds(ind)}.avg_fitness_values;
        
        EM_PR_combined.timepoints(:, ind) = ...
            EM_PR_all{passed_inds(ind)}.timepoints;
        
        if isfield(EM_PR_all{passed_inds(ind)}, 'ode_warn_flags')
            EM_PR_combined.ode_warn_flags(:, ind) = ...
                EM_PR_all{passed_inds(ind)}.ode_warn_flags;
        end
        
    end
        
case 'patternsearch'
    
    patternsearch_options = optimoptions('patternsearch');

    % Patternsearch requires initial points - no creation fcn
    if ~isempty(best_n_ksets_as_init)
        pattsearch_xo = InitialSwarmMatrix(:, 1)';
    elseif ~isempty(n_init_ksets)
        fprintf("Warning: Randomly sampling for patternsearch initial points\n")
       
        creationfcn_options.PopulationSize = 1;
        [Population, ~] = create_param_set_population(n_vars, [], creationfcn_options, ...
            model, Options, []);
        pattsearch_xo = Population;
        
    else
        error("Must specify either # of initial points or best prev points\n")        
    end
    
    if UseParallel
        patternsearch_options.UseParallel = true;
    end
    
    patternsearch_options.MaxTime = Options.MaxTime;
    patternsearch_options.TimeLimit = Options.MaxTime;
    
    if ~isempty(Options.MaxFunEvals)
        patternsearch_options.MaxFunEvals = Options.MaxFunEvals;
        patternsearch_options.MaxFunctionEvaluations = Options.MaxFunEvals;
    end
    
    if ~isempty(Options.MaxIterations)
        patternsearch_options.MaxIterations = Options.MaxIterations;
        patternsearch_options.MaxIter = Options.MaxIterations;
    end
    
    starttime = tic;
    
    % Uncomment to show plotted output
%     patternsearch_options.PlotFcn = {'psplotbestf', 'psplotmeshsize'};
    
    if run_fake_optimization % TODO - delete before publication
        
        [fit, PR] = obj_fxn(pattsearch_xo(:));
        opt_params = pattsearch_xo(:);
        exitflag = NaN;
        output = NaN;
        
    else
    
        [opt_params, opt_fitness, exitflag, output] = ...
            patternsearch(obj_fxn, pattsearch_xo', [], [], A_eq, b_eq, ...
            Options.lb, Options.ub, nonlcon, patternsearch_options);

    end
    
    options_struct = patternsearch_options;
   
case 'PSO' % Particle Swarm Optimization (MATLAB toolbox)
    
    % Standard PSO (only bound constraints)
    pso_options = optimoptions('particleswarm');
    
    if ~isempty(n_init_ksets)
        pso_options.SwarmSize = n_init_ksets;
    end

    if ~isempty(best_n_ksets_as_init)
        pso_options.InitialSwarmMatrix = InitialSwarmMatrix';
    end
    
    if UseParallel
        pso_options.UseParallel = true;
    end
    
    pso_options.MaxStallIterations = 1e3;
    
    pso_options.MaxTime = Options.MaxTime;
    
    starttime = tic;
    
    pso_options.PlotFcn = {'pswplotbestf'};
    
    [opt_params, opt_fitness, exitflag, output] = ...
        particleswarm(obj_fxn, n_vars, Options.lb, Options.ub, pso_options);
    
    options_struct = pso_options;
    
case 'PSO_const' % Constrained (linear/nonlinear ) Particle Swarm (user-written)    

    if exist('pso') ~= 2
        error("To use this, install https://www.mathworks.com/matlabcentral/fileexchange/25986-constrained-particle-swarm-optimization")
    end
    % This user-made function uses GA settings
    pso_options = optimoptions('ga');

    if ~isempty(n_init_ksets)
        pso_options.PopulationSize = n_init_ksets;
    end

    if ~isempty(best_n_ksets_as_init)
        pso_options.InitialPopulation = InitialSwarmMatrix';
    end
    
    if UseParallel
        pso_options.UseParallel = 'always';
    else
        pso_options.UseParallel = 'never';
    end
    
    % Increase number of stall generations
    pso_options.StallGenLimit = 200;

    starttime = tic;
    
    pso_options.CreationFcn = @(GenomeLength, FitnessFxn, ga_options) ...
        create_param_set_population(GenomeLength, FitnessFxn, ga_options, model, Options, []);
    
    % Set up plot functions - just leave on by default
    pso_options.PlotFcn = {'pswplotbestf'};
    pso_options.PlotFcn = {'psoplotbestf'};
    
    [opt_params, opt_fitness, exitflag, output] = ...
        pso(obj_fxn, n_vars, [], [], A_eq, b_eq, ...
        Options.lb, Options.ub, nonlcon, pso_options);
    
    options_struct = pso_options;
    
    pso_savestruct.opt_params = opt_params;
    pso_savestruct.opt_fitness = opt_fitness;
    pso_savestruct.output = output;
    pso_savestruct.runtime = toc(starttime);
    save(savename, '-struct','pso_savestruct')
    
case 'GA'
    % Option to do genetic algorithm
    ga_options = optimoptions('ga');
        
    if ~isempty(best_n_ksets_as_init)
        % Save into options struct
        ga_options.InitialPopulationMatrix = InitialSwarmMatrix';
        
        % Add in scores because why not
        if exist('init_ksets_fits','var')
            ga_options.InitialScoresMatrix = init_ksets_fits(:);
        end
    end
    
    if ~isempty(n_init_ksets)
    	ga_options.PopulationSize = n_init_ksets;
    end
    
    if UseParallel
        ga_options.UseParallel = true;
    end
    
    % Set crossover and creation fxns
    ga_options.CreationFcn = @(GenomeLength, FitnessFxn, ga_options) ...
        create_param_set_population(GenomeLength, FitnessFxn, ga_options, model, Options, []);
    
    ga_options.CrossoverFcn = @(parents, ga_options, nvars, FitnessFcn,...
        unused, thisPopulation) ...
        ga_crossover_fxn (parents, ga_options, nvars, FitnessFcn, unused, thisPopulation,...
        model, xover_ratio, Options);
    
    % Set timeouts
    ga_options.MaxTime = Options.MaxTime;
    ga_options.TimeLimit = Options.MaxTime;
    
    starttime = tic;
    
    % Show plots
    ga_options.PlotFcns = {@gaplotbestf,@gaplotscores,@gaplotstopping,@gaplotbestindiv,@gaplotrange};
    ga_options.Display = 'iter';
    
    % Set arbitrary limits
    ga_options.MaxStallGenerations = 200;
    ga_options.StallGenLimit = 200;
    % Limits from varargin
    ga_options.NonlinearConstraintAlgorithm = Options.NonlinConAlgorithm;
    ga_options.NonlinConAlgorithm = Options.NonlinConAlgorithm;
    if ~isempty(Options.PenaltyFactor)
        % Only used for auglag (I think) % default is 100
        ga_options.PenaltyFactor = Options.PenaltyFactor; 
    end
    if ~isempty(Options.InitialPenalty)
        % Only used for auglag (I think) % default is 10
        ga_options.InitialPenalty = Options.InitialPenalty;
    end

    % Call GA
    [opt_params, opt_fitness, exitflag, output] = ...
        ga(obj_fxn, n_vars, [], [], A_eq, b_eq,...
        Options.lb, Options.ub, nonlcon, ga_options);
    
    options_struct = ga_options;
    
    ga_savestruct.opt_params = opt_params;
    ga_savestruct.opt_fitness = opt_fitness;
    ga_savestruct.output = output;
    ga_savestruct.runtime = toc(starttime);
    save(savename, '-struct','ga_savestruct')
    
otherwise

    error("This algorithm is not valid")
    
end

%% If log-transformed variables were used, re-transform and check for slack vars
rerun_Options = Options;
if Options.use_log_params
        
    % If some parameters were not optimized, add those back in
    if (isfield(Options, 'ignore_unused_params') && Options.ignore_unused_params) || ...
        (isfield(Options, 'param_inds_for_optimization') && ...
            ~isempty(Options.param_inds_for_optimization))
        
        % Add in to correct indices
        full_param_vector = Options.fixed_param_values;
        variable_param_inds = Options.variable_param_inds;
        full_param_vector(variable_param_inds) = opt_params;
        % Reassign full param vector
        opt_params = full_param_vector;

        % Turn off options to use again
        rerun_Options.fixed_param_values = [];
        
    end
    
    % Revert from log values
    opt_params = exp(opt_params);
    
    % If log & "nonlinear params", remove slack variables
    if Options.use_nonlincon
        n_net_rxns = size(model.S, 2);
        n_params = length(opt_params);
        n_slack_vars = n_net_rxns;
        slack_var_inds = n_params - n_slack_vars + 1 : n_params;
        opt_params(slack_var_inds) = [];
    end
    rerun_Options.use_log_params = false;
end

%% Rerun single best/optimized model
[opt_fitness, opt_PR] = calc_rates_and_fitness(...
    model, Experimental_Data, opt_params, rerun_Options);

%% Store primary fields into output struct
output_struct.input_parser = p;
output_struct.varargin = varargin;
output_struct.model = model;
output_struct.Options = Options;
output_struct.load_file = load_file;
output_struct.all_bounds = all_bounds;
output_struct.opt_params = opt_params;
output_struct.opt_fitness = opt_fitness;
output_struct.exitflag = exitflag;
output_struct.output = output;
output_struct.opt_PR = opt_PR;
output_struct.optim_options_struct = options_struct;
output_struct.runtime = toc(starttime);
output_struct.Experimental_Data = Experimental_Data;

%% Check which fields were used and save relevant ones in output struct
if ~isempty(best_n_ksets_as_init)
    output_struct.init_swarm = InitialSwarmMatrix;
end

if exist('optim_output_history', 'var')
    output_struct.optim_output_history = optim_output_history;
end

if exist('PR_litvars','var')
    output_struct.PR_litvars = PR_litvars;
end

if exist('EM_fits','var') && (keep_all_fits || keep_only_fits)
   output_struct.EM_fits = EM_fits;
end

if exist('EM_PR_combined','var')
   output_struct.EM_PR_combined = EM_PR_combined;
end

if exist('kept_kset_inds','var')
    output_struct.kept_kset_inds = kept_kset_inds;
end

if exist('kset_runtimes','var') && keep_all_fits && ~keep_only_fits
   output_struct.EM_kset_times = kset_runtimes;
end

if exist('Population','var') && save_EM_population
    output_struct.EM_population = Population;
end

if exist('best_kset_inds','var')
    output_struct.best_kset_inds = best_kset_inds;
end

if exist('best_kset_rng_seed','var')
    output_struct.best_kset_rng_seed = best_kset_rng_seed;
end

if exist('EM_fits_all','var') && keep_all_fits
    output_struct.EM_fits_all = EM_fits_all;
end
    
%% Save if optioned
if ~isempty(savename)
    save(savename, '-struct', 'output_struct','-v7')
end



end


