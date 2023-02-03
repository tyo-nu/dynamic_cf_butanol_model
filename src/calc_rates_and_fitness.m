% CALC_RATES_AND_FITNESS Solve ODEs for a single k-set and intial
% conditions and then calculate fitness to Experimental Data
% Meant to be used with self-contained optimization functions

% TODO

function [fitness_score, PR] = calc_rates_and_fitness(...
    model, Experimental_Data, kset_and_init_concs, Options)

% Set up Initial_Ensemble
Initial_Ensemble = parse_param_vector(model, Experimental_Data, Options,...
    kset_and_init_concs);
% Save vars into model struct that need to be there
model.lit_values.sat_constants = Initial_Ensemble.sat_constants;

% Set up any Options needed
Options.ode_solver = str2func('ode15s');

% Fields that could be changed
if ~isfield(Options,'fitness_fxn') || isempty(Options.fitness_fxn)
    Options.fitness_fxn = 'RMSE';
end

% Just call calculate_perturbed_steady_states to get everything
[fvals, all_complete_times, ...
    all_ode_warn_flags, last_metab_concs, all_time_points] = ... 
    calculate_perturbed_steady_states(model, Experimental_Data, ...
    Initial_Ensemble, Options);

% Make a struct ("Perturbation Results") to pass back
PR = struct();
PR.last_metab_concs = last_metab_concs;
PR.fitness_values = fvals;
PR.avg_fitness_values = mean(fvals);
PR.timepoints = all_time_points;
PR.complete_times = all_complete_times;
PR.ode_warn_flags = all_ode_warn_flags;

% Take mean fitness across conditions for this parameter set
fitness_score = mean(fvals);

end