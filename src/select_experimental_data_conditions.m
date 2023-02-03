% SELECT_EXPERIMENTAL_DATA_CONDITIONS Trim down which experimental data
% conditions are used from the larger set.  Called a few different places
% so useful to have a function
%
% [trimmed_Exp_Data] = select_experimental_data_conditions(...
%   Experimental_Data, new_conditions)

% Revision History:
%{
    
2019-07-03: jpm
    Created script
2019-10-08
    Default to passing all variables over to new struct
2019-12-05
    Added flux_reps and flux_error fields if present
2020-03-13
    Added initial_metab_concs if present
%}

function [trimmed_Exp_Data, trimmed_Options] = select_experimental_data_conditions(...
    Experimental_Data, Options, new_conditions)

% save old struct
Experimental_Data.All_Exp_Data = Experimental_Data;

% If new_conditions is explicitly passed in, use that, otherwise use field
% in Options - should really just use Options field
if ~exist('new_conditions','var')
    new_conditions = Options.experimental_conditions;
end

trimmed_Exp_Data = Experimental_Data;

trimmed_Exp_Data.ref_fluxes = Experimental_Data.ref_fluxes(new_conditions);
trimmed_Exp_Data.flux_rxns = Experimental_Data.flux_rxns(new_conditions);


trimmed_Exp_Data.constant_metab_indices = ...
    Experimental_Data.constant_metab_indices(new_conditions);
trimmed_Exp_Data.constant_metab_levels = ...
    Experimental_Data.constant_metab_levels(new_conditions);

if isfield(Experimental_Data,'flux_error')
   trimmed_Exp_Data.flux_error = ...
       Experimental_Data.flux_error(new_conditions);
end

if isfield(Experimental_Data,'flux_reps')
   trimmed_Exp_Data.flux_reps = ...
       Experimental_Data.flux_reps(new_conditions);
end

if isfield(Experimental_Data,'enzyme_levels')
    trimmed_Exp_Data.enzyme_levels = ...
        Experimental_Data.enzyme_levels(new_conditions);
    trimmed_Exp_Data.enzymes_changed = ...
        Experimental_Data.enzymes_changed(new_conditions);
elseif isfield(Experimental_Data,'expression_levels')
    
    fprintf('WARNING: Exp_Data struct field names are outdated\n')
    
    trimmed_Exp_Data.enzyme_levels = ...
        Experimental_Data.expression_levels(new_conditions);
    trimmed_Exp_Data.enzymes_changed = ...
        Experimental_Data.flux_KOs(new_conditions);
else
    error('Exp_Data struct doesnt have enzyme level information')
end

if isfield(Experimental_Data,'initial_metab_concs')
    trimmed_Exp_Data.initial_metab_concs = ...
        Experimental_Data.initial_metab_concs(new_conditions);
end

if isfield(Experimental_Data,'metab_timecourse')
    trimmed_Exp_Data.metab_timecourse = ...
        Experimental_Data.metab_timecourse(new_conditions);
end

if isfield(Experimental_Data, 'condition_names')
    trimmed_Exp_Data.condition_names = ...
        Experimental_Data.condition_names(new_conditions);
end

if isfield(Options,'ref_from_ED_experimental_condition')
    new_ref_condition = find(isequal(...
        new_conditions,Options.ref_from_ED_experimental_condition));
        
    if ~isempty(new_ref_condition)
        
        Options.ref_from_ED_experimental_condition = ...
            new_ref_condition;
    else
        
%         error("Options.ref_from_ED_experimental_condition is not in kept conditions")
        fprintf("WARNING: Options.ref_from_ED_experimental_condition is not in kept conditions\n")
        % Just make into the first (new) experimental condition
        Options.ref_from_ED_experimental_condition = 1;
    end
    
end

% Note: ext_metab_index may be only one column - no re-ordering
% If it is a cell or matrix (not just a single vector for all conditions)
if isfield(Experimental_Data,'ext_metab_indices')
    
    if iscell(Experimental_Data.ext_metab_indices)
        
        trimmed_Exp_Data.ext_metab_indices = ...
            Experimental_Data.ext_metab_indices(new_conditions);
        trimmed_Exp_Data.ext_metab_levels = ...
            Experimental_Data.ext_metab_levels(new_conditions);
    
    elseif size(Experimental_Data.ext_metab_indices,2) > 1
    
        trimmed_Exp_Data.ext_metab_indices = ...
            Experimental_Data.ext_metab_indices(:,new_conditions);
        trimmed_Exp_Data.ext_metab_levels = ...
            Experimental_Data.ext_metab_levels(:,new_conditions);

    else
        % TODO: this shouldn't happen - fix
        error("You shouldn't be here");
        trimmed_Exp_Data.ext_metab_indices = ...
            Experimental_Data.ext_metab_levels;
    end
end

% Check some others
if isfield(Experimental_Data,'fixed_metab_timecourse')
    trimmed_Exp_Data.fixed_metab_timecourse = ...
        Experimental_Data.fixed_metab_timecourse(new_conditions);
end
if isfield(Experimental_Data,'ref_metab_ranges')
    trimmed_Exp_Data.ref_metab_ranges = ...
        Experimental_Data.ref_metab_ranges(new_conditions);    
end
if isfield(Experimental_Data,'mdf_opt_concs')
    trimmed_Exp_Data.mdf_opt_concs = ...
        Experimental_Data.mdf_opt_concs(new_conditions);   
    
end
if isfield(Experimental_Data, 'condition_names')
    trimmed_Exp_Data.condition_names = ...
        Experimental_Data.condition_names(new_conditions);
end

% Things that are not cells
if isfield(Experimental_Data,'CV')
    trimmed_Exp_Data.CV = Experimental_Data.CV;   
end

if isfield(Experimental_Data,'weighted_rxn')
    trimmed_Exp_Data.weighted_rxn = Experimental_Data.weighted_rxn;   
end

% If screen_threshold has already been made, trim that too
trimmed_Options = Options;
if isfield(Options,'screen_threshold') &&...
        (length(Options.screen_threshold) > length(new_conditions))

    Options.untrimmed_screen_threshold = Options.screen_threshold;
    
    trimmed_Options.screen_threshold = Options.screen_threshold(new_conditions);

end

end

