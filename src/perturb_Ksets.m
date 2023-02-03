% PERTURB_KSETS Calculate flux solutions for each kset at a given
%   condition
%
%  [solutions, complete_times, ode_warn_flags, slope_norms, v_elem, ...
%     x_final, all_x, ME_row, perturbed_stoich, Vinit] = ... 
%     perturb_Ksets(model, Experimental_Data, All_K, All_fractions,...
%     yth_perturb, enz_perturbed_rxns, ...
%     enz_comps_per_pert_rxn, expression_levels, Options)
%
% model                   - Model structure
% Experimental_Data       - Struct with Experimental Data
% yth_perturb             - Integer number of condition currently being
%                           used, relative to the index of cells 
%                           (flux_rxns, ref_fluxes, or
%                           constant_metab_levels, etc.) within local
%                           Experimental_Data struct
% yth_enzymes_changed     - Indices in model.rxns of reactions involved in
%                           over/underexpression in this condition.
%                           Represented as column vector.
% enz_comps_per_pert_rxn  - OBSOLETE - removed on 2019-02-26.
%                           List of indices in
%                           model.metabs_and_enzyme_complexes for all
%                           complexes involved in over/underexpression in
%                           this condition.  For now, only does a single
%                           over/underexpression, but will need to be
%                           updated to do multiple in a single condition -
%                           this will require changes in
%                           POST_sensitivity_analysis.m as well - likely
%                           need to make this variable a cell array
% yth_enzyme_levels       - Fold change for all complexes listed in
%                           enz_comps_per_pert_rxn.  Column vector with
%                           same number of elemens as yth_enzymes_changed
% TODO

function [complete_times, ode_warn_flags, all_x, all_t] = ...
    perturb_Ksets(model, Experimental_Data, Initial_Ensemble,...
    yth_perturb, yth_enzymes_changed, yth_enzyme_levels, Options)

All_K = Initial_Ensemble.K_sets;

n_ksets = size(All_K,2);
complete_times = zeros(1,n_ksets);
ode_warn_flags = zeros(n_ksets,1);
all_x = cell(1,n_ksets);
all_t = cell(1,n_ksets);

%% Account for constant metabolites
% Get indices and levels for metabolites held constant in this perturbation
yth_constant_metab_levels = ...
    Experimental_Data.constant_metab_levels{yth_perturb};
yth_constant_metab_indices = ...
    Experimental_Data.constant_metab_indices{yth_perturb};

% alter K's to include concentrations of constant metabolites
% (fixed-timecourse is done within ODE solver)
% Only do for "mass-action" reactions (not saturation)
mass_action_net_rxns = transpose(find(ismember(model.rxn_type, [1,2,6,7,8,11])));
mass_action_elem_rxns = [];
for ma = mass_action_net_rxns
    mass_action_elem_rxns = [mass_action_elem_rxns, ...
        model.rxn_indices(ma, 1) : model.rxn_indices(ma, 2)];
end
for i=1:length(yth_constant_metab_levels)
    % Get indices of elementary reactions where constant metabs are
    % reactants
    reactant_elem_rxns = ...
        find(model.S_f_b(yth_constant_metab_indices(i), :) < 0);
    % Only keep those from mass action elem rxns
    reactant_elem_rxns = intersect(reactant_elem_rxns, mass_action_elem_rxns);
    % Adjust K's
    for j=1:length(reactant_elem_rxns)
        All_K(reactant_elem_rxns(j),:) = ...
            All_K(reactant_elem_rxns(j),:) .* yth_constant_metab_levels(i);
    end
end

% Get fields for fixing metab timecourses across timepoints (i.e., pH)
if ~isempty(Experimental_Data.fixed_metab_timecourse{yth_perturb})
    
    fixed_met_inds = Experimental_Data.fixed_metab_timecourse{yth_perturb}(:, 1);
    if iscell(fixed_met_inds)
        [~, fixed_met_inds] = ismember(fixed_met_inds, model.mets);
    end

    % Get timepoints and concs
    fixed_met_timepoints = cell2mat(Experimental_Data.fixed_metab_timecourse...
        {yth_perturb}(:, 2));
    fixed_met_concs = cell2mat(Experimental_Data.fixed_metab_timecourse...
        {yth_perturb}(:, 3));
else
    fixed_met_inds = [];
    fixed_met_timepoints = [];
    fixed_met_concs = [];
end


%% Get list of enzyme and enzyme complexes associated with chnaged rxn y 

% Cell array of which metab/enz/complexes correspond with each reaction
enz_enz_comp_indices = cell(1,length(yth_enzymes_changed));

% If a cell of enzyme_changes has multiple enzyme levels to
% change in this perturbation, go through each of them
n_enz_perturbed = length(yth_enzymes_changed);
% Keep track of which are adjusted from Experimental_Data fields (if doing
% approx model & converting kcat to vmax, don't want to adjust twice)
elem_rxns_already_adjusted = [];

% Don't perturb for regulation reactions affecting approx rate forms - get
% those
reg_rxn_net_inds = find(ismember(model.rxn_type, [6,7,8]));
% Figure out which of these are for approx rate form rxns
reg_rxns_on_approx_rates = [];
for reg_rxn = reg_rxn_net_inds(:)'
    shared_rxns = find(model.enzyme_rxn_specificity(:, reg_rxn));
    if any(ismember(model.rxn_type(shared_rxns), [5,50,51]))
        reg_rxns_on_approx_rates = [reg_rxns_on_approx_rates; reg_rxn];
    end
end

for enz = 1:n_enz_perturbed
    
    added_indices = [];
    yth_enz_enz_comp_indices = ...
        find(model.enz_enzComplex(:,yth_enzymes_changed(enz)) ~=0);

    % If this is empty (will happen if irreversible/approximate form
    % rxn was adjusted), then adjust K's for this reaction instead
    if isempty(yth_enz_enz_comp_indices)

        % get elementary reaction indices associated with this reaction (or
        % others with same enzyme)
        enz_net_rxns = find(model.enzyme_rxn_specificity(:,yth_enzymes_changed(enz)));
        % Remove reg rxns for approx rate form rxns
        enz_net_rxns = setdiff(enz_net_rxns, reg_rxns_on_approx_rates);
        % Track which elementary reactions
        elem_rxns = [];
        for enz_rxn = enz_net_rxns(:)'
            elem_rxns = [elem_rxns, model.rxn_indices(enz_rxn,1) : ...
                 model.rxn_indices(enz_rxn,2)];
        end
        
        % Adjust enzymes changed relative to reference level 
        ref_enz_levels = Experimental_Data.enzyme_levels...
            {Options.ref_from_ED_experimental_condition};
        % Get relative level
        ref_relative_adj = yth_enzyme_levels(enz) ./ ref_enz_levels(enz);
        % If NaN (ref level was zero), just set K to zero
        if isnan(ref_relative_adj)
            ref_relative_adj = 0;
        end
        % Adjust k's
        All_K(elem_rxns,:) = All_K(elem_rxns,:) .* ref_relative_adj;                
        
        % Track which were adjusted
        elem_rxns_already_adjusted = [elem_rxns_already_adjusted, elem_rxns];

        % Still need to pass in enz_enz_comp_indices
        enz_enz_comp_indices{enz} = [];
        
    else
    
        %first element will be index of free enzyme
        free_enz_ind_y = yth_enz_enz_comp_indices(1);

        %find all reactions which use that free enzyme
        rxns_with_enz_y = find(model.enz_enzComplex(free_enz_ind_y,:) ~= 0);

        %for each of those reactions, find the metab/enz/complex indices
        %then get only unique values

        for rxn = rxns_with_enz_y
            new_indices = find(model.enz_enzComplex(:,rxn) ~= 0);
            added_indices = [added_indices; new_indices];
        end
        
        %save unique sorted elements into column of that reaction
        enz_enz_comp_indices{enz} = unique(added_indices);
        
    end

end
    
% Set init concentrations
all_initial_metab_concs = Initial_Ensemble.initial_metab_concs;


%% Adjust initial concentrations based on condition-specific info in
%% Experimental_Data
% Get values from `initial_metab_concentrations`
ED_init_concs = Experimental_Data.initial_metab_concs{yth_perturb};
ED_met_inds = ED_init_concs(:, 1);
ED_met_concs = ED_init_concs(:, 2);
if iscell(ED_met_inds)
    [~,ED_met_inds] = ismember(ED_met_inds, model.mets);
    ED_met_concs = cell2mat(ED_met_concs);
end
all_initial_metab_concs(ED_met_inds) = ED_met_concs;

% overwrite any t=0 values in `metab_timecourse`
met_tc = Experimental_Data.metab_timecourse{yth_perturb};
for m = 1:size(met_tc, 1)
    t_zero_tspan_ind = find(met_tc{m, 2}(1, :) == 0);
    % If there is a t=0
    if ~isempty(t_zero_tspan_ind)
        [~,t_zero_met_ind] = ismember(met_tc{m, 1}, model.mets);
        t_zero_met_conc = met_tc{m,2}(2, t_zero_tspan_ind);
        all_initial_metab_concs(t_zero_met_ind) = t_zero_met_conc;
    end
end

% Get relative levels of enzymes in this condition compared to reference
yth_relative_enz_levels = yth_enzyme_levels ./ ref_enz_levels;
% Replace nans with zeros
yth_relative_enz_levels( isnan(yth_relative_enz_levels) ) = 0;

%% Loop over k-sets and call calculate_rates

% for x = 1:n_ksets % For debugging
parfor x = 1:n_ksets

    %initial concentration of all metabolites is 1; fractions are
    % different for each kset
    initial_conc_perturbed = all_initial_metab_concs(:, x);

    % assign indices for enzyme over/underexpression
    for enz = 1:n_enz_perturbed

        % multiply all elements at indices in column vector of
        % [enz_enz_comp_indices{enz}] by the single value at
        % [yth_enzyme_levels(enz)]
        initial_conc_perturbed(enz_enz_comp_indices{enz}) = ...
            initial_conc_perturbed(enz_enz_comp_indices{enz}).*yth_relative_enz_levels(enz);

    end

    [complete_times(x), ode_warn_flags(x), all_x{1,x}, all_t{x}] = ...
    calculate_rates(model, All_K(:,x), initial_conc_perturbed, Options,...
        yth_constant_metab_levels, fixed_met_inds, fixed_met_timepoints, fixed_met_concs);
    
end %end parfor loop
   
    
end
