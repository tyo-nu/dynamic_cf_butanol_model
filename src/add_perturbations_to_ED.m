% Take in cell array of perturbed initial concentrations and add to
% experimental data

function new_ED = add_perturbations_to_ED(...
        model, Experimental_Data, Options)
   
% Make sure that 
conc_pert = Options.initial_concentration_perturbations;

assert(ismember(size(conc_pert, 2), [2,3]), ...
    "Field `initial_concentration_perturbations` must have 2 columns (or 3 if including names)")

[n_perturbations, n_cols] = size(conc_pert);

all_perturbed_vars = conc_pert(:, 1);
all_perturbed_levels = conc_pert(:, 2);

if n_cols == 3
    perturbation_names = conc_pert(:, 3);
else
    perturbation_names = sprintfc("Perturbation %i", 1 : n_perturbations);
end

% Get existing condition(s) - replicate all existing conditions for all new
%   initial concentration perturbations
n_existing_conditions = length(Experimental_Data.initial_metab_concs);

% Make a copy of the old ED condition to move this initial conc
% info into
% `copy_conditions` field copies conditions from the left column into the
% right column
%   Want to take all orig conditions (left column is
%   1:n_existing_conditions), right column will alternate with new
%   perturbations
copy_conditions = {1, 1:n_perturbations};
for i=2:n_existing_conditions
    new_cond_nums = (i-1)*n_perturbations + 1 : i*n_perturbations;
    copy_conditions = [copy_conditions; {i, new_cond_nums}];
end
    
new_ED = zz_add_rm_experimental_conditions(Experimental_Data, ...
    'copy_conditions', copy_conditions);

% Loop through new perturbations
for p = 1:n_perturbations
   
    % Replicate for each existing condition
    for c = 1:n_existing_conditions
        
        % Get list of metabs for this perturbation
        pert_mets_or_rxns = string(all_perturbed_vars{p});
        pert_levels = all_perturbed_levels{p};
        assert(length(pert_mets_or_rxns) == length(pert_levels), ...
            "These need to be equal")
        
        % For each metab in this (new) list, see if it's in this condition
        starting_init_concs = Experimental_Data.initial_metab_concs{c};
        new_init_concs = starting_init_concs;
        
        % Do same for enzyme levels
        new_enz_changed = Experimental_Data.enzymes_changed{c};
        new_enz_levels = Experimental_Data.enzyme_levels{c};
        
        % For each perturbed species in this perturbation
        for m = 1:length(pert_mets_or_rxns)
            
            % Check if perturbation is concentration or enzyme level
            metab_hit = ismember(pert_mets_or_rxns(m), model.mets);
            [enz_hit, rxn_index] = ismember(pert_mets_or_rxns(m), model.rxns);
            
            if metab_hit && enz_hit
                error("This shouldn't happen - ambiguity in metab vs rxn names")
            
            % If metab concentration changed
            elseif metab_hit
            
                % Check if already in existing ED field
                [~,row_in_prev_conc_list] = ismember(pert_mets_or_rxns(m), starting_init_concs(:, 1));

                % If found
                if row_in_prev_conc_list ~= 0
                    % Make sure only one hit?
                    assert(length(row_in_prev_conc_list) == 1, "How did I do this?")
                    % Add in concentration
                    new_init_concs{row_in_prev_conc_list, 2} = pert_levels(m);
                else
                    % If not found, add to end
                    new_init_concs = [new_init_concs; {char(pert_mets_or_rxns(m)), pert_levels(m)}];
                end
                
            % If enzyme level changed
            elseif enz_hit
                               
                % Check if already in existing ED field
                [~,index_in_prev_enz_changed] = ismember(pert_mets_or_rxns(m), ...
                    model.rxns(Experimental_Data.enzymes_changed{c}));

                % If found
                if index_in_prev_enz_changed ~= 0
                    % Make sure only one hit?
                    assert(length(index_in_prev_enz_changed) == 1, "How did I do this?")
                    % Add in concentration
                    new_enz_levels(index_in_prev_enz_changed) = pert_levels(m);
                else
                    % If not found, add to end (both inds and levels)
                    new_enz_changed = [new_enz_changed, rxn_index];
                    new_enz_levels = [new_enz_levels, pert_levels(m)];
                    
                end
                
            end
            
            
        end
        
        
        
        % Add into correct cell in expaneded/copied ED struct
        new_ED_cell_num = (c-1)*n_perturbations + p;
        
        new_ED.initial_metab_concs{new_ED_cell_num, 1} = new_init_concs;
        new_ED.enzymes_changed{new_ED_cell_num, 1} = new_enz_changed;
        new_ED.enzyme_levels{new_ED_cell_num, 1} = new_enz_levels;
        
        % If no changes, keep metab timecourse - otherwise, discard
        if isempty(pert_mets_or_rxns)
            new_metab_tc = Experimental_Data.metab_timecourse{c};
        else
            new_metab_tc = [];
        end
        new_ED.metab_timecourse{new_ED_cell_num, 1} = new_metab_tc;
        
        % Add in name
        new_name = perturbation_names{p};
        if isfield(Experimental_Data, 'condition_names')
            new_name = sprintf('%s - %s',...
                Experimental_Data.condition_names{c}, new_name);
        end
        new_ED.condition_names{new_ED_cell_num, 1} = string(new_name);
        
        
    end
        
    
end
    
end