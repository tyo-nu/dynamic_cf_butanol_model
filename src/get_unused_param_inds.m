% Get the indices of parameters not used for this model

function unused_inds = get_unused_param_inds(model, Experimental_Data, Options)

%% Get "reverse" kcats for inhibition rxns of approx rate form rxns
reg_rxn_net_inds = find(ismember(model.rxn_type, [6,7,8]));
% Figure out which of these are for approx rate form rxns
reg_rxns_on_approx_rates = [];
for reg_rxn = reg_rxn_net_inds(:)'
    shared_rxns = find(model.enzyme_rxn_specificity(:, reg_rxn));
    if any(ismember(model.rxn_type(shared_rxns), [5,50,51]))
        reg_rxns_on_approx_rates = [reg_rxns_on_approx_rates; reg_rxn];
    end
end
% Convert to elementary rxns for rate constant inds, get reverse (second)
% steps
rev_elem_reg_inds = model.rxn_indices(reg_rxns_on_approx_rates, 2);

%% Get rate/saturation constants for any reactions KOd in all conditions 
%   (plus rate constants of mass-action "inhibition" elementary steps)
n_conditions = length(Experimental_Data.metab_timecourse);
zero_enz_level_rxns = cell(n_conditions, 1);
for c = 1:n_conditions
    zero_enz_level_rxns{c} = Experimental_Data.enzymes_changed{c}(...
        Experimental_Data.enzyme_levels{c} == 0);    
    
end

% If log-params & constrained (slack variables), include those too
slack_var_inds = [];
n_rate_consts = size(model.S_f_b, 2);
sat_const_matrix = model.lit_values.sat_constants_matrix_locs;
n_rate_sat_consts = n_rate_consts + size(sat_const_matrix, 1);
use_slack_vars = Options.use_log_params && Options.use_nonlincon;

% Check how many conditions each reaction has zero enzymes for
n_zero_conditions_per_rxn = zeros(length(model.rxns), 1);
for r = 1:length(model.rxns)
    for c = 1:n_conditions
        if ismember(r, zero_enz_level_rxns{c})
            n_zero_conditions_per_rxn(r) = n_zero_conditions_per_rxn(r) + 1;
        end
    end
end
% Find any that are in all
zero_enz_rxn_inds = find(n_zero_conditions_per_rxn == n_conditions);
% Get any reactions that use the same enzymes
zero_enz_rxn_inds_all = [];
enz_spec = model.enzyme_rxn_specificity;
for zr = zero_enz_rxn_inds(:)'
    zero_enz_rxn_inds_all = [zero_enz_rxn_inds_all; find(enz_spec(:, zr))];
    
    % If using slack variables, add on the indices in the param vector for
    % those
    if use_slack_vars
        slack_var_inds = [slack_var_inds; n_rate_sat_consts + zr];
    end
    
end
zero_enz_rxn_inds_all = unique(zero_enz_rxn_inds_all);

% Get rate constant inds (elementary) for each of these
unused_rate_const_inds = [];
for zra = zero_enz_rxn_inds_all(:)'
    unused_rate_const_inds = [unused_rate_const_inds; ...
        [model.rxn_indices(zra, 1) : model.rxn_indices(zra, 2)]' ];
end

% Get associated sat constants for each unused elementary reaction index
sat_const_matrix_rows = [];
for er = unused_rate_const_inds(:)'
    locs_rows = find(sat_const_matrix(:, 1) == er);
    % Add to list
    sat_const_matrix_rows = [sat_const_matrix_rows; locs_rows];
end
% Convert to indices in whole param vector
unused_sat_const_inds = sat_const_matrix_rows + n_rate_consts;

% Put all together
unused_inds = [rev_elem_reg_inds; unused_rate_const_inds; unused_sat_const_inds; slack_var_inds];
unused_inds = unique(unused_inds);

end