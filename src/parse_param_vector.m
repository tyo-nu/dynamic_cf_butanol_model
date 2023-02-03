% Take a vector of parameters (i.e., for PSO) and parse into needd fields

function Initial_Ensemble = parse_param_vector(model, Experimental_Data,...
    Options, param_vector, condition_number)

% If "constant"/"set" params are passed in, set those up
if isfield(Options, 'fixed_param_values') && ~isempty(Options.fixed_param_values) 
    
    full_param_vector = Options.fixed_param_values;
    % Assume it has this field - what are indices in 'full' param vector for
    % those variables that are being optimize
    variable_param_inds = Options.variable_param_inds;
    full_param_vector(variable_param_inds) = param_vector;
    % Reassign full param vector for parsing
    param_vector = full_param_vector;
end

% Check if log values
if isfield(Options,'use_log_params') && Options.use_log_params
    
    % convert to normal values
    param_vector = exp(param_vector);
    
    % Remove slack variables if present (only done if log-transformed)
    if isfield(Options, 'use_nonlincon') && Options.use_nonlincon 
        n_net_rxns = size(model.S, 2);
        n_params = length(param_vector);
        n_slack_vars = n_net_rxns;
        slack_var_inds = n_params - n_slack_vars + 1 : n_params;
        param_vector(slack_var_inds) = [];
    end
end

% Set up Initial_Ensemble
[n_metabs_and_fractions, n_ks] = size(model.S_f_b);
kset = param_vector(1 : n_ks);
Initial_Ensemble.K_sets = kset(:);
% Save param vector as easy way to convert back log params, if used
Initial_Ensemble.param_vector = param_vector;

n_metabs = size(model.S, 1);

% % % if ~isfield(Options, 'vary_saturation_constants') || Options.vary_saturation_constants

    n_sat_consts = length(param_vector) - n_ks;

    sat_const_vals = param_vector(n_ks + 1 : n_ks + n_sat_consts);
% % % else
% % %     sat_const_vals = model.sat_consts;
% % %     n_sat_consts = length(sat_const_vals);
% % % end

% Parse out saturation constants
% Get indices in sat_constants field (n_rxns_f_b x n_mets) where each
% constant belongs

% % % sat_constants_matrix_locs = [];
% % % mult_sat_rxns = find(model.rxn_type == 5);
% % % mult_sat_rxns = mult_sat_rxns(:)';
% % % for r = mult_sat_rxns
% % %     % For forward and reverse
% % %     elem_rxn_inds = model.rxn_indices(r, 1) : model.rxn_indices(r, 2);
% % %     for e = elem_rxn_inds
% % %         substrate_inds = find(model.S_f_b(:, e) < 0);
% % %         for s = 1:length(substrate_inds)
% % %             sat_constants_matrix_locs = [sat_constants_matrix_locs; ...
% % %                 [e, substrate_inds(s)] ];
% % %         end
% % %     end
% % % end

% Don't get these from stoich - less reliable (doesn't work with Varner
% model)
sat_constants_matrix_locs = model.lit_values.sat_constants_matrix_locs;

% Make sure right number of values are there
n_subs = size(sat_constants_matrix_locs, 1);
if length(sat_const_vals) < n_subs
    error("Warning: not enough saturation constants given for all substrates in model\n")
    % Replace remaining (missing) saturation constants with 1's 
    % TODO - make this an option
    sat_const_vals(end+1 : n_subs) = 1;
    
end
sat_constants_matrix = NaN(n_ks, n_metabs);
for row = 1:n_subs
   
    % Get index in flattened matrix
    matrix_row = sat_constants_matrix_locs(row,1);
    matrix_col = sat_constants_matrix_locs(row,2);
    flat_ind = (matrix_col - 1)*n_ks + matrix_row;
    
    
    % Add values into matrix
%     sat_constants_matrix(sat_constants_matrix_locs(row, :)) = sat_const_vals(row);
    sat_constants_matrix(flat_ind) = sat_const_vals(row);
    
end
   
% get initial concentrations from model struct (set earlier)
Initial_Ensemble.initial_metab_concs = model.init_concs(:);

init_fractions = param_vector(n_ks + n_sat_consts + n_metabs + 1 : end);
Initial_Ensemble.fractions = init_fractions(:);

Initial_Ensemble.sat_constants = sat_constants_matrix;

end