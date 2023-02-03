% Set up everything needed for optimization and pass back as struct

function optimization_struct = setup_optimization(model, Experimental_Data, ...
    Options, params)

% Function is occasionally used without an input parameter vector - allow that
if ~exist('params','var')
    params = [];
end

% Adjust bounds based on varargin
if ~isempty(Options.optim_rate_const_bounds)
    Options.lb(Options.rate_const_inds) = Options.optim_rate_const_bounds(1);
    Options.ub(Options.rate_const_inds) = Options.optim_rate_const_bounds(2);
end

if ~isempty(Options.optim_sat_const_bounds)
    Options.lb(Options.sat_const_inds) = Options.optim_sat_const_bounds(1);
    Options.ub(Options.sat_const_inds) = Options.optim_sat_const_bounds(2);
end

% Get separate bounds for inhibition constants - use sat constant
% bounds if none specifically given for inhibition constants
if ~isempty(Options.optim_inhib_const_bounds)
    Options.lb(Options.reg_elem_rxns) = Options.optim_inhib_const_bounds(1);
    Options.ub(Options.reg_elem_rxns) = Options.optim_inhib_const_bounds(2);

elseif ~isempty(Options.optim_sat_const_bounds)
    fprintf("Warning: using saturation constant bounds for inhibition constants\n")
    Options.lb(Options.reg_elem_rxns) = Options.optim_sat_const_bounds(1);
    Options.ub(Options.reg_elem_rxns) = Options.optim_sat_const_bounds(2);
end

lb = Options.lb;
ub = Options.ub;

% If using previous points, set bounds to include that point
if ~isempty(params)
    
    % First check if manually imposing new init point
    if isfield(Options, 'optim_init_params') && ~isempty(Options.optim_init_params)
        params(Options.optim_init_params(:,1)) = Options.optim_init_params(:, 2);
    end
    
    x0_min = min(params, [], 2);
    x0_max = max(params, [], 2);
    Options.lb = min(lb, x0_min);
    Options.ub = max(ub, x0_max); 
end

% If optioned, add final constraints on optimization problem (even if it
% excludes the initial point)
if isfield(Options, 'optim_param_bounds') && ~isempty(Options.optim_param_bounds)
    param_bound_inds = Options.optim_param_bounds(:, 1);
    param_bound_lb = Options.optim_param_bounds(:, 2);
    param_bound_ub = Options.optim_param_bounds(:, 3);
    
    Options.lb(param_bound_inds) = param_bound_lb;
    Options.ub(param_bound_inds) = param_bound_ub;
    
    kcat_inds = 1 : size(model.S_f_b, 2);
    
    % Make init points fit within these bounds - throw warning if needed
    n_params = size(params, 2);
    for i = 1:n_params
        lb_viol_inds = find(params(:, i) < Options.lb);
        ub_viol_inds = find(params(:, i) > Options.ub);

        % Loop over lb violations
        for v = 1:length(lb_viol_inds)
            fprintf("Warning: lb constraint is violated in kset %i.\n", i)

            % Also get kcat in other direction and keep ratio the same
            if ismember(lb_viol_inds(v), kcat_inds)
                [net_rxn_ind, forward_rev_ind] = ...
                    find(model.rxn_indices == lb_viol_inds(v));

                if forward_rev_ind == 1
                    opp_fr_ind = 2;
                else
                    opp_fr_ind = 1;
                end

                opp_kcat_ind = model.rxn_indices(net_rxn_ind, opp_fr_ind);

                main_opp_param_ratio = params(lb_viol_inds(v), i) ./ ...
                    params(opp_kcat_ind, i);

                params(lb_viol_inds(v), i) = Options.lb(lb_viol_inds(v));

                params(opp_kcat_ind, i) = Options.lb(lb_viol_inds(v)) / main_opp_param_ratio;

            else
                error("Implement setting after-the-fact bounds on Km params")
            end

        end
            
        for v = 1:length(ub_viol_inds)
            fprintf("Warning: ub constraint is violated in kset %i.\n", i)

            % Also get kcat in other direction and keep ratio the same
            if ismember(ub_viol_inds(v), kcat_inds)
                [net_rxn_ind, forward_rev_ind] = ...
                    find(model.rxn_indices == ub_viol_inds(v));

                if forward_rev_ind == 1
                    opp_fr_ind = 2;
                else
                    opp_fr_ind = 1;
                end

                opp_kcat_ind = model.rxn_indices(net_rxn_ind, opp_fr_ind);

                main_opp_param_ratio = params(ub_viol_inds(v), i) ./ ...
                    params(opp_kcat_ind, i);

                params(ub_viol_inds(v), i) = Options.ub(ub_viol_inds(v));

                params(opp_kcat_ind, i) = Options.ub(ub_viol_inds(v)) / main_opp_param_ratio;

            else
                error("Implement setting after-the-fact bounds on Km params")
            end

        end % Loop over ub constraints
    end % Loop over ksets

end

% Check if using log-transformed parameters
if isfield(Options,'use_log_params') && Options.use_log_params
    % Don't do for EM
    if strcmpi(Options.opt_alg, 'EM')
        fprintf("Warning: can't do log params for EM\n")
    else
        % Won't use nonlcon if log params
        nonlcon = [];
        
        if Options.use_nonlincon
            % Get constraints
            [A, b, lb, ub, params] = get_log_param_haldane_constraints(...
                model, Options, params);
            % Keq values become slack params
        else
            A = [];
            b = [];
            lb = log(Options.lb);
            ub = log(Options.ub);
            params = log(params);
            
        end
    end
else
    % Set placeholder values
    A = [];
    b = [];
    % Nonlinear constraints
    if Options.use_nonlincon
        nonlcon = @(vars) calc_thermo_penalty(vars, model, Options);
    else
        nonlcon = [];
    end
    
end

% Calculate slack variable indices
n_final_params = length(params);

% If fixing/removing some variables, or only optimizing certain variables, do that here
if (isfield(Options, 'ignore_unused_params') && Options.ignore_unused_params) || ...
        (isfield(Options, 'param_inds_for_optimization') && ...
         ~isempty(Options.param_inds_for_optimization))
     
    % If automatically getting "unused" param inds
    if (isfield(Options, 'ignore_unused_params') && Options.ignore_unused_params)
    
        unused_param_inds = get_unused_param_inds(model, Experimental_Data, Options);
        all_param_inds = 1:n_final_params;
        used_param_inds = setdiff(all_param_inds, unused_param_inds);
        % Make sure not also trying to manually select param indices for
        % optimization
        if isfield(Options, 'param_inds_for_optimization') && ...
         ~isempty(Options.param_inds_for_optimization)
            fprintf("Warning: field `param_inds_for_optimization` will not be used\n")
        end
        
    elseif isfield(Options, 'param_inds_for_optimization') && ...
         ~isempty(Options.param_inds_for_optimization)
        % Otherwise, manually geting inds - if doing log params & slack
        % vars, those need to be included
        fprintf("Warning: must include slack variable indices in `param_inds_for_optimization\n")
        used_param_inds = Options.param_inds_for_optimization;
    end
    
    % Save things into Options
    Options.fixed_param_values = params;
    Options.variable_param_inds = used_param_inds;
    
    % Trim down param fields
    params = params(used_param_inds);
    lb = lb(used_param_inds);
    ub = ub(used_param_inds);
    A = A(:, used_param_inds);
end


optimization_struct.params = params;
optimization_struct.lb = lb;
optimization_struct.ub = ub;
optimization_struct.A = A;
optimization_struct.b = b;
optimization_struct.nonlcon = nonlcon;
optimization_struct.Options = Options;

end