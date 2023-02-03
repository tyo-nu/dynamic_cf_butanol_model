% If optimizing using log values of parameters, return A & b matrices for
% linear constraints for optimization methods

function [A, b, log_lb, log_ub, log_params_and_slack_vars] = get_log_param_haldane_constraints(...
    model, Options, params)

% Get info from model
[n_params, n_ksets] = size(params);
[n_species, n_rate_consts] = size(model.S_f_b);
n_net_rxns = size(model.S, 2);
sat_const_locs = model.lit_values.sat_constants_matrix_locs;
n_sat_consts = size(sat_const_locs, 1);
sat_rxns = transpose(find(ismember(model.rxn_type, [4,5,50])));
fast_eq_rxns = transpose(find(model.rxn_type == 11));

%% Get bounds on actual parameters
log_lb_params = log(Options.lb);
log_ub_params = log(Options.ub);

%% Get dG range bounds
dG_mean_std = model.dG_rxn_phys_unitless;
% Add some small wiggle room to any uncertainties less than some limit
min_dG_uncertainty = 1e-3;
dG_mean_std(dG_mean_std(:, 2) < min_dG_uncertainty, 2) = min_dG_uncertainty;
dG_min = dG_mean_std(:, 1) - dG_mean_std(:, 2);
dG_max = dG_mean_std(:, 1) + dG_mean_std(:, 2);

% Ignore for other reactions
thermo_constrained_rxns = [sat_rxns, fast_eq_rxns];
ignored_thermo_rxns = setdiff(1:n_net_rxns, thermo_constrained_rxns);
dG_min(ignored_thermo_rxns) = 0;
dG_max(ignored_thermo_rxns) = 0;

% For "irreversible" reactions, just make sure that calculated dG is also
% irreversible in same direction - don't care about how much
irrev_dG_limit = 20; % dG/RT - corresponds to 12 kcal/mol or 50 kJ/mol
irrev_forward_rxns = find(dG_max < -1*irrev_dG_limit);
irrev_reverse_rxns = find(dG_min > irrev_dG_limit);
dG_min(irrev_forward_rxns) = -1e3;
dG_max(irrev_forward_rxns) = -1.*irrev_dG_limit;
dG_min(irrev_reverse_rxns) = irrev_dG_limit;
dG_max(irrev_reverse_rxns) = 1e3;

%% Set up linear inequality A*x <= b

% Make part of A for rate constants so that each row is a net reaction with
% a +1 for the forward rate constant index and -1 for the reverse rate constant index
A_rate_consts = zeros(n_net_rxns, n_rate_consts);
A_sat_consts = zeros(n_net_rxns, n_sat_consts);
for rxn_ind = sat_rxns(:)'
    % Fill in rate constants
    rxn_elem_inds = model.rxn_indices(rxn_ind, :);
    A_rate_consts(rxn_ind, rxn_elem_inds) = [1, -1];
    
    % Fill in sat constants
    sat_const_locs_rows = find(ismember(sat_const_locs(:,1), rxn_elem_inds));
    rxn_substrate_inds = sat_const_locs(sat_const_locs_rows, 2);
    rxn_stoich = model.S(rxn_substrate_inds, rxn_ind);
    A_sat_consts(rxn_ind, sat_const_locs_rows) = rxn_stoich';
    
end

% Do fast equilibrium reactions
for fe_net_rxn_ind = fast_eq_rxns
    % Get elem indices for this reaction
    fe_elem_inds = model.rxn_indices(fe_net_rxn_ind, :);
    % Get dG and Keq
    fe_dG = dG_mean_std(fe_net_rxn_ind, 1);
    fe_Keq = exp(-1 * fe_dG);
    
    % Use Keq to set ratio of forward to reverse - 
    %   Keq = kf / kr
    kf_ind = fe_elem_inds(1);
    kr_ind = fe_elem_inds(2);
    % Set into A
    A_rate_consts(fe_net_rxn_ind, kf_ind) = 1;
    A_rate_consts(fe_net_rxn_ind, kr_ind) = -1;
    
    % Get minimum (forward vs reverse) rate constant allowed - use for
    % bounds (give an order of magnitude below, nothing above - for now)
    if isfield(Options, 'min_equilibrium_kinetic_rate')
        min_eq_rate = Options.min_equilibrium_kinetic_rate;
    else
        min_eq_rate = 1e5;
    end
    if fe_Keq > 1
        % If kr is the smaller, set that only
        log_lb_params(kr_ind) = log(min_eq_rate) / 10;
        log_lb_params(kf_ind) = -Inf;
    else
        log_lb_params(kf_ind) = log(min_eq_rate) / 10;
        log_lb_params(kr_ind) = -Inf;
    end
    % Set upper bounds to Inf
    log_ub_params(kr_ind) = Inf;
    log_ub_params(kr_ind) = Inf;
end

A_eye = -1.*eye(n_net_rxns);
A = [A_rate_consts, A_sat_consts, A_eye];

% Convert to log Keq bounds
Keq_min = exp(-1 * dG_max);
Keq_max = exp(-1 * dG_min);

% Make zeros for non-approx reactions
A(ignored_thermo_rxns, :) = 0;
log_Keq_min = log(Keq_min);
log_Keq_max = log(Keq_max);
log_Keq_min(ignored_thermo_rxns) = -Inf;
log_Keq_max(ignored_thermo_rxns) = Inf;

% Get Keq_calc from these params if given, & add onto params (as slack
% variables)
if ~isempty(params)
    log_params = log(params);
    n_ksets = size(params, 2);
    log_Keq_all = [];
    for i = 1:n_ksets
        [~,~,Keq_0] = calc_thermo_penalty(params(:,i), model, Options);
        log_Keq = log(Keq_0);
        % Set all non-approx rxns to 0 - don't care
        log_Keq(ignored_thermo_rxns) = 0;
        log_Keq_all = [log_Keq_all, log_Keq];
    end
    % Put together full vector
    log_params_and_slack_vars = [log_params; log_Keq_all];
    
else
    log_params_and_slack_vars = log(params);
end

% Put together full bounds
log_lb = [log_lb_params; log_Keq_min];
log_ub = [log_ub_params; log_Keq_max];

if any(log_params_and_slack_vars < log_lb) || any(log_params_and_slack_vars > log_ub)
    error("One of the initial points violates thermodynamics - please fix (somehow)")
end

b = zeros(size(A, 1), 1);

end