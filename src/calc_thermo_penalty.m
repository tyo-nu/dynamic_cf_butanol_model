% Function to calculate thermo penalty - use in nonlinear constraints
function [c, ceq, Keq_calc, rxn_violations] = calc_thermo_penalty(params, model, Options)

% takes in vector `params`, tries to satisfy:
%   ceq(params) = 0
%   c(params) <= 0

% First, if "constant"/"set" params are passed in, set those up
if isfield(Options, 'fixed_param_values') && ~isempty(Options.fixed_param_values)
    full_param_vector = Options.fixed_param_values;
    % Assume it has this field - what are indices in 'full' param vector for
    % those variables that are being optimize
    variable_param_inds = Options.variable_param_inds;
    full_param_vector(variable_param_inds) = params;
    % Reassign full param vector for parsing
    params = full_param_vector;
end

params = params(:);

% Don't use any strict equality constraints
ceq = 0;

% Ignore reactions with zero kcat for forward and reverse
zero_kcat_fr_rxns = [];

% Save calculated dG for each reaction
n_rxns = size(model.S, 2);
dG_calc = NaN(n_rxns, 1);
Keq_calc = zeros(n_rxns, 1);

% Get reactions with thermo constraints
% If specified, use only requested reactions
if isfield(Options, 'rxns_for_nonlcon') && ~isempty(Options.rxns_for_nonlcon)
    sat_rxns = Options.rxns_for_nonlcon;
else
    %Otherwise use all approx rate forms
    sat_rxns = transpose(find(ismember(model.rxn_type, [4,5,50])));
end

% Get dG range
dG_mean_std = model.dG_rxn_phys_unitless;
dG_min = dG_mean_std(:, 1) - dG_mean_std(:, 2);
dG_max = dG_mean_std(:, 1) + dG_mean_std(:, 2);

% Ignore for other reactions
ignored_thermo_rxns = setdiff(1:n_rxns, sat_rxns);
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

for r = sat_rxns
           
    % Get index of kcat forward and reverse
    kcat_inds = model.rxn_indices(r, :);
    
    % Get indices of substrates/products in this reaction
    subs_prod_inds = find(model.S(:, r));
    
    % Dont include species without thermo
    excluded_species = ["h_c","h_e"];
    excluded_inds = find(ismember(model.mets, excluded_species));
    subs_prod_inds = setdiff(subs_prod_inds, excluded_inds);
    
    % Set up which rows in matrix_locs each has
    subs_prod_matrix_locs_row = zeros(size(subs_prod_inds));
        
    % Get stoich
    subs_prod_stoich = model.S(subs_prod_inds, r);
    
    
    % Get rows of each substrate/product in sat_constants_matrics_locs
    matrix_locs = model.lit_values.sat_constants_matrix_locs;
    for i = 1:length(subs_prod_inds)
        subs_prod_matrix_locs_row(i) = find(...
            matrix_locs(:, 2) == subs_prod_inds(i) & ...
            ismember(matrix_locs(:, 1), kcat_inds) );
    end
    
    % Get indices in param vector
    n_rate_consts = size(model.S_f_b, 2);
    subs_prod_param_vector_inds = subs_prod_matrix_locs_row + n_rate_consts;
    
    % Calculate Keq from Haldane
    Km_prods_over_subs = prod(...
        params(subs_prod_param_vector_inds) .^ subs_prod_stoich );
    if params(kcat_inds(2)) == 0
        if params(kcat_inds(1)) == 0
            zero_kcat_fr_rxns = [zero_kcat_fr_rxns; r];
        end
        kcat_f_over_r = 0;
    else
       kcat_f_over_r =  params(kcat_inds(1)) / params(kcat_inds(2));
    end
    
    if kcat_f_over_r == 0 || ~isfinite(kcat_f_over_r)
        dG_calc(r) = 0;
        
    else
        Keq_calc(r) = kcat_f_over_r * Km_prods_over_subs;

        % Convert to dG
        dG_calc(r) = -1 .* log(Keq_calc(r));
    end   
    
    % TESTING
    if dG_calc(r) < dG_min(r) || dG_calc(r) > dG_max(r)
        a=1;       
        
    end
    
end

% Do fast equilibrium reactions too
fast_eq_rxns = transpose(find(model.rxn_type == 11));
for fe_net_rxn_ind = fast_eq_rxns
    % Get elem indices for this reaction
    fe_elem_inds = model.rxn_indices(fe_net_rxn_ind, :);
    
    % Calc Keq from ratio of forward to reverse - 
    %   Keq = kf / kr
    kf_ind = fe_elem_inds(1);
    kr_ind = fe_elem_inds(2);
    Keq_calc(fe_net_rxn_ind) = params(kf_ind) / params(kr_ind);
end

% Ignore reactions with zero in forward and reverse
dG_calc(zero_kcat_fr_rxns) = NaN;

% Across all reactions, get the value that is most out of bounds - set that
% absolute value to c
lb_violation = dG_min - dG_calc;
ub_violation = dG_calc - dG_max;

[lb_max_viol, lb_max_loc] = max(lb_violation);
[ub_max_viol, ub_max_loc] = max(ub_violation);

% Get max of either
max_violation = max(lb_max_viol, ub_max_viol);

c = max_violation;

% Also get the max violation (in either direction) for each reaction, &
% preserve sign
rxn_violations = zeros(n_rxns, 1);
for r = 1:n_rxns
    % If lb in violation, add in and flip sign
    if lb_violation(r) > 0
        rxn_violations(r) = lb_violation(r) * -1;
    
    % If ub in violation
    elseif ub_violation(r) > 0
        rxn_violations(r) = ub_violation(r);
    end
       
    % (If both are inside, leave as zero)
end
    
% TESTING
if c < 1.9
    a=1;
end

% Check for NaN's
if isnan(c) || isinf(c)
%     error("Nonlcon thermo check returning NaN - fix")
    % For some reason GA just makes these sometimes? Just give high penalty
    c = 1e0;
%     fprintf("NaN/Inf value\n")
end

% fprintf('%0.3e\n',c)

% % Removing for now - debugging something and wasnt using anyway
% if isfield(Options, 'thermo_fitness_penalty_multiplier')
%     thermo_penalty = Options.thermo_fitness_penalty_multiplier;
% else
%     thermo_penalty = 1;
% end
% if c > 0
%     penalty = c * thermo_penalty;
% else
%     penalty = 0;
% end



end