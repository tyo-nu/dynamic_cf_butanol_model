% CALCULATE_RATES Calculate fluxes and metabolite concentrations for a given
% set of parameters and initial condition
%
% [Vnet, complete_time, ode_warn_flag, norm_slope, vuni, x_final,...
%     x, t, ME, Vinit] = ...
%     calculate_rates(model, K, initial_conc, Options, ...
%     mass_balance_ode_conserved,constant_metab_levels,kset_num,cond_num)
%
% kset_num        - number of kset from perturb_Ksets being calculated,
%                   will likely not match with the index of this kset in
%                   the initial ensemble, since the ensemble may have
%                   already been screened
% cond_num        - number of the condition currently being claculated
% TODO


%}

function [complete_time, ode_warn_flag, x, t_return] = ...
    calculate_rates(model, K, initial_conc, Options, ...
    constant_metab_levels, fixed_met_inds, fixed_met_timepoints, fixed_met_concs)
    
% Get fields from model
[n_mets, n_net_rxns] = size(model.S);
S_f_b = model.S_f_b;

% Get volume of the compartment of each species - needed for "equilibrium"
% reactions
gas_vol = 1.5e-3; % 1.5 mL
liq_vol = 30e-6; % 30 uL
gas_species_inds = find(contains(model.mets, "_g"));

% Adjust stoichiometry to capture differences in compartment volumes for
% gas species - default to not adjusting if field isn't present
% Only change stoich of gas species
S_f_b(gas_species_inds, :) = S_f_b(gas_species_inds, :) .* ...
    (liq_vol ./ gas_vol);

t_interval = Options.t_interval;

% preassign ode_warn_flag and norm_slope values; if these values are
% not-reassigned, then the parameter set integration was not
% successful
ode_warn_flag = 0;
norm_slope = 1000;

% track how long integration is running
start_time = tic; 

% Get type of ode solver (e.g., 'ode15s')
ode = Options.ode_solver;
                   
ode_options = odeset();

% Set tolerances (Check if in Options)
if isfield(Options, 'AbsTol') && ~isempty(Options.AbsTol)
    ode_options.AbsTol = Options.AbsTol;
elseif isfield(Options, 'RelTol') && ~isempty(Options.RelTol)
    ode_options.RelTol = Options.RelTol;
else
    ode_options.RelTol = 1e-30;
end

% Set non-negative if optioned
if Options.ode_nonneg
    ode_options.NonNegative = true;
end

% Add Event Function
if Options.use_event_fxn
    if ~Options.event_fxn_ignore_dxdt
        error("Haven't implemented this")
    end

    ode_options.Events = @(t, x) event_function_conserved(...
        t,x,start_time,[],[],[],[],[],[],[],[],Options,[],...
        [],[],[],[],[],[],[],[]);
end

% Get fields
sat_constants = model.lit_values.sat_constants;
S_f_b_trimmed = S_f_b;
init_concs_trimmed = initial_conc(1 : n_mets);

% Remove reg reactions from stoich (just don't want concs
% changing from those reactions when getting dx/dt = S_f_b * v)
if isfield(Options, 'reg_elem_rxns')
    S_f_b_trimmed(:, Options.reg_elem_rxns) = 0;
end

const_met_inds = [model.constant_metabs_list(:); fixed_met_inds(:)];

sat_constants(:, const_met_inds) = [];
S_f_b_trimmed(const_met_inds, :) = [];
init_concs_trimmed(const_met_inds) = [];

% Get indices for additive saturation term-type reactions
% (random-order Michaelis-Menten ternary complex rates)
add_sat_net_rxns = transpose(find(model.rxn_type == 50));

for as = add_sat_net_rxns
    add_sat_elem_inds = ...
        model.rxn_indices(as, 1) : model.rxn_indices(as, 2);

    % Duplicate Km values onto both forward and reverse
    % reactions
    Km_substrate_inds = find(~isnan(...
        sat_constants(add_sat_elem_inds(1), :) ) );

    Km_product_inds = find(~isnan(...
        sat_constants(add_sat_elem_inds(2), :) ) );

    % Make sure there aren't conflicting values that were
    % incorrectly manually entered
    Kms_in_both_dir_inds = intersect(Km_substrate_inds, Km_product_inds);
    % If there are, make sure they're the same value, then
    % ignore product indices
    for i=1:length(Kms_in_both_dir_inds)
        Km_both_concs = sat_constants(add_sat_elem_inds, Kms_in_both_dir_inds(i));
        if diff(Km_both_concs) > 1e-6
            error("Conflicting Km values given for rxn %s",...
                model.rxns{as})
        end
        fprintf("Warning: Km for rxn %s listed in both directions",...
            model.rxns{as})
    end

    % Duplicate values
    % Add product Km values into forward reaction
    sat_constants(add_sat_elem_inds(1), Km_product_inds) = ...
        sat_constants(add_sat_elem_inds(2), Km_product_inds);
    % Add substrate Km values into reverse reaction
    sat_constants(add_sat_elem_inds(2), Km_substrate_inds) = ...
        sat_constants(add_sat_elem_inds(1), Km_substrate_inds);
end


elast_fxn = str2func(strcat('elasticity_coeff_',model.model_name));
ode_options.Jacobian = @(t,x) jacobian_base(...
    t, x, K, sat_constants, S_f_b_trimmed, elast_fxn);

% Set up mex and non-mex mb functions
mb_fxn_std = str2func(strcat('mass_balance_ode_', ...
                             model.model_name));
mb_fxn_mex = str2func(strcat('mass_balance_ode_', ...
                             model.model_name, '_mex'));

if isfield(Options,'use_mex') && Options.use_mex
    % Make mex primary, non-mex alt (for testing)
    mb_fxn = mb_fxn_mex;     

else
    % Make non-mex primary, mex alt (for testing)
    mb_fxn = mb_fxn_std;     
end

[t_return, x_all] = ode(@(t,x) mb_fxn(t, x, K, sat_constants, S_f_b_trimmed),...
    t_interval,init_concs_trimmed,ode_options);

% If it fails, lower AbsTol (maybe RelTol) and try again
if size(x_all, 1) < length(t_interval)
    ode_options.AbsTol = ode_options.AbsTol / 1e20;

    if Options.use_event_fxn

        ode_options.Events = @(t, x) event_function_conserved(...
            t,x,start_time,[],[],[],[],[],[],[],[],Options,[],...
            [],[],[],[],[],[],[],[]);
    end
    [t_return, x_all] = ode(@(t,x) mb_fxn(t, x, K, sat_constants, S_f_b_trimmed),...
        t_interval,init_concs_trimmed,ode_options);

end

x = zeros(length(t_return), n_mets);
nonfixed_mets = setdiff(1:n_mets, const_met_inds);
x(:, nonfixed_mets) = x_all;
% Add constant metabs
x(:, model.constant_metabs_list) = repmat(constant_metab_levels, length(t_return), 1);

% If integration stops before reaching full timespan, 
% Assume it reached a steady-state and just replicate the last row
if size(x, 1) < length(t_interval)
    if size(x, 1) == 1
        ode_warn_flag = -1; % Integration failed
        fprintf("Warning: Integration stopped at initial condition\n")
    else
        ode_warn_flag = 2;
        fprintf("Warning: Integration stopped after initial condition but before full tspan\n")
	end

    new_rows = length(t_interval) - size(x, 1);
    x = [x; repmat(x(end, :), new_rows, 1)];

end

t_return = t_interval;

% Add fixed timecourse metabs
for i = 1:length(fixed_met_inds)
    met_tc_concs = interp1(...
        fixed_met_timepoints(i, :), fixed_met_concs(i, :), t_return);

    x(:, fixed_met_inds(i)) = met_tc_concs;
end

% Stop timer
complete_time = toc(start_time);

% Throw warning if any species are "significantly" negative
if any(x < -1e-6, 'all')
    ode_warn_flag = -2;
    fprintf("Warning: negative concentration in integration\n")
end
    
end
