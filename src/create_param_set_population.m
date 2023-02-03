%% Fucntion to generate an initial ensemble for EM/Patternsearch/PSO/GA

% TODO - documentation

% Make initial populations
function [Population, warningMsgs] = create_param_set_population(...
    GenomeLength, FitnessFxn, ga_options, ...
    model, Options, rng_seed)
% The input arguments to the function are:
% 
% Genomelength — Number of independent variables for the fitness function
% 
% FitnessFcn — Fitness function
% 
% options — Options
% 
% The function returns Population, the initial population for the genetic algorithm.

pop_size = ga_options.PopulationSize;

Population = zeros(pop_size, GenomeLength);
warningMsgs = cell(pop_size, 1);

% Increase rng seed by init_rng_seed
rng_seed = rng_seed + Options.init_rng_seed;

% for i = 1:pop_size % Testing
parfor i = 1:pop_size
    
    rng('default')
    
    % Option to reproduce results from specific ksets in large run 
    if ~isempty(rng_seed)
        rng(rng_seed, 'twister')
    elseif isfield(Options,'use_spec_rng_seeds') && ~isempty(Options.use_spec_rng_seeds)
        rng(Options.use_spec_rng_seeds(i), 'twister');
    else
        rng(i + Options.init_rng_seed,'twister');
    end
    
    [param_set, warningMsgs{i}] = sample_params_from_ranges(model, Options);
    
    Population(i, :) = param_set(:)';
    
end

end



% Helper function to get single param set for model
function [param_set, warningMsg] = sample_params_from_ranges(model, Options)

% For any type 50 (random M-M - additive) or type 5 (multiple saturation,
% multiplicative) reactions, randomly sample within ranges
sat_rxns = transpose(find(ismember(model.rxn_type, [5,50])));

n_net_rxns = length(model.rxns);
n_Kms = size(model.lit_values.sat_constants_matrix_locs, 1);

% Get sampling ranges from model struct
all_kcat_ranges = model.all_kcat_ranges;
all_Km_ranges = model.all_Km_ranges;
    
sampled_kcats = zeros(length(model.rxns_f_b), 1);
sampled_Kms = NaN(length(model.rxns_f_b), length(model.mets));

% Allocate warning message
warningMsg = [];

% Get indices of excluded species for thermo - really just hydrogen 
excluded_species = ["h_c","h_e"];
excluded_inds = find(ismember(model.mets, excluded_species));

sat_const_locs = model.lit_values.sat_constants_matrix_locs;

for r = sat_rxns
    % Get indices of forward and reverse rxns
    elem_f = model.rxn_indices(r, 1);
    elem_r = model.rxn_indices(r, 2);
    
    % Get indices of substrates and products (net rxn)
    net_subs = find(model.S(:, r) < 0);
    net_prods = find(model.S(:, r) > 0);
    
    % Remove h20 and h+ from these (not included in dG_r/Keq and therefore
    % not included in Haldane relationship)
    net_subs = setdiff(net_subs, excluded_inds);
    net_prods = setdiff(net_prods, excluded_inds);
        
    net_subs_and_prods_inds = [net_subs; net_prods];
    net_subs_and_prods_stoich = model.S(net_subs_and_prods_inds, r);
    % Sample Km's first
    % Km ranges correspond to all substrate/products, NOT reactions, so
    % loop through those
    Km_subs_and_prods = zeros(length(net_subs), 1);
    
    % Testing - make sure all Km inds are used
    sp_inds_used = [];
  
    for sp = 1:length(net_subs_and_prods_inds)
    
        % Get index of this substrate/product in model.mets
        sp_ind = net_subs_and_prods_inds(sp);
        % Add to list of those used (for debugging)
        sp_inds_used = [sp_inds_used; sp_ind];
        
        % Get index of this substrate-elem rxn pair in Km ranges
        % (sat_const_matrix_locs)
        elem_rxn_matches = find(ismember(sat_const_locs(:, 1), [elem_f, elem_r]));
        sp_ind_matches = find(sat_const_locs(:, 2) == sp_ind);
        sp_Km_ranges_ind = intersect(elem_rxn_matches, sp_ind_matches);
        if length(sp_Km_ranges_ind) ~= 1
            error("fix this")
        end
        
        switch Options.sat_const_distribution
            case 'log'
                
                Km_subs_and_prods(sp) = exp(log(all_Km_ranges(sp_Km_ranges_ind,1)) + ...
                    diff(log(all_Km_ranges(sp_Km_ranges_ind,:))) .* rand());
                
            case 'uniform'
                
                Km_subs_and_prods(sp) = all_Km_ranges(sp_Km_ranges_ind,1) + ...
                    diff(all_Km_ranges(sp_Km_ranges_ind,:)) .* rand();

            case 'lognorm'
                
                min_max = all_Km_ranges(sp_Km_ranges_ind, :);

                % Make mu the geometric mean (log mean)
                mu = mean(log(min_max));
                % Just make the "edges" of the bins within 1 std
                sigma = diff(log(min_max)) / 2;
                
                Km_subs_and_prods(sp) = lognrnd(mu, sigma);                
                
            otherwise
                error("Distribution is not valid value")
        end
        
    end
    % Pull apart substrates and products
    Km_subs = Km_subs_and_prods(1:length(net_subs));
    Km_prods = Km_subs_and_prods(length(net_subs)+1 : end);
    % Get total term
    Km_prods_over_subs = prod(Km_subs_and_prods(:) .^ net_subs_and_prods_stoich(:));
    
    %% Get Keq value
    % Sample Keq from range
    dG_rxn_mean = model.dG_rxn_phys_unitless(r, 1);
    dG_rxn_err = model.dG_rxn_phys_unitless(r, 2);
    % Going to do strictly within ranges instead of normal distribution -
    % errors are acually 95% CIs
    dG_min = dG_rxn_mean - dG_rxn_err;
    dG_max = dG_rxn_mean + dG_rxn_err;
    dG_sampled = dG_min + rand * (dG_max - dG_min);
    Keq_sampled = exp(-1.*dG_sampled);

    
    % Sample *larger* kcat - will depend on Keq and sampled Km's
    kcat_f_over_r = Keq_sampled / Km_prods_over_subs;
       
    if kcat_f_over_r >= 1
        larger_kcat_ind = elem_f;
        smaller_kcat_ind = elem_r;
        calc_dir = "reverse";
    else
        larger_kcat_ind = elem_r;
        smaller_kcat_ind = elem_f;
        calc_dir = "forward";
        % Take reciprocal of ratio so we can pretend like forward is largest
        kcat_f_over_r = 1 / kcat_f_over_r;
    end
        
    % Get kcat range, set to min values
    rxn_kcat_range = all_kcat_ranges(larger_kcat_ind,:);
    rxn_kcat_range(isnan(rxn_kcat_range)) = 0;
    rxn_kcat_range(rxn_kcat_range < Options.min_abs_kcat) = Options.min_abs_kcat;
    
    % Sample rate constant value - I'm being inconsistent/meaningless with
    % rate_constant vs kcat vs vmax here, so I'm sorry 
    switch Options.rate_const_distribution
        
        case 'log'
            larger_kcat_val = exp( log(rxn_kcat_range(1)) + ...
                diff(log(rxn_kcat_range)) * rand() );
        case 'uniform'
            larger_kcat_val = rxn_kcat_range(1) + ...
                diff(rxn_kcat_range) * rand();
        otherwise
            error("Distribution is not valid value")
    end
    
    % Use Haldane to back out reverse 
    % Keq = (kcat_f / kcat_r) / (Km_products / Km_substrates) = 
    %       (kcat_f / kcat_r) / Km_tot
    smaller_kcat_val = larger_kcat_val / kcat_f_over_r;
    
    % If smaller kcat value from Keq & Haldane isn't within "lit value"
    % ranges, add a warning message
    if smaller_kcat_val < all_kcat_ranges(smaller_kcat_ind, 1) 
        new_warning = sprintf(...
            "Sampled rate constant in the %s direction for rxn %i is less than lb of %0.2e\n",...
            calc_dir, r, all_kcat_ranges(smaller_kcat_ind, 1) );
        warningMsg = [warningMsg; new_warning];        
    elseif smaller_kcat_val > all_kcat_ranges(smaller_kcat_ind, 2) 
        new_warning = sprintf(...
            "Sampled rate constant in the %s direction for rxn %i is greater than ub of %0.2e\n",...
            calc_dir, r, all_kcat_ranges(smaller_kcat_ind, 2) );
        warningMsg = [warningMsg; new_warning];
    end
    
    % Save kcat values
    sampled_kcats([larger_kcat_ind, smaller_kcat_ind]) = ...
        [larger_kcat_val, smaller_kcat_val];
    
    
    % For type 4 reaction, save Km's in sat_constants for both directions
    if model.rxn_type(r) == 4
    
        sampled_Kms(elem_f, [net_subs; net_prods]) = [Km_subs; Km_prods];
        sampled_Kms(elem_r, [net_subs; net_prods]) = [Km_subs; Km_prods];
        
    else
        % For type 5 reaction, only save Km's for substrates (negative stoich)
        % for the relative forward reaction (and products (in the net rxn) for
        % the reverse elementary reaction)
        sampled_Kms(elem_f, net_subs) = Km_subs;
        sampled_Kms(elem_r, net_prods) = Km_prods;    
    
    end
    
end

%% Get fast_eq or other rxns
fast_eq_net_rxns = transpose(find(model.rxn_type == 11));
fast_eq_elem_rxns = [];
for fe = fast_eq_net_rxns
    fast_eq_elem_rxns = [fast_eq_elem_rxns, ...
        model.rxn_indices(fe, 1) : model.rxn_indices(fe, 2)];
end
% sample from known Keq 
% Get from elementary fluxes - set blanks for other rxns
Rref_blank = zeros(length(model.rxns_f_b), 1);
ref_net_flux_blank = zeros(length(model.rxns), 1);
[v_ik] = calculate_elementary_fluxes(model, Rref_blank, Options, ref_net_flux_blank);

% Assign into "kcats"
sampled_kcats(fast_eq_elem_rxns) = v_ik(fast_eq_elem_rxns);

%% Set regulation reaction ks
% Competitive inhibition reactions - only sample a single value (use
% forward direction index), which gets used in both directions
comp_net_rxns = find(model.rxn_type == 6)';
comp_elem_rxns = model.rxn_indices(comp_net_rxns, 1)';

% Uncompetitive inhibition reactions
uncomp_net_rxns = find(model.rxn_type == 7)';
% Option to add uncompetitive inhibition to some of the reverse rates
uncomp_elem_rxns = model.rxn_indices(uncomp_net_rxns, 1)';

all_sampling_reg_elem_rxns = [comp_elem_rxns(:)', uncomp_elem_rxns(:)'];

% Also just get indices of all regulation elementary reactions - set
% default values to 1e6 (needs to be positive to avoid divide-by-zero)
all_reg_net_rxns = [comp_net_rxns, uncomp_net_rxns];
all_reg_elem_rxns = reshape(model.rxn_indices(...
    all_reg_net_rxns, :), 2*length(all_reg_net_rxns), 1);
sampled_kcats(all_reg_elem_rxns) = 1e3;

for rr = all_sampling_reg_elem_rxns
    % Decide on sampling distribution
    switch Options.sat_const_distribution
        case 'log'

            sampled_kcats(rr) = exp(log(all_kcat_ranges(rr,1)) + ...
                diff(log(all_kcat_ranges(rr,:))) .* rand());

        case 'uniform'

            sampled_kcats(rr) = all_kcat_ranges(rr,1) + ...
                diff(all_kcat_ranges(rr,:)) .* rand();

        case 'lognorm'

            min_max = all_kcat_ranges(rr, :);
            % Make mu the geometric mean (log mean)
            mu = mean(log(min_max));
            % Just make the "edges" of the bins within 1 std
            sigma = diff(log(min_max)) / 2;

            sampled_kcats(rr) = lognrnd(mu, sigma);                

        otherwise
            error("Distribution is not valid value")
    end
end

%% Other rxn types
% Haven't implemented other rxn types yet (elementary-rate mass action (1)
% or irreversible non-enzymatic (2) )
if any(ismember(model.rxn_type, [1, 2]))
    error("TODO: implement other rxns in creation fcn")
end

% flatten out saturation constants into a vector along rows first
flattened_Kms = reshape(sampled_Kms', [], 1);
used_Kms = flattened_Kms(~isnan(flattened_Kms));
% Combine params to return

param_set = [sampled_kcats; used_Kms];

% If using log-transformed vars & 'nonlcon', account for slack vars
if isfield(Options,'use_log_params') && Options.use_log_params && ...
        isfield(Options,'use_nonlincon') && Options.use_nonlincon
    % Slack variables should be Keq values
    n_slackvars = n_net_rxns;
    slackvar_inds = length(param_set) + 1 : length(param_set) + n_slackvars;
    
    % Sample in log uniform dist
    slackvar_samples = exp(log(Options.lb(slackvar_inds)) + ...
        (log(Options.ub(slackvar_inds)) - log(Options.lb(slackvar_inds))) .* ...
        rand(n_slackvars, 1) );
    
    % Add onto param vector
    param_set = [param_set; slackvar_samples];
    
end

% If optioned to set a specific value
if isfield(Options, 'manual_k_values') && ~isempty(Options.manual_k_values)
    assert(size(Options.manual_k_values, 2) == 2)
    param_set(Options.manual_k_values(:, 1)) = Options.manual_k_values(:, 2);
end

end

