% Get crossover for GA
function xoverKids = ga_crossover_fxn(parents, ga_options, nvars, FitnessFcn,...
    thisScore, thisPopulation, model, xover_ratio, Options)

% Function to take parents and return child for genetic algorithm
% 
% The arguments to the function are
% 
% parents — Row vector of parents chosen by the selection function
% 
% options — options
% 
% nvars — Number of variables
% 
% FitnessFcn — Fitness function
% 
% unused — Placeholder not used
% 
% thisPopulation — Matrix representing the current population. 
%   The number of rows of the matrix is PopulationSize and the number of columns is nvars.

% Taking some of this from default crossover<x> fcns - mostly crossoverheuristic.m
nKids = length(parents)/2;
% Allocate children
xoverKids = zeros(nKids, nvars);

% Get number of rxn "groups" that are constrained by thermo (here, net rxns)
n_net_rxns = length(model.rxns);

% For each child, use the next two parents
index = 1;
for i = 1:nKids
    % get parents
    parent1 = thisPopulation(parents(index),:);
    score1 = thisScore(parents(index));
    index = index + 1;
    parent2 = thisPopulation(parents(index),:);
    score2 = thisScore(parents(index));
    index = index + 1;

    % Get % of reactions to keep for each parent
    % If no ratio is given, use ratio of scores
    if isempty(xover_ratio)
        parent1_ratio = score2 / (score1 + score2);
    else
        % Xover_ratio should be amount to which better parent is biased
        % So xover_ratio of 1 is 50/50 split, xover_ratio of 9 is 90/10 split
        
        % If parent 1 is better, give it more genes (xover_ratio > 1)
        if score1 < score2 
            parent1_ratio = xover_ratio / (1 + xover_ratio);
        % If nearly the same, give same ratio
        elseif abs(score1 - score2) < 1e-3
            parent1_ratio = 0.5;
        else            
            parent1_ratio = 1 - xover_ratio / (1 + xover_ratio);
        end
    end
    
    % Make sure scores weren't NaNs
    if isnan(score2) && isnan(score1)
        parent1_ratio = 0.50;
    elseif isnan(score1)
        parent1_ratio = 1 - xover_ratio / (1 + xover_ratio);
    elseif isnan(score2)
        parent1_ratio = xover_ratio / (1 + xover_ratio);
    end
    
    % Assign individually by reaction type (so each type has same % of 1v2)
    approx_rxn_inds = find(ismember(model.rxn_type, [5,50]));
    n_approx_rxns = length(approx_rxn_inds);
    n_approx_1 = round(n_approx_rxns * parent1_ratio);
    n_approx_2 = n_approx_rxns - n_approx_1;
    
    mass_act_rxn_inds = setdiff(1:n_net_rxns, approx_rxn_inds);
    n_mass_act_rxns = length(mass_act_rxn_inds);
    n_mass_act_1 = round(n_mass_act_rxns * parent1_ratio);
    n_mass_act_2 = n_mass_act_rxns - n_mass_act_1;
    
    % Set base order 
    approx_base_order = zeros(n_approx_rxns, 1);
    mass_act_base_order = zeros(n_mass_act_rxns, 1);
    approx_base_order(1:n_approx_1) = 1;
    approx_base_order(n_approx_1 + 1 : end) = 2;
    mass_act_base_order(1:n_mass_act_1) = 1;
    mass_act_base_order(n_mass_act_1 + 1 : end) = 2;
    
    % Shuffle up to decide which reactions come from which parent
    shuffled_approx_parents = approx_base_order(randperm(n_approx_rxns));
    shuffled_mass_act_parents = mass_act_base_order(randperm(n_mass_act_rxns));
    % Combine to get all new parents
    net_rxn_parents = zeros(n_net_rxns, 1);
    net_rxn_parents(approx_rxn_inds) = shuffled_approx_parents;
    net_rxn_parents(mass_act_rxn_inds) = shuffled_mass_act_parents;

    % Get indices of kcats, Kms (if present), and init concs (if present)
    n_kcats = length(model.rxns_f_b);
    rate_const_inds = 1 : n_kcats;
    n_Kms = 0;
    if isfield(model,'lit_values') && isfield(model.lit_values, 'sat_constants')
        Km_matrix = model.lit_values.sat_constants;
        Km_locs = model.lit_values.sat_constants_matrix_locs;
        n_Kms = size(Km_locs, 1);
    end
    Km_inds = n_kcats + 1 : n_kcats + n_Kms;
    if Options.vary_init_concs
        n_concs = size(model.S_f_b, 1);
        init_conc_inds = (n_kcats + n_Kms + 1) : (n_kcats + n_Kms + n_concs);
    else
        init_conc_inds = [];
    end
    n_init_concs = length(init_conc_inds);
    
    %% Do approx rate rxns    
    % Track # of Kms for saturation-type rxns found
    n_Kms_found = 0;
    % For each reaction
    for r = 1:n_net_rxns
    
        % Get indices of kcats/Kms used for this reaction
        rxn_kcat_inds = model.rxn_indices(r, 1) : model.rxn_indices(r, 2);
        
        rxn_Km_inds = [];
        % Get sat constant params
        if ismember(model.rxn_type(r), [4,5,50])
            elem_rxns = model.rxn_indices(r, 1) : model.rxn_indices(r, 2);
            elem_rxns = elem_rxns(:)';

            % For saturation constants, just need to know how many are in
            % this reaction, and keep a counter
            n_rxn_Kms = length(find(ismember(Km_locs(:, 1), elem_rxns)));
            
            % Set those inds
            rxn_Km_inds = (n_Kms_found + 1) : (n_Kms_found + n_rxn_Kms);
            % Update counter
            n_Kms_found = n_Kms_found + n_rxn_Kms;
        end
        % Add those params to child
        if net_rxn_parents(r) == 1
            rxn_parent = parent1;
        else
            rxn_parent = parent2;
        end
        xoverKids(i, rxn_kcat_inds) = rxn_parent(rxn_kcat_inds);
        xoverKids(i, rxn_Km_inds) = rxn_parent(rxn_Km_inds);
    end

    % Rearrange slack vars if present - keep each slack var (Keq)
    %   with its corresponding net rxn
    if isfield(Options,'use_log_params') && Options.use_log_params && ...
            isfield(Options,'use_nonlincon') && Options.use_nonlincon
        n_slack_vars = n_net_rxns;
        n_actual_params = nvars - n_slack_vars;
        
        for sind = 1:n_slack_vars
            slack_var_ind = n_actual_params + sind;
            if net_rxn_parents(sind) == 1
                slackvar_parent = parent1;
            else
                slackvar_parent = parent2;
            end
            xoverKids(i, slack_var_ind) = slackvar_parent(slack_var_ind);
        end
    end
        
    % For init_concs, assume all are separate
    n_concs_1 = round(n_init_concs * parent1_ratio);
    n_concs_2 = n_init_concs - n_concs_1;
    
    base_order = [repmat(1, n_concs_1, 1); repmat(2, n_concs_2, 1)];
    conc_parents = base_order(randperm(n_init_concs));
    for c = 1 : n_init_concs
        if conc_parents(c) == 1
            conc_parent = parent1;
        else
            conc_parent = parent2;
        end
        % Assign concs
        xoverKids(i, init_conc_inds(c)) = conc_parent(init_conc_inds(c));
        
    end
    
end
    
end