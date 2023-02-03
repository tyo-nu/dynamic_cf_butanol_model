% ZZ_CREATE_ELASTICITY_COEFFICIENT_MATRIX Write cells of strings that can
% be used to write pre-written functions to calculate elasticity
% coefficients, as well as rates/mass balance and parameter-based
% elasticity

% Revision history
%{
2021-10-21: jpm
    History started (previously made)
    Added parameter-based elasticity

%}

function [dv_dx, v, dv_dk] = zz_create_elasticity_coefficient_matrix(model, Options)

% IF RUNNING FROM WRITE_COMPILED_MASS_BAL_ELAST_FXN, THIS NOT THE 
% CONSERVATION-ANALYSIS VERSION - ACTUALLY THE STANDARD S_f_b (trimmed of
% constant species)
% % % S_f_b = model.conserved_model_info.S_f_b_conserved; % Why did i do this
S_f_b_trimmed = model.S_f_b_trimmed; % Use orig order
S_trimmed = model.S_trimmed;

sat_constants = model.lit_values.sat_constants;

% Create vars and counters
dv_dx_row = 1;
dv_dk_row = 1;

dv_dx = {};
dv_dk = {};

v = {};
v_row = 1;

mass_action_net_rxns = transpose(find(ismember(model.rxn_type,...
    [1, 2, 11])));

% Only add reg rxns (6 & 7) if acting on mass-action
reg_rxns = find(ismember(model.rxn_type, [6,7,8]));
reg_on_mass_action_net_rxns = [];
for rr = reg_rxns(:)'
    if any(model.rxn_type(find(model.enzyme_rxn_specificity(:, rr))) == 1)
        reg_on_mass_action_net_rxns = [reg_on_mass_action_net_rxns, rr];
    end
end
mass_action_net_rxns = [mass_action_net_rxns(:)', reg_on_mass_action_net_rxns(:)'];


% Add commenting
dv_dx(dv_dx_row) = {'% Mass action reactions'};
v(v_row) = {'% Mass action reactions'};
dv_dk(dv_dk_row) = {'% Mass action reactions'};

mass_action_elem_rxns = [];
for ma = mass_action_net_rxns
    
        new_ma_elem_rxns = model.rxn_indices(ma, 1) : model.rxn_indices(ma, 2);
        
        mass_action_elem_rxns = [mass_action_elem_rxns, new_ma_elem_rxns];
    
    %% Add symbolic within loop for net rxns (if MA / reg on MA)
    n_elem_ks = length(new_ma_elem_rxns);
    kcat = sym('kcat', [n_elem_ks, 1]);
    n_k_pairs = n_elem_ks / 2;
    assert(floor(n_k_pairs) == n_k_pairs, "This should be a whole number")
    k_stoich = repmat([1; -1], n_k_pairs, 1);
    % Set Keq relationship
    Keq = sym('Keq');
    
    haldane_reln = Keq == prod(kcat .^ k_stoich);
    % Get expression for final k
    last_k_exp = isolate(haldane_reln, kcat(end));
    
    for ma_e_rxn = new_ma_elem_rxns(:)'
        ma_e_ind = find(new_ma_elem_rxns == ma_e_rxn);
        elem_rxn_substrate_inds = find(S_f_b_trimmed(:, ma_e_rxn) < 0);
        elem_rxn_subs_stoich = S_f_b_trimmed(elem_rxn_substrate_inds, ma_e_rxn);
        n_elem_cpds = length(elem_rxn_substrate_inds);
        X = sym('x', [length(elem_rxn_substrate_inds), 1]);
        
        % Get symbolic expression
        v_sym_ma_elem = kcat(ma_e_ind) * prod(X .^ abs(elem_rxn_subs_stoich));
        
        %% For now, I'm only going to use this symbolic expression for dv/dk
        % Replace last k with expression using Keq
        v_sym_ma_elem = subs(v_sym_ma_elem, kcat(end), rhs(last_k_exp));
        
        %% Loop through all k's (except last) and get derivative wrt that k
        for ma_k_ind = 1 : n_elem_ks - 1
            % Get index of this k in the model
            k_model_index = new_ma_elem_rxns(1) + ma_k_ind - 1;
            
            % Take analytic derivative
            dv_dk_sym = diff(v_sym_ma_elem, kcat(ma_k_ind));
            % If this is dv/d(first k), save expression for use in
            % calculating dv/d(last k)
            if ma_k_ind == 1
                dv_dk_first = dv_dk_sym;
            end
            % Make into string
            dv_dk_str = string(dv_dk_sym);
            % Parse symbolic indices into model indices
            dv_dk_str = symb_to_model_x_k_inds(dv_dk_str, Options,...
                elem_rxn_substrate_inds, ma_e_rxn, ...
                [], [], [], [], [], [],...
                new_ma_elem_rxns, ma);
            
            % Finish and add to cell
            dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',...
                ma_e_rxn, k_model_index, dv_dk_str);
            dv_dk_row = dv_dk_row + 1;
            dv_dk(dv_dk_row) = cellstr(dv_dk_new);
            
            
        end
        %% For dv/d(final k), get expression for d(first_k)/d(final_k) and
        % multiply that by dv/d(first k)
        kf_exp = isolate(haldane_reln, kcat(1));
        dk_first_dk_final = diff(rhs(kf_exp), kcat(end));
        % Check if this expression explicitly uses final k (it shouldn't)
        kfinal_tf = has(dk_first_dk_final, kcat(end));
        assert(~kfinal_tf, "Fix me")
        dv_dk_final = dv_dk_first * dk_first_dk_final;
        dv_dk_final_str = string(dv_dk_final);
        dv_dk_final_str = symb_to_model_x_k_inds(dv_dk_final_str, Options,...
                elem_rxn_substrate_inds, ma_e_rxn, ...
                [], [], [], [], [], [],...
                new_ma_elem_rxns, ma);
        % Finish and add to cell
        dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',...
            ma_e_rxn, new_ma_elem_rxns(end), dv_dk_final_str);
        dv_dk_row = dv_dk_row + 1;
        dv_dk(dv_dk_row) = cellstr(dv_dk_new);
        
            
        
        
    end
    
end
mass_action_elem_rxns = mass_action_elem_rxns(:)';



for ma_elem = mass_action_elem_rxns
    
    % Get fluxes
    v_row = v_row + 1;
    v_new = strcat('v(', num2str(ma_elem), ') = k(',...
        num2str(ma_elem), ')');
    
    % Get substrates used in this reaction
    elem_rxn_substrate_inds = find(S_f_b_trimmed(:,ma_elem)<0);
        
    
    %% Loop through metabs to take derivative of each for dv/dx and to add
    % into flux equation / dv/dk
    for met_ind = 1:length(elem_rxn_substrate_inds)
        dv_dx_row = dv_dx_row + 1;
        dv_dx_new = strcat('dv_dx(', num2str(ma_elem), ',', ...
            num2str(elem_rxn_substrate_inds(met_ind)), ') =');      
        
        if length(elem_rxn_substrate_inds) > 1
            other_metab = elem_rxn_substrate_inds(...
                elem_rxn_substrate_inds ~= elem_rxn_substrate_inds(met_ind));
            
            if S_f_b_trimmed(elem_rxn_substrate_inds(met_ind),ma_elem) ~= -1
                rxn_new = strcat(num2str(abs(S_f_b_trimmed(elem_rxn_substrate_inds(met_ind),ma_elem))), ...
                    '*k(',num2str(ma_elem), ')*','x(',num2str(other_metab), ...
                    ')*x(', num2str(elem_rxn_substrate_inds(met_ind)),')^(',...
                    num2str(abs(S_f_b_trimmed(elem_rxn_substrate_inds(met_ind),ma_elem))-1),')');
         
            elseif S_f_b_trimmed(other_metab,ma_elem) ~= -1
                 rxn_new = strcat(' k(',num2str(ma_elem), ')*','x(', ...
                     num2str(other_metab),')^(', ...
                     num2str(abs(S_f_b_trimmed(other_metab,ma_elem))),')');
                
            else
                rxn_new = strcat(' k(', num2str(ma_elem), ')*', ...
                    'x(',num2str(other_metab),')');
                
            end
            
            
        else
            rxn_new = strcat(' k(', num2str(ma_elem), ')');
        end
        
        dv_dx_new = strcat(dv_dx_new, rxn_new, ';');

        dv_dx(dv_dx_row) = cellstr(dv_dx_new);
        
        % Finish fluxes & dv_dk
        if abs(S_f_b_trimmed(elem_rxn_substrate_inds(met_ind), ma_elem)) ~= 1
            v_new = sprintf("%s * x(%i)^%i",...
                v_new, elem_rxn_substrate_inds(met_ind), abs(S_f_b_trimmed(elem_rxn_substrate_inds(met_ind), ma_elem)));

        else
            v_new = sprintf("%s * x(%i)", v_new, elem_rxn_substrate_inds(met_ind));
            
        end
                
    end % Loop over substrates in this elem rxn
    
    % Add in flux row
    v_new = strcat(v_new, ';');
    
    v(v_row) = cellstr(v_new);
    
    
end % Loop over mass action elem rxns

%% Do M-M (additive saturation rates)
add_sat_net_rxns = transpose(find(model.rxn_type == 50));
add_sat_elem_rxns = [];

%% Write v and elasticity for approx rates
% Add commenting
dv_dx_row = dv_dx_row + 1;
v_row = v_row + 1;
dv_dk_row = dv_dk_row + 1;
dv_dx(dv_dx_row) = {'% M-M reactions'};
v(v_row) = {'% M-M reactions'};
dv_dk(dv_dk_row) = {'% M-M reactions'};

% Get fields from model
sat_constants_matrix_locs_trimmed = model.lit_values.sat_constants_matrix_locs_trimmed;

inhib_rxn_inds = find(ismember(model.rxn_type, [6,7]));
enz_specificity = model.enzyme_rxn_specificity;
if any(model.rxn_type == 8)
    error("Need to convert type 8 inhibition (mixed) into componenent comp/uncomp")
end

% Testing - save denominators - should be same within each net rxn
denoms_elem = {};

% Replicate values for forward & reverse rxns
for as = add_sat_net_rxns
    add_sat_elem_inds = ...
        model.rxn_indices(as, 1) : model.rxn_indices(as, 2);

    add_sat_elem_rxns = [add_sat_elem_rxns, add_sat_elem_inds];

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
    for Km_ind = 1:length(Kms_in_both_dir_inds)
        Km_both_concs = sat_constants(add_sat_elem_inds, Kms_in_both_dir_inds(Km_ind));
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

% Loop over elementary reactions
for as_elem = add_sat_elem_rxns
    
    % Get fluxes as symbolic
    subs_inds = find(S_f_b_trimmed(:, as_elem) < 0);
    prod_inds = find(S_f_b_trimmed(:, as_elem) > 0);
    
    % Also get elementary rxn index for reverse of this rxn
    net_rxn_ind = intersect(find(as_elem >= model.rxn_indices(:, 1)), ...
            find(as_elem <= model.rxn_indices(:, 2)));
    % If this is the forward rxn ,reverse is one index after (and vice
    % versa)
    if as_elem == model.rxn_indices(net_rxn_ind, 1)
        rev_elem_ind = as_elem + 1;
    else
        rev_elem_ind = as_elem - 1;
    end
    
    % Get substrates and products for this elementary rxn
    sub_and_prod_inds = find(S_f_b_trimmed(:, as_elem));
    n_subs_and_prods = length(sub_and_prod_inds);
    all_cpds_for_derivatives_inds = sub_and_prod_inds(:);
    
    % Need to check if any inhibitors in this reaction
    if Options.use_reg_rxns       
        % Get net rxn index for this elementary reaction
        
        if as_elem == model.rxn_indices(net_rxn_ind, 1)
            elem_dir = 1;
        else
            elem_dir = 2;
        end
        
        % See if other rxns use this enzyme
        other_rxns_same_enz = setdiff(find(enz_specificity(:, net_rxn_ind)), ...
            net_rxn_ind);
        
        % Find inhib rxn(s) that acts on this full (net) reaction (if any)
        comp_inhib_inds = other_rxns_same_enz(find(...
            model.rxn_type(other_rxns_same_enz) == 6));

        uncomp_inhib_inds = other_rxns_same_enz(find(...
            model.rxn_type(other_rxns_same_enz) == 7));
        
        inhib_rxns = [comp_inhib_inds(:); uncomp_inhib_inds(:)];
        
        % For each competitive rxn
        Ic = sym('Ic', [length(comp_inhib_inds), 1]);
        Kic = sym('Kic', [length(comp_inhib_inds), 1]);
        ci_met_inds = zeros(length(comp_inhib_inds), 1);
        ci_k_inds = zeros(length(comp_inhib_inds), 1);
        
        %% Get metab inds of comp/uncomp inhibitors first
        % need to know which are unique, then just save along with symb 'X'
        for ci = 1:length(comp_inhib_inds)
        	% Get ind of metab inhibiting it
            ci_met_inds(ci) = find(S_trimmed(:, comp_inhib_inds(ci))); % TODO: error
        end
        
        [ci_found_in_X, ci_dup_ind_in_X] = ismember(...
            ci_met_inds, sub_and_prod_inds);
        
        unique_ci_met_inds = ci_met_inds(ci_found_in_X == 0);
        n_unique_ci_mets = sum(ci_found_in_X == 0);
        % Get array of where the new (unique) ci metabs will be in X
        unique_ci_inds_in_X = transpose(...
            n_subs_and_prods + 1 : n_subs_and_prods + n_unique_ci_mets);
        
        ci_inds_in_X = [ci_dup_ind_in_X(find(ci_dup_ind_in_X)); unique_ci_inds_in_X];
            
        % For each uncompetitive rxn
        ui_met_inds = zeros(length(uncomp_inhib_inds), 1);
        for ui = 1:length(uncomp_inhib_inds)
            
            % Get ind of metab inhibiting it
            ui_met_inds(ui) = find(S_trimmed(:, uncomp_inhib_inds(ui)));
        end
        
        [ui_found_in_X, ui_dup_ind_in_X] = ismember(...
            ui_met_inds, [sub_and_prod_inds; unique_ci_met_inds]);
        
        unique_ui_met_inds = ui_met_inds(ui_found_in_X == 0);
        n_unique_ui_mets = sum(ui_found_in_X == 0);
        
        % Get array of where the new (unique) ci metabs will be in X
        unique_ui_inds_in_X =  transpose(n_subs_and_prods + n_unique_ci_mets + 1 : ...
            n_subs_and_prods + n_unique_ci_mets + n_unique_ui_mets);
                
        ui_inds_in_X = [ui_dup_ind_in_X(find(ui_dup_ind_in_X)); unique_ui_inds_in_X];
        
        % Make sure each metab is only listed once - sign of a bug
        if length(ci_met_inds) ~= length(unique(ci_met_inds)) || ...
                length(ui_met_inds) ~= length(unique(ui_met_inds))
            error("Regulation stoich is inconsistent/redundant")
        end
        
        
        all_cpds_for_derivatives_inds = [all_cpds_for_derivatives_inds; ...
            unique_ci_met_inds; unique_ui_met_inds];
        
    end
    
    n_cpds_for_derivatives = length(all_cpds_for_derivatives_inds);
    
    
    rS = S_f_b_trimmed(sub_and_prod_inds, as_elem);
    
    X = sym('x', [n_cpds_for_derivatives, 1]); % X is only 
    Km = sym('Km', [n_subs_and_prods, 1]);
    
       
    v_numerator_sym = prod((X(rS<0)./Km(rS<0)).^abs(rS(rS<0)));
    
    % Get each term for the denominator separately - the [A][B]/(KmA*KmB)
    % term is needed by itself for uncompetitive inhibition
    substrate_concs_over_Kms_term = prod( (X(rS<0) ./ Km(rS<0)) .^ abs(rS(rS<0)) );
    
    % Only add both a*[A]/Kma term AND ([A]/Kma)^a term if there are
    % multiple substrates or multiple products
    % Start with denominator of 1 
    v_denom_part = 1;
    if length(subs_inds) ~= 1 || S_f_b_trimmed(subs_inds, as_elem) ~= -1
        % If only 1 substrate w/ stoich of 1, don't add any terms - only
        % going to use the `substrate_concs_over_Kms_term
        
        % Only add the a*[A]/KmA terms (this statement) if multiple
        % substrates or non-1 stoich
        v_denom_part = v_denom_part + ...
            sum( abs(rS(rS<0)) .* X(rS<0) ./ Km(rS<0) );
    end
    % Then check products - 
    % Always add power term
    products_concs_over_Kms_term = prod( (X(rS>0) ./ Km(rS>0)) .^ abs(rS(rS>0)) );
    % Only add sum term if more than 1 product
    if length(prod_inds) ~= 1 || S_f_b_trimmed(prod_inds, as_elem) ~= 1
        v_denom_part = v_denom_part + ...
            sum( abs(rS(rS>0)) .* X(rS>0) ./ Km(rS>0) );  
    end
    
    % Check for regulation/inhibition on this reaction
    if Options.use_reg_rxns
        
        % See if other rxns use this enzyme
        other_rxns_same_enz = setdiff(find(enz_specificity(:, net_rxn_ind)), ...
            net_rxn_ind);
        
        % Find inhib rxn(s) that acts on this full (net) reaction (if any)
        comp_inhib_inds = other_rxns_same_enz(find(...
            model.rxn_type(other_rxns_same_enz) == 6));

        uncomp_inhib_inds = other_rxns_same_enz(find(...
            model.rxn_type(other_rxns_same_enz) == 7));
        
        inhib_rxns = [comp_inhib_inds(:); uncomp_inhib_inds(:)];
        
        % For each competitive rxn
        Kic = sym('Kic', [length(comp_inhib_inds), 1]);
        ci_k_inds = zeros(length(comp_inhib_inds), 1);
        for ci = 1:length(comp_inhib_inds)
        
            % Add terms - 
            v_denom_part = v_denom_part + X(ci_inds_in_X(ci)) / Kic(ci);

            % Track index of Ki (index of elementary rxn of inhib rxn, not
            % base rxn)
            % For competitive inhibition, should use the same Ki for
            % forward and reverse directions - just grab the index for the
            % forward direction
            ci_k_inds(ci) = model.rxn_indices(comp_inhib_inds(ci), 1);
        end        
            
        % For each uncompetitive rxn
        Kiu = sym('Kiu', [length(uncomp_inhib_inds), 1]);
        ui_k_inds = zeros(length(uncomp_inhib_inds), 1);
        for ui = 1:length(uncomp_inhib_inds)
            
            % Add terms
            substrate_concs_over_Kms_term = substrate_concs_over_Kms_term * ...
                (1 + X(ui_inds_in_X(ui)) / Kiu(ui) );
            products_concs_over_Kms_term = products_concs_over_Kms_term * ...
                (1 + X(ui_inds_in_X(ui)) / Kiu(ui) );
            
            % Track index of Ki (index of elementary rxn of inhib rxn, not
            % base rxn)
            % Standard option will be to use same sampled value for
            % uncompetitive binding to substrate complex as uncompetitive
            % binding to product complex, in which case only use single k
            ui_k_inds(ui) = model.rxn_indices(uncomp_inhib_inds(ui), 1); 
            
        end
                
        % Combine symbolic expression
        v_denom = v_denom_part + ...
            substrate_concs_over_Kms_term + products_concs_over_Kms_term;
        
        v_sym =  v_numerator_sym / v_denom;
        
        % Convert symbolic expression to string and replace Ki indices
        v_str_inds = string(v_sym);
        % Do for comp then uncomp
        for ci = flip(1:length(comp_inhib_inds))
            
            v_str_inds = strrep(v_str_inds, strcat('Kic', num2str(ci)), ...
                strcat('k(', num2str(ci_k_inds(ci)), ')'));

        end
        for ui = flip(1:length(uncomp_inhib_inds))
            
            v_str_inds = strrep(v_str_inds, strcat('Kiu', num2str(ui)), ...
                strcat('k(', num2str(ui_k_inds(ui)), ')'));

        end
        
        % Get fields needed for elasticity
        all_cpds_sym = X;
            
    else
        
        % Just convert symbolic expression to string
        v_denom = v_denom_part + ...
            substrate_concs_over_Kms_term + products_concs_over_Kms_term;
        v_sym =  v_numerator_sym / v_denom;
        v_str_inds = string(v_sym);
        
        % Get fields needed for elasticity
        all_cpds_sym = X;
        all_cpds_for_derivatives_inds = sub_and_prod_inds;
        
        % Make blanks for private fxn
        comp_inhib_inds = [];
        uncomp_inhib_inds = [];
        ci_k_inds = [];
        ui_k_inds = [];
        ci_met_inds = [];
        ui_met_inds = [];
            
    end
    
    denoms_elem{as_elem, 1} = v_denom;
    
    if mod(as_elem,2) == 0 && ~isequal(denoms_elem{as_elem}, denoms_elem{as_elem - 1})
        error("This probably shouldn't be like this")
    end
    
       
    for cpd_ind = flip(1:n_cpds_for_derivatives)
        
        %% For dvdx string
        % Get derivative wrt each cpd INCLUDING inhibitors
        dv_dx_sym = diff(v_sym, all_cpds_sym(cpd_ind));

        % Again change over indices
        dv_dx_str = string(dv_dx_sym);
        dv_dx_str_inds = dv_dx_str;

        
        % Replace substrates and products (x1, x2, x3)
        for xi = flip(1:n_cpds_for_derivatives)
        	% Convert x1, x2, etc into x(ind1), x(ind2),
            dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('x',num2str(xi)), ...
                strcat('x(',num2str(all_cpds_for_derivatives_inds(xi)), ')'));
            % And same for Kms
            dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('Km',num2str(xi)), ...
                strcat('Km(',num2str(as_elem),',', ...
                num2str(all_cpds_for_derivatives_inds(xi)), ')'));
        end
        
        % If using regulatio reactions, convert Ki's to elem k's
        if Options.use_reg_rxns
            for ci = flip(1:length(comp_inhib_inds))

                dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('Kic', num2str(ci)), ...
                    strcat('k(', num2str(ci_k_inds(ci)), ')'));
                dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('Ic', num2str(ci)), ...
                    strcat('x(', num2str(ci_met_inds(ci)), ')'));
            end
            for ui = flip(1:length(uncomp_inhib_inds))

                dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('Kiu', num2str(ui)), ...
                    strcat('k(', num2str(ui_k_inds(ui)), ')'));
                dv_dx_str_inds = strrep(dv_dx_str_inds, strcat('Iu', num2str(ui)), ...
                    strcat('x(', num2str(ui_met_inds(ui)), ')'));
            end

        end     

        % Index dvdx row and add into cell array
        dv_dx_row = dv_dx_row + 1;
        dv_dx_new = sprintf('dv_dx(%i,%i) = k(%i) * (%s);',...
            as_elem, all_cpds_for_derivatives_inds(cpd_ind), as_elem, dv_dx_str_inds);

        dv_dx(dv_dx_row) = cellstr(dv_dx_new);

    end
    
    
    %% Convert symb inds to string (model indices) for v
    for xi = flip(1:n_cpds_for_derivatives)
        
        % Convert x1, x2, etc into x(ind1), x(ind2),
        v_str_inds = strrep(v_str_inds, strcat('x',num2str(xi)), ...
            strcat('x(',num2str(all_cpds_for_derivatives_inds(xi)), ')'));
        % And same for Kms - there shouldn't be any Km's for inhibitors,
        % but it won't hurt to keep this in the same loop
        v_str_inds = strrep(v_str_inds, strcat('Km',num2str(xi)), ...
            strcat('Km(',num2str(as_elem),',', ...
            num2str(all_cpds_for_derivatives_inds(xi)), ')'));
        
    end
    
    % Add in fluxes
    v_row = v_row + 1;
    v_new = sprintf('v(%i) = k(%i) * %s;', ...
        as_elem, as_elem, v_str_inds);
    
    v(v_row) = cellstr(v_new);
    
    %% Get dv_dk
    % Need to make new symbolic rate using kcat values
    kcat = sym('kcat', [2,1]);
    v_with_k_sym = v_sym * kcat(elem_dir);
    
    % Substitute out k_r anywhere in this rate with Keq
    Keq = sym('Keq');
    kcat_stoich = [1; -1];
    Km_forward_stoich = S_trimmed(sub_and_prod_inds, net_rxn_ind);
    kcat_model_inds = model.rxn_indices(net_rxn_ind, :);
    haldane_reln = Keq == prod(kcat .^ kcat_stoich) * prod(Km .^ Km_forward_stoich);
        
    % Get k_r isolated 
    kr_exp = isolate(haldane_reln, kcat(2));
    
    % Substitute out k_r if present
    v_with_k_sym = subs(v_with_k_sym, kcat(2), rhs(kr_exp));
    
    % Get derivative wrt k_f
    dv_dkf_sym = diff(v_with_k_sym, kcat(1));
    dv_dkf_str = string(dv_dkf_sym);
    % replace indices in string
    dv_dkf_str = symb_to_model_x_k_inds(dv_dkf_str, Options,...
        all_cpds_for_derivatives_inds, as_elem, ...
        comp_inhib_inds, uncomp_inhib_inds, ...
        ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds,...
        kcat_model_inds, net_rxn_ind);
    
    % Save into cell
    dv_dk_row = dv_dk_row + 1;
    dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',... 
        as_elem, kcat_model_inds(1), dv_dkf_str);
    dv_dk(dv_dk_row) = cellstr(dv_dk_new);
        
    % Get d(kf)/d(kr)
    kf_exp = isolate(haldane_reln, kcat(1));
    dkf_dkr = diff(rhs(kf_exp), kcat(2));
    
    % Set d(velem)/d(kr) = d(velem)/d(kf) * d(kf)/d(kr)
    dv_dkr_sym = dv_dkf_sym * dkf_dkr;
    dv_dkr_str = string(dv_dkr_sym);
    % replace indices in string
    dv_dkr_str = symb_to_model_x_k_inds(dv_dkr_str, Options,...
        all_cpds_for_derivatives_inds, as_elem, ...
        comp_inhib_inds, uncomp_inhib_inds, ...
        ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds,...
        kcat_model_inds, net_rxn_ind);
    
    % Save into cell
    dv_dk_row = dv_dk_row + 1;
    dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',... 
        as_elem, kcat_model_inds(2), dv_dkr_str);
    dv_dk(dv_dk_row) = cellstr(dv_dk_new);

    
    % Loop through all Km's and Ki's in this elementary reaction - get
    % dv/d(rate_const) later (use same str as v eq)
    for Km_sym_ind = 1:length(Km)
        
        % Get the index of this Km in the full parameter vector
        metab_ind = sub_and_prod_inds(Km_sym_ind);
        sat_const_rxn_hit = find(ismember(sat_constants_matrix_locs_trimmed(:, 1), ...
            [as_elem, rev_elem_ind]));
        sat_const_met_hit = find(sat_constants_matrix_locs_trimmed(:, 2) == metab_ind);
        sat_const_locs_row = intersect(sat_const_rxn_hit, sat_const_met_hit);
        
        Km_param_ind = sat_const_locs_row + size(S_f_b_trimmed, 2);
       
        % Get derivative
        dv_dk_sym = diff(v_with_k_sym, Km(Km_sym_ind));
        % Convert to string
        dv_dk_str = string(dv_dk_sym);
        
        dv_dk_str = symb_to_model_x_k_inds(dv_dk_str, Options,...
            all_cpds_for_derivatives_inds, as_elem, ...
            comp_inhib_inds, uncomp_inhib_inds, ...
            ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds,...
            kcat_model_inds, net_rxn_ind);
        
        dv_dk_row = dv_dk_row + 1;
        dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',... 
            as_elem, Km_param_ind, dv_dk_str);
        dv_dk(dv_dk_row) = cellstr(dv_dk_new);
        
    end
    
    % Check if any Ki's for dv/dk
    if Options.use_reg_rxns
       
        for Kic_sym_ind = 1:length(comp_inhib_inds)
            % Get the index of this Ki in the full parameter vector
            Kic_param_ind = ci_k_inds(Kic_sym_ind);
            
            % take derivative
            dv_dk_sym = diff(v_with_k_sym, Kic(Kic_sym_ind));
            dv_dk_str = string(dv_dk_sym);
            
            % Replace inds
            dv_dk_str = symb_to_model_x_k_inds(dv_dk_str, Options,...
                all_cpds_for_derivatives_inds, as_elem, ...
                comp_inhib_inds, uncomp_inhib_inds, ...
                ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds,...
                kcat_model_inds, net_rxn_ind);
            
            % Save in new row
            dv_dk_row = dv_dk_row + 1;
            dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',...
                as_elem, Kic_param_ind, dv_dk_str);
            dv_dk(dv_dk_row) = cellstr(dv_dk_new);
            
        end
        for Kiu_sym_ind = 1:length(uncomp_inhib_inds)
            % Get the index of this Ki in the full parameter vector
            Kiu_param_ind = ui_k_inds(Kiu_sym_ind);
            
            % take derivative
            dv_dk_sym = diff(v_with_k_sym, Kiu(Kiu_sym_ind));
            dv_dk_str = string(dv_dk_sym);
            
            % Replace inds
            dv_dk_str = symb_to_model_x_k_inds(dv_dk_str, Options,...
                all_cpds_for_derivatives_inds, as_elem, ...
                comp_inhib_inds, uncomp_inhib_inds, ...
                ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds,...
                kcat_model_inds, net_rxn_ind);
            
            % Save in new row
            dv_dk_row = dv_dk_row + 1;
            dv_dk_new = sprintf('dv_dk(%i,%i) = %s;',...
                as_elem, Kiu_param_ind, dv_dk_str);
            dv_dk(dv_dk_row) = cellstr(dv_dk_new);
            
        end
    end
    
    
    
end % Loop over approx elementary rxns

dv_dx = dv_dx';

v = v';

dv_dk = dv_dk';

end

%% Private function to replace symbolic indices in string with model indices
function model_ind_string = symb_to_model_x_k_inds(symb_ind_string, Options,...
        all_cpds_for_derivatives_inds, as_elem, ...
        comp_inhib_inds, uncomp_inhib_inds, ...
        ci_k_inds, ui_k_inds, ci_met_inds, ui_met_inds, ...
        rate_const_model_inds, net_rxn_ind)
    
    model_ind_string = symb_ind_string;

    % Loop through x's, Km's, replace Sym indices with model indices
    for xil = flip(1:length(all_cpds_for_derivatives_inds))
        % Convert x1, x2, etc into x(ind1), x(ind2),
        model_ind_string = strrep(model_ind_string, strcat('x',num2str(xil)), ...
            strcat('x(',num2str(all_cpds_for_derivatives_inds(xil)), ')'));
        % And same for Kms
        model_ind_string = strrep(model_ind_string, strcat('Km',num2str(xil)), ...
            strcat('Km(',num2str(as_elem),',', ...
            num2str(all_cpds_for_derivatives_inds(xil)), ')'));
    end
    if Options.use_reg_rxns
        for cil = flip(1:length(comp_inhib_inds))

            model_ind_string = strrep(model_ind_string, strcat('Kic', num2str(cil)), ...
                strcat('k(', num2str(ci_k_inds(cil)), ')'));
            model_ind_string = strrep(model_ind_string, strcat('Ic', num2str(cil)), ...
                strcat('x(', num2str(ci_met_inds(cil)), ')'));
        end
        for uil = flip(1:length(uncomp_inhib_inds))

            model_ind_string = strrep(model_ind_string, strcat('Kiu', num2str(uil)), ...
                strcat('k(', num2str(ui_k_inds(uil)), ')'));
            model_ind_string = strrep(model_ind_string, strcat('Iu', num2str(uil)), ...
                strcat('x(', num2str(ui_met_inds(uil)), ')'));
        end
    end
    
    % Add in rate constant model indices for rate constants (kcats) if
    % those are present in this expression (will be if doing dv/dk and have
    % substituted a k with an expression for other k's)
    for rc_ind = 1 : length(rate_const_model_inds)
        % Convert kcat1, kcat2, etc into k(ind1), k(ind2), etc
        model_ind_string = strrep(model_ind_string, strcat('kcat',num2str(rc_ind)), ...
            strcat('k(',num2str(rate_const_model_inds(rc_ind)), ')'));
        
    end
    
    % Do same thing if Keq is present (will be in same scenario as rate
    % constants)
    % substituted a k with an expression for other k's)
    if ~isempty(net_rxn_ind)
        % Give a net rxn index to Keq
        model_ind_string = strrep(model_ind_string, strcat('Keq'), ...
            strcat('Keq(',num2str(net_rxn_ind), ')'));
        
    end

end