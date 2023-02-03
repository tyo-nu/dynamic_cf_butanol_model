% CALCULATE_ELEMENTARY_FLUXES Take in model structure and reversibilities,
% and calculate the elementary fluxes of each reaction for the reference
% condition.
%
% [v_ik] = calculate_elementary_fluxes(Network_Data, Rref, Options)

% Revision history:
%{

2019-04-16: jpm
    Added header
    Reformatted for consistency
    Renamed 'Network_Data' to 'model'
    Fixed bugs that assume all irreversible reactions are at end of list

2019-05-15
    Fixed bug where rev_rxn_count wasnt updating - was indexing the wrong
        reversibilities

%}

function [v_ik] = calculate_elementary_fluxes(model, Rref, Options, ...
    ref_flux)
     
% Define Model Characteristics

rxn_indices = model.rxn_indices;
n_tot_rxns = max(...
    [size(model.S,2),length(model.rxns),size(model.rxn_indices,1)]);
n_rxn_steps = ((rxn_indices(:,2) - rxn_indices(:,1))+ 1)./2;
n_rev_rxns = sum(model.rxn_type == 1);
rev_rxns = find(model.rxn_type == 1);
n_elem_rxns = length(model.rxns_f_b);
n_reg_rxns = sum(model.rxn_type == 3);
n_irrev_rxns = sum(model.rxn_type == 2);
irrev_rxns = find(model.rxn_type == 2);
rxn_type = model.rxn_type;   

%set within main Options body
min_flux_threshold = abs(Options.min_flux_threshold); 
% set small rxn flux values to 0 for v_ik calculations
ref_flux(ref_flux > -min_flux_threshold & ref_flux < min_flux_threshold) = 0;                                  

v_ik = zeros(n_elem_rxns,1);     
rev_rxn_count = 1;

%%%%%%%%%%%%%%% Solve for elementary reaction rates, v_i,k %%%%%%%%%%%%%%%%
  % v_i,2j-1 = V_i,ref / (1 - R_i,j^(sign(V_i,ref)))
  % v_i,2j = V_i,ref * R_i,j^(sign(V_i,ref)) / (1 - R_i,j^(sign(V_i,ref)))

for i = 1:n_tot_rxns     
    
    if rxn_type(i) == 1
    
        n_steps = n_rxn_steps(i);
        v_index = rxn_indices(i,1) : rxn_indices(i,2);
        
        % Get indices of reversibilities used in this reaction
        revs_used = rev_rxn_count : (rev_rxn_count + n_steps - 1);
        if length(revs_used) ~= n_steps
            error("Using wrong number of reversibilities")
        end

        % update counter
        rev_rxn_count = rev_rxn_count + n_steps;
        
        % Get sign of this reaction in reference/WT
        sign_wt = sign(ref_flux(i));

        if ref_flux(i) ~= 0                                                        

            for j = 1:n_steps

                R_step = Rref(revs_used(j));
                v_ik(v_index(2*j-1)) = ref_flux(i)/(1 - R_step^(sign_wt));     
                v_ik(v_index(2*j)) = ...                                       
                    ref_flux(i)*R_step^(sign_wt)/(1 - R_step^(sign_wt));

                % if reaction is reversed, switch forward and backward reaction rates
                if v_ik(v_index(2*j-1)) < 0                                    
                    switch_forward = v_ik(v_index(2*j-1));
                    switch_backward = v_ik(v_index(2*j));
                    v_ik(v_index(2*j-1)) = -switch_backward;
                    v_ik(v_index(2*j))= -switch_forward;
                end  

            end

        % if net reaction flux is 0, randomly sample elementary reaction fluxes
        % within range of net reaction flux distribution
        else                                                                    
            for j = 1:n_steps

                v_ik(v_index(2*j-1)) = rand(1,1)*abs(max(ref_flux));               
                % set reverse elementary reaction flux to forward value so net 
                % is equal to 0
                v_ik(v_index(2*j)) = v_ik(v_index(2*j-1));                      

            end       

        end
        
    elseif rxn_type(i) == 2
        
        % set elementary fluxes for transport reactions to net rxn flux value
        v_ik(rxn_indices(i,1)) = ref_flux(i);
        
    elseif rxn_type(i) == 11 % fast (near equilibrium) but kinetic
        % Set so kf/kr = Keq, and min(kf, kr) = 1000 (or something large)
        
        if isfield(Options, 'min_equilibrium_kinetic_rate')
            min_eq_rate = Options.min_equilibrium_kinetic_rate;
        else
            min_eq_rate = 1e5;
        end
        dG_i = model.dG_rxn_phys_unitless(i, 1);
        Keq_i = exp(-1 * dG_i);
        % Check if forward or reverse is smaller, set that to min rate
        if Keq_i < 1
            kf = min_eq_rate;
            kr = kf / Keq_i;
        else
            kr = min_eq_rate;
            kf = Keq_i * kr;
        end
            
        % Assuming substrate concs are 1 (if not, corrected in
        % `adjust_params_for_reference_concs`), rates are just k's
        v_ik(rxn_indices(i, 1)) = kf;
        v_ik(rxn_indices(i, 2)) = kr;
        
    elseif ismember(rxn_type(i), [3,6,7]) % Equivalent to x == 3 || x == 6...

        % Sample forward elementary regulation reaction flux within 
        % range of reference fluxes
        forward_steps = rxn_indices(i,1) : 2 : rxn_indices(i,2);
        reverse_steps = forward_steps + 1;
        v_ik(forward_steps) = rand(n_rxn_steps(i),1) * abs(max(ref_flux));                  
        % Set reverse elementary regulation reaction flux to forward 
        % value so net is equal to 0
        v_ik(reverse_steps) = v_ik(forward_steps);
            
    end
    
    % For type 10 (fast equilibrium), keep v_ik at zero so that K ends up
    % being zero - don't want kinetics for these
    
end


end
