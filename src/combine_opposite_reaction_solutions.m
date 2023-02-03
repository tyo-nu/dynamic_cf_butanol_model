% COMBINE_OPPOSITE_REACTION_SOLUTIONS Take in experimental data and
% predicted fluxes for opposite irreversible reactions generated from
% approximate rate forms, and return a single reaction with net rate.
%
% [ ] = combine_opposite_reaction_solutions( )
%
%

% Revision history:
%{

2019-04-17: jpm
    Created script

%}

function [combined_predicted_flux,combined_exp_flux,removed_rxns] = ...
    combine_opposite_reaction_solutions(model,yth_rxns,...
    predicted_flux,exp_flux)
    
% Check if any measured reactions are now approximate and treated as
% two irreversible reactions - if so, combine and treat as one
n_flux_rxns = length(yth_rxns);
irrev_pairs = {};
rxns_to_remove = [];
for i=1:n_flux_rxns

    % for all other reactions
    for j=1:n_flux_rxns

        % if they're not the same reaction but have equal but
        % opposite stoichiometry
        if i<j && isequal(model.S(:,yth_rxns(i)),...
                -1*model.S(:,yth_rxns(j)))

            % save that pair as a list, if it's not there yet
            irrev_pairs = [irrev_pairs; [i,j]];

            % and save index of reaction to remove
            rxns_to_remove = [rxns_to_remove; j];

        end
    end
end

% Combine by making the first flux equal to first minus second, and
%  just removing the second element from ref_fluxes, 
combined_predicted_flux = predicted_flux;
combined_exp_flux = exp_flux;
n_irrev_pairs = length(irrev_pairs);
for i=1:n_irrev_pairs
    combined_predicted_flux(irrev_pairs{i}(1),:) = ...
        predicted_flux(irrev_pairs{i}(1),:) - ...
        predicted_flux(irrev_pairs{i}(2),:);
    % Experimental fluxes were just duplicated - don't need to
    % subtract
    combined_exp_flux(irrev_pairs{i}(1)) = ...
        exp_flux(irrev_pairs{i}(1));
end
% And then remove the second element from each
combined_exp_flux(rxns_to_remove) = [];
combined_predicted_flux(rxns_to_remove,:) = [];

removed_rxns = rxns_to_remove;
        
end