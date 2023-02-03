% GET_FRACTIONS_PER_ENZYME Gets a version of the model.enz_enzComplex that
% is specific by enzyme instead of by reaction, so that each rxn (column)
% with the same enzyme will be identical and include a 1 for all fractions
% associated with that enzyme, not just those in that reaction
% Also returns a bool vector of whether this reaction was the first to use
% this enzyme
% Also also passes back a matrix to multiply 

% Revision history
%{
2020-09-28: jpm
    Script created

%}

function [fractions_per_enzyme, first_rxn_yn, normalized_enz_fractions, ...
    free_enz_ind] = ...
    get_fractions_per_enzyme(model, absolute_levels);

complexes = model.enz_enzComplex;

n_rxns = size(complexes, 2);
% Set up outputs
fractions_per_enzyme = zeros(size(complexes));
first_rxn_yn = zeros(n_rxns, 1);
free_enz_ind = zeros(n_rxns, 1);
if exist('absolute_levels','var')
    normalized_enz_fractions = zeros(size(absolute_levels));
else
    normalized_enz_fractions = NaN;
end

% Get rxns with enzyme
rev_rxns = find(model.rxn_type == 1);
rev_rxns = rev_rxns(:)';

for r = rev_rxns
    
    % Get indices of all fractions in this reaction
    fraction_inds = find(complexes(: ,r));
    
    % Get index of free enzyme
    free_enz_ind(r) = fraction_inds(1);
    
    % See if any other reactions use this enzyme
    enz_rxn_inds = find(complexes(free_enz_ind(r), :));
    
    % Check if this is first
    if r == enz_rxn_inds(1)
        first_rxn_yn(r) = 1;
    end
    
    % Get a list of all fractions for this enzyme
    all_enz_fractions = [];
    for rxn = enz_rxn_inds
        all_enz_fractions = [all_enz_fractions; find(complexes(:, rxn))];
    end
    all_enz_fractions = unique(all_enz_fractions);
        
    fractions_per_enzyme(all_enz_fractions, r) = 1;
    
    % If given, convert absolute fractions (using all levels includign
    % metabs to have consistent dimensions) into fractions that sum to 1
    if exist('absolute_levels','var')
       
        % Get totals for all complexes of this enzyme
        enz_total = sum(absolute_levels(all_enz_fractions));
        
        if first_rxn_yn(r)
            
            normalized_enz_fractions(all_enz_fractions) = ...
                absolute_levels(all_enz_fractions) ./ enz_total;

        end
        
    end
    
end


end