% ADJUST_INITIAL_METAB_CONCS Get adjusted initial metabolite concentrations
% 
% inputs:
%   model - struct
%   initial_metab_inds - indices, w.r.t. original (non-conserved) Stoich
%                        matrix
%   initial_metab_concs - concentration vector, same size as
%                         initial_metab_inds
%

% Revision history:
%{
2020-03-01: jpm
    Script created

2020-03-13
    Now pulls involved reactions from S_f_b instead of S and
        enz_enzComplex, 
    Also make option whether to always use passed-in concentrations (vs
        potentially lower concentrations if related to another given metab)
    Also check reactions through model.rxn_indices - can have 
    


TODO: might want to make new approach :
to change all levels associated with conserved totals,
    instead of by looking at fractions/species shared across reactions

%}

function [adjusted_metab_inds, adjusted_metab_concs] = ...
    adjust_init_enzyme_fractions(...
    model,initial_metab_inds,initial_metab_concs)

% option whether to always use passed-in concentrations (vs
%  potentially lower concentrations if related to another given metab)
use_given_levels = 0;

new_frac_inds = [];
new_frac_concs = [];
new_rxns_taken = [];
new_associated_metab = [];
% For each adjusted metab, get rxns it is in
for m = 1:length(initial_metab_inds)

    elem_rxns_involved = find(model.S_f_b(initial_metab_inds(m),:));

    all_enz_frac_inds = [];
    
    %Debuging
    rxn_taken = [];

    % Get list of enzyme fraction indices for each reaction, add to
    % list
    for r = 1:length(elem_rxns_involved)

        enz_inds = find(model.S_f_b(:,elem_rxns_involved(r)));

        all_enz_frac_inds = [all_enz_frac_inds; enz_inds];
        
        rxn_taken = [rxn_taken; repmat(elem_rxns_involved(r),size(enz_inds))];
                
    end
    
    % Don't want to make unique at this point - can lose info about which
    % was lower concentration
% % %     [all_enz_frac_inds, frac_inds_unique,~] = unique(all_enz_frac_inds);
% % %     rxn_taken = rxn_taken(frac_inds_unique);
    
    % Add this list to new_metab_inds/concs
    new_frac_inds = [new_frac_inds; all_enz_frac_inds];

    new_frac_concs = [new_frac_concs; ...
        repmat(initial_metab_concs(m),length(all_enz_frac_inds),1)];
    
    %Debugging
    new_rxns_taken = [new_rxns_taken; rxn_taken];
    
    new_associated_metab = [new_associated_metab; ...
        repmat(initial_metab_inds(m),size(rxn_taken) ) ];

end

% Make sure there are no replicates - if there are, take lower
% value (do this by sorting low to high first, then taking first)

% [B,i] = sort(a) gives B=A(i)
[sorted_frac_concs,sorted_inds] = sort(new_frac_concs,'ascend');
sorted_frac_inds = new_frac_inds(sorted_inds);

% [C, ia, ic] = unique(A) will give C = A(ia), A = C(ic)
[unique_frac_inds, inds_to_keep, inds_new_to_old] = ...
    unique(sorted_frac_inds,'first');

unique_frac_concs = sorted_frac_concs(inds_to_keep);


%%% Since I'm grabbing from S_f_b, this should now include indices of
%%% metabs themselves - don't need to add on

% Make sure all initial metabs are included in this final list
assert(all(ismember(initial_metab_inds, unique_frac_inds)), ...
    'One if the initial metabolites passed in was not found in S_f_b');

% If always using the passed-in levels for the initial metabolites, 
% adjust here
if use_given_levels
    
    % Get indices of initial metabs in the new list and remove
    inds_in_unique_fracs = ...
        find(ismember(unique_frac_inds, initial_metab_inds));
    
    unique_frac_inds(inds_in_unique_fracs) = [];
    unique_frac_concs(inds_in_unique_fracs) = [];

    % Add back on at initial concentrations
    
    adjusted_metab_concs = [initial_metab_concs; unique_frac_concs];
    adjusted_metab_inds = [initial_metab_inds; unique_frac_inds];
else
    
    adjusted_metab_concs = unique_frac_concs;
    adjusted_metab_inds = unique_frac_inds;
    
end

%DEBUGGING
sorted_rxns = new_rxns_taken(sorted_inds);
unique_rxns_taken = sorted_rxns(inds_to_keep);

%For testing, get list of all species names and initial levels
species = cell(length(adjusted_metab_inds),3);
for i=1:length(adjusted_metab_inds)
    
    species{i,1} = adjusted_metab_concs(i);
    species{i,3} = model.metabs_and_enzyme_complexes{adjusted_metab_inds(i)};
    if i > length(initial_metab_inds)
        species{i,2} = unique_rxns_taken(i - length(initial_metab_inds));
    end
    
end
   
end