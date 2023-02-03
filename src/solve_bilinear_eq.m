% FUNCTION SOLVE_BILINEAR_EQ Formulate the system of ODEs in the
% elementary rate form as two systems of lienar algebraic equations, and
% solve them iteratively until a steady solution is achieved
%
% vars:
%
%

% Revision history:
%{
2019-07-22: jpm
    Created script
2019-07-25
    bug fixes
2019-07-26
    Trying to speed up some of the loops - 'find' statements at lines
        156 and 158
2019-07-30
    New way of changing order (after CA) to account for constant
        metabolites
%}

function [all_ss_concs,error_diag_code,iter] = solve_bilinear_eq(model,K,...
    initial_conc,all_options,ode_options,denom_metabs,denom_Kms)

% exit criteria checker
exit_yn = 0;

converge_tol = all_options.bilinear_converge_tolerance;
max_iterations = all_options.bilinear_max_iterations;

error_diag_code = [];

iter = 0;
convergence_reached = 0;
last_tries = 0;

S_f_b = model.conserved_model_info.S_f_b_conserved;

n_free_metabs = length(model.mets) - length(model.constant_metabs_list);
% Find new index of all metabolites and fractions
%  i.e., if a species was orginally (before CA) at index i, it is now at
%  all_metab_inds(i)
all_metab_inds = cell2mat(model.conserved_model_info.metab_index_old_to_new);
% Since free metabolties were (pre-CA) listed first, the first
% n_free_metabs of the list of new locations are the new indices of the
% free metabolites
free_metab_inds = all_metab_inds(1:n_free_metabs);

n_metabs_and_fractions = size(S_f_b,1);

% Get the new indices (in initial_conc and S_f_b_conserved) of which 
% species are enzyme fractions - NOT related to change in order
frac_inds = setdiff(1:n_metabs_and_fractions,free_metab_inds);

n_indep_species = size(model.conserved_model_info.Sr,1);

n_fracs = length(frac_inds);

% Split elementary stoich matrix into metab and fraction matrices
% (still otherwise maintaining the same order of rows)
S_metabs = S_f_b(free_metab_inds,:);
S_fracs = S_f_b(frac_inds,:);

% Also define stoich of just substrates and products - needed in loop
S_metabs_prods = S_metabs;
S_metabs_prods(S_metabs < 0) = 0;

S_metabs_subs = S_metabs;
S_metabs_subs(S_metabs > 0) = 0;

% Track what rows of S_fracs are independent species
conserved_frac_inds = find(frac_inds <= n_indep_species);

n_elem_rxns = size(S_f_b,2);
rev_rxn_inds = model.rxn_indices(model.rxn_type == 1,:);
n_rev_rxns = size(rev_rxn_inds,1);
rev_elem_rxns = [];
for i=1:n_rev_rxns
    new_inds = rev_rxn_inds(i,1) : rev_rxn_inds(i,2);
    rev_elem_rxns = ...
        [rev_elem_rxns, new_inds];
end

metab_conc = initial_conc(free_metab_inds);
frac_conc = initial_conc(frac_inds);
% Initialize vector to hold last tries in case of issues
last_10_concs = zeros(n_metabs_and_fractions,10);
% Initialize final output vector to zeros - if fails, this will be default
% and give zero flux for all reactions
all_ss_concs = zeros(n_metabs_and_fractions,1);

% Need number of free enzymes to account for enzyme fraction relationships
n_net_rxns = length(model.rxns);
% Track which rxns use NEW free enzymes (1's)
new_free_enzymes = zeros(n_net_rxns,1);
for rxn=1:n_net_rxns
    % if the first element of the enzyme specificity matrix is along the 
    %  diagonal, then that is a new (unique) enzyme - mark with a 1
    rxns_shared = find(model.enzyme_rxn_specificity(rxn,:));
    if model.rxn_type(rxn) == 1 && rxns_shared(1) == rxn
        new_free_enzymes(rxn) = 1;
    end
end

free_enzymes = find(new_free_enzymes);
n_free_enzymes = length(free_enzymes);

% R is which enzymes and fractions are conserved in their totals - not just
% enz_enzComplex, because only want one row if there are promiscuous
% enzymes
R = zeros(n_free_enzymes,n_fracs);

% enz_comps = model.enz_enzComplex;
% enz_comps(all_metab_inds,:) = enz_comps;
new_to_old_order = ...
    cell2mat(model.conserved_model_info.metab_index_new_to_old);
enz_comps = model.enz_enzComplex(new_to_old_order,:);
for enz = 1:n_free_enzymes
    % Find all reactions that use each free enzyme
    enz_comp_inds = find(enz_comps(:,free_enzymes(enz)));
    free_enz_ind = enz_comp_inds(1);  
    R(enz,:) = double(...
        any(enz_comps(frac_inds, enz_comps(free_enz_ind,:) ~= 0),2));
        
end

% Get conserved enzyme totals for each fraction
enz_totals = R*frac_conc;
% Add that to zeros to get rhs of equation
rhs = [zeros(length(conserved_frac_inds),1); enz_totals];

% Before loop, find the index of the substrate fraction/metab in each reaction
substrate_frac_ind = cell(n_elem_rxns,1);
metab_subs = zeros(n_elem_rxns,1);
for rxn=1:n_elem_rxns
    substrate_frac_ind{rxn} = find(S_fracs(:,rxn) < 0);
    metab_subs_ind = find(S_metabs(:,rxn) < 0);
    if ~isempty(metab_subs_ind)
        metab_subs(rxn) = metab_subs_ind;
    end
end

rxns_involved = cell(n_fracs,1);
for row = 1:n_fracs
    rxns_involved{row} = find(S_fracs(row,:));
end

while ~exit_yn
    % save old metab/fraction values to check against later for convergence
    prev_metab_conc = metab_conc;
    prev_frac_conc = frac_conc;    
    
    %% hold metabolite concentrations constant and update enzyme fractions
    %   Make square matrix X of all terms needed to be multiplied by enzyme
    %   fractions to get dx/dt
    X = zeros(n_fracs);
    for row = 1:n_fracs
                
        for rxn = rxns_involved{row}
           
            % Add the k from this rxn to the column of the enzyme fraction
            % which is the substrate (negative stoich) in this reaction
            % along with metabolites with negative stoich for this reaction
            % (if any)
            % Now indices are found outside of loop
            stoich = S_fracs(row,rxn);
            if metab_subs(rxn) ~= 0
                stoich = stoich * metab_conc(metab_subs(rxn));
            end
            
            X(row,substrate_frac_ind{rxn}) = X(row,substrate_frac_ind{rxn}) + ...
                stoich * K(rxn);
            
        end
        
    end
    
%     *Might* be able to make this even faster through a few linear
%       algebra operations
% % %     k_times_metab_subs = K;
% % %     for rxn = rev_elem_rxns
% % %         %Multiply concentration of substrate metabolite into k
% % %         k_times_metab_subs(rxn) = ...
% % %             k_times_metab_subs(rxn) * metab_conc(S_metabs(:,rxn) < 0);
% % %         
% % %     end
    
    
    
    % Only keep independent rows (from CA)
    X = X(conserved_frac_inds,:);   
    
    % Add on rows for enzyme total level relationships
    X_full = [X; R];
    
%     e = sym('e',[n_fracs, 1]);
%     eq = X_full * e == ...
%         [zeros(length(conserved_frac_inds),1); ones(n_free_enzymes,1)];
%     soln_cell = cell(n_fracs,1);
%     [soln_cell{:}] = solve(eq);
%     solver_frac_conc = zeros(n_fracs,1);
%     for i=1:n_fracs
%         solver_frac_conc(i) = double(soln_cell{i});
%     end    
    
    % This is giving the same answer as using the solver with symbolic eq
    % Solve the system
    frac_conc = X_full \ rhs;
    
    %% hold fractions constant and update metabolites
    %   Use eq. 57 from supplement of Maranas paper
    k_times_frac_subs = K;
    for rxn = rev_elem_rxns
        
        % Every reaction should have *exactly one* enzyme fraction being consumed
%         frac_substrate_inds = find(S_fracs(:,rxn) < 0);
%         
%         % Multiply concentration of substrate enzyme fraction into k
%         k_times_frac_subs(rxn) = ...
%             k_times_frac_subs(rxn) * frac_conc(frac_substrate_inds);

        %%% This is faster
        
        %Multiply concentration of substrate enzyme fraction into k
        k_times_frac_subs(rxn) = ...
            k_times_frac_subs(rxn) * frac_conc(substrate_frac_ind{rxn});
        
    end
    
    % Numerator will be terms generating each metabolite
    numerator = S_metabs_prods * k_times_frac_subs;
    % Denominator will be terms consuming each metabolite - flux will be
    % multiplied by concentration of that metabolite
    denominator = S_metabs_subs * k_times_frac_subs;
    
    metab_conc = -1 .* numerator ./ denominator;
    
%     dxdt = numerator + denominator .* metab_conc; % Debugging
    
    %% Check if exit criteria has been reached
    if max(abs(metab_conc - prev_metab_conc)) <= converge_tol && ...
            max(abs(frac_conc - prev_frac_conc)) <= converge_tol
        
        exit_yn = 1;
        convergence_reached = 1;
        error_diag_code = 0;
        
    end
    
    %% advance loop counter and check if limit reached
    iter = iter + 1;
    
    % After max_iterations, try something different
    if iter > max_iterations
        
        % Save current concentrations and try 10 more - if it doesn't
        % change appreciably (order of magnitude less stringent than the 
        % tolerance) in those 10, call it good enough
        %   Metabs wont be saved in correct order, but just need to see 
        %       if it's changing               
        last_tries = last_tries + 1;
        
        if last_tries <= 10
            last_10_concs(:,last_tries) = [metab_conc; frac_conc];
        end
        
        % After 10, check to see if it's 
        if last_tries > 10
           
            first_try = last_10_concs(:,1);
            diff_from_first = last_10_concs - first_try;
            if max(max(abs(diff_from_first))) <= converge_tol * 10
                
                convergence_reached = 1;
                error_diag_code = 1;
            else
                % If not reached, exit with message saying so
                convergence_reached = 0;
                error_diag_code = 2;
            end
            
            % exit regardless at this point
            exit_yn = 1;
            
        end
                
    end
    
end

% save final metab/fraction concentrations
%  If convergence isn't reached, all_ss_concs is already initialized to
%  zeros, and error message is set above
if convergence_reached == 1
    
    all_ss_concs(free_metab_inds) = metab_conc;
    all_ss_concs(frac_inds) = frac_conc;
    
end
    
% output of all_ss_concs is expected as a row vector - transpose that
all_ss_concs = all_ss_concs(:)' ;

end % End function