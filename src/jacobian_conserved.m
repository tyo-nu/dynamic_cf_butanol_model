% JACOBIAN Calculate and return jacobian matrix (df/dx, where
% f(x) = dx/dt) at a given metabolite concentration.
% Takes in independent metabolite concentration vector xi, kinetic 
% paramters k, reduced stoichiometric matrix Sr, full link matrix L, 
% conserved metabolite totals T, dependent link matrix Lo, and elasticity
% coefficient function handle, and returns the sparse jacobian matrix df/dx

% Revision history:
%{
2019-03-13: jpm
    added header and reformatted for style guidelines
    Changed name of function to jacobian_conserved, so as to not confuse 
    with built-in function 'jacobian' included in Matlab
2021-01-27
    Option to use elast_coeff_general - now evaluates a string passed in
2021-01-28
    Now takes in fixed_met info, adjusts k's
%}

function df_dx = jacobian_conserved(t,xi,k,Sr,L,T,Lo,...
    elasticity_coeff_w_inputs, subs_per_rxn, ...
    fixed_met_subs_rxns, fixed_met_timepoints, fixed_met_concs)

% Calculate dependent metabolite concentrations
xd = T + Lo*xi;

% Combine independent and dependent metabolites into a single vector
x = [xi; xd];                                                           

% Update k with any fixed timecourse metabolites
for i = length(fixed_met_subs_rxns)
    k(fixed_met_subs_rxns{i}) = k(fixed_met_subs_rxns{i}) .* ...
        interp1(fixed_met_timepoints(i, :), fixed_met_concs(i, :), t);
end
% % % for i = length(fixed_met_subs_rxns)
% % %     x(fixed_met_subs_rxns(i)) = interp1(...
% % %         fixed_met_timepoints(i, :), fixed_met_concs(i, :), t);
% % % end

% Calculate elasticity matrix needed to calculate jacobian of reduced model
% fxn handle `elasticity_coeff_w_inputs` will include inputs
e_ij = eval(elasticity_coeff_w_inputs);

% Calculate Jacobian matrix
df_dx = Sr * e_ij * L;

% Convert Jacobian to sparse matrix for memory and computation time 
% improvements
df_dx = sparse(df_dx);                                                  

end
