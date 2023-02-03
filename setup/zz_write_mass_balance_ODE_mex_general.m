% ZZ_WRITE_MASS_BALANCE_ODE_MEX_FILE Use codegen toolbox to write the
% function in mass_balance_ode_conserved into a codegen file.
%
% This'general' script is specifically for mass-balance files that were
% pre-written and include saturation (approx/M-M) kinetics
%
% Arguments to be included:
%   

% Revision history:
%{
2019-03-28: jpm
    Added header
    Add lines to include arguments (str7)

%}

function zz_write_mass_balance_ODE_mex_general(model)

[n_independent_metabs, n_rxns_f_b] = size(model.conserved_model_info.Sr);
n_metabs = length(model.conserved_model_info.metab_index_new_to_old);
n_rxns = size(model.S, 2);
n_dependent_metabs = size(model.conserved_model_info.Lo,1);
model_name = model.model_name;
[n_species, n_rate_consts] = size(model.S_f_b);

% Need to exclude constant or fixed-timecourse species
if isfield(model, 'fixed_tc_metabs_list')
    fixed_tc_species = model.fixed_tc_metabs_list;
else
    fixed_tc_species = [];
end 

if isfield(model, 'constant_metabs_list')
    const_species = model.constant_metabs_list;
else
    const_species = [];
end

% Get number of used species only
unused_species_inds = unique([fixed_tc_species(:); const_species(:)]);
n_unused_species = length(unused_species_inds);
n_used_species = n_species - n_unused_species;

% str1 = strcat('codegen mass_balance_ode_conserved -args {zeros(1,1), zeros(');
% str2 = strcat(num2str(n_independent_metabs),',1),zeros(',num2str(n_rxns_f_b));
% str3 = strcat(',1),zeros(',num2str(n_metabs),',',num2str(n_rxns_f_b));
% str4 = strcat('),zeros(',num2str(n_independent_metabs),',',num2str(n_rxns_f_b), ')');
% str5 = strcat(',zeros(',num2str(n_dependent_metabs),',1)');
% str6 = strcat(',zeros(', num2str(n_dependent_metabs),',',num2str(n_independent_metabs), ')'); % Lo
% str7 = strcat(',cell(',num2str(n_rxns_f_b), ',1)'); % demon_metabs
% str8 = strcat(',cell(',num2str(n_rxns_f_b), ',1)'); % demon_Kms
% str9 = strcat(',zeros(',num2str(n_independent_metabs),',1)'); % ref_metab_concs
% str10 = strcat(',zeros(', num2str(n_rxns), ',1)'); % Keq
% str_end = strcat('} -o mass_balance_ode_',model_name,'_mex');
% 
% string = strcat(str1, str2, str3, str4, str5, str6, str7, str8, str9,...
%     str10, str_end);
% 
% eval(string)

fxn_name = strcat('mass_balance_ode_', model.model_name);

mex_file_name = strcat(fxn_name, '_mex');

codegen(fxn_name, '-args', {...
    zeros(1,1),... % t
    zeros(n_used_species, 1),... % x
    zeros(n_rate_consts, 1),... % k
    zeros(n_rate_consts, n_used_species),... % Km
    zeros(n_used_species, n_rate_consts),... % S_f_b
    });

end