function rerun_final_ensemble(savename_date_version, n_best)

savename = '_em_rerun_test';
if exist('savename_date_version', 'var')
    savename = strcat(savename_date_version, savename);
end

EM_pert_output = run_dynamic_EM(...
'ecoli_butanol_inhib_fix_model',...
'Options_struct_loadfile','2021_04_01_ecoli_approx_loadfile',...
'best_ks_loadfile','2021_11_18_Quest_ps_init_20210927opt_inhib-bux-fix_250cpu_4hr_best112',...
'savename', savename,...
'k_range',[0, 1e4],...
'unknown_sat_constant_range',[0.01, 1],...
'rate_constant_bounds',[0.1, 10],...
'sat_constant_bounds',[0.1, 10],...
'n_init_ksets',n_best,... %%
'opt_alg','EM',...
'best_n_ksets',[],...
'unknown_Ki_range', [5e-1, 5e1],...
'dataset_for_sampling',1,...
'optim_rate_const_bounds', [0, 1e4],...
'optim_sat_const_bounds',[1e-3, 1e2],...
'optim_inhib_const_bounds',[1e-1, 1e2],...
'best_n_ksets_as_init',n_best,... %%
'best_init_kset_inds',[],...
'metab_fitness_multipliers',{'succ_c',6;'1btoh_c',20},...
'exp_conditions',[3,4],...
'ode_nonneg', false,...
'use_mex',false);