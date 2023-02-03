function run_EM_screen(savename_date_version, n_best)

% savename = '_em_screen_test';
% if exist('savename_date_version', 'var')
%     savename = strcat(savename_date_version, savename);
% end

k_ranges = [444,1,40;445,0.100000000000000,2;449,0.100000000000000,1;40,10000,1000000;19,200,4000;685,0.00500000000000000,0.500000000000000;611,0.100000000000000,0.200000000000000;11,100,10000;447,0.0200000000000000,0.600000000000000;448,0.0500000000000000,1;682,0.200000000000000,10;683,10,200;684,0.120000000000000,0.240000000000000;661,0.100000000000000,10;404,0.0500000000000000,5;25,1,20;117,0.100000000000000,5;118,0.100000000000000,5;35,0.100000000000000,2;61,0.100000000000000,10;544,0.0700000000000000,0.140000000000000;547,0.0400000000000000,0.0800000000000000;374,0.400000000000000,0.800000000000000;375,0.500000000000000,1;373,0.350000000000000,0.700000000000000;603,0.100000000000000,0.200000000000000;604,0.150000000000000,0.300000000000000;393,3.50000000000000,7;474,0.0600000000000000,0.120000000000000;483,0.0400000000000000,0.0800000000000000;488,0.0180000000000000,0.0360000000000000;495,0.00900000000000000,0.0180000000000000;500,0.0150000000000000,0.0300000000000000;626,0.0500000000000000,0.100000000000000;677,0.00200000000000000,0.00400000000000000;371,0.00100000000000000,0.100000000000000;394,0.00100000000000000,0.100000000000000;392,0.00100000000000000,0.100000000000000;503,0.00100000000000000,0.100000000000000;508,0.00100000000000000,0.100000000000000;473,0.00100000000000000,0.100000000000000;368,0.180000000000000,0.360000000000000;386,0.00700000000000000,0.0140000000000000;430,0.0300000000000000,0.0600000000000000;482,0.180000000000000,0.360000000000000;486,0.0450000000000000,0.0900000000000000];


run_dynamic_EM(...
    'ecoli_butanol_inhib_fix_model',...
    'Options_struct_loadfile','2021_04_01_ecoli_approx_loadfile',...
    'savename', '2022_04_07_local_2e2_inhib_fix_unscreened_no_mex',...
    'k_range',[0, 1e4],...
    'unknown_sat_constant_range',[.01, 1],...
    'rate_constant_bounds',[0.1, 10],...
    'sat_constant_bounds',[.1, 10],...
    'n_init_ksets',200,...
    'opt_alg','EM',...
    'best_n_ksets',200,...
    'unknown_Ki_range', [5e-1, 5e1],...
    'dataset_for_sampling',1,...
    'manual_k_ranges',k_ranges,...
    'metab_fitness_multipliers',{'succ_c',6; '1btoh_c',20},...
    'init_rng_seed',0,...
    'exp_conditions',[3,4],...
    'keep_all_fits', false,...
    'keep_only_fits',true,...
    'use_mex',false);

run_dynamic_EM(...
    'ecoli_butanol_inhib_fix_model',...
    'Options_struct_loadfile','2021_04_01_ecoli_approx_loadfile',...
    'savename', '2022_04_07_local_2e2_inhib_fix_unscreened_yes_mex',...
    'k_range',[0, 1e4],...
    'unknown_sat_constant_range',[.01, 1],...
    'rate_constant_bounds',[0.1, 10],...
    'sat_constant_bounds',[.1, 10],...
    'n_init_ksets',200,...
    'opt_alg','EM',...
    'best_n_ksets',200,...
    'unknown_Ki_range', [5e-1, 5e1],...
    'dataset_for_sampling',1,...
    'manual_k_ranges',k_ranges,...
    'metab_fitness_multipliers',{'succ_c',6; '1btoh_c',20},...
    'init_rng_seed',0,...
    'exp_conditions',[3,4],...
    'keep_all_fits', false,...
    'keep_only_fits',true,...
    'use_mex',true);

