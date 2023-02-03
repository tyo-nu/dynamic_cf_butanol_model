n_best = 1;


savename_base = '2022_02_24_refactor_opt_test_best1';
loadfile = '2022_02_25_refactor_test_1e1best2';

jobs = cell(n_best, 1);

counter = 1;

while counter <= n_best
        
    try
        
        fprintf("Starting run %i\n\n\n\n\n\n", counter)
        
        warning off

        savename = strcat(savename_base, '_init_', num2str(counter));

        jobs{counter} = run_dynamic_EM(...
        'ecoli_butanol_trim_approx_model',...
        'Options_struct_loadfile','2021_04_01_ecoli_approx_loadfile',...
        'best_ks_loadfile',loadfile,...
        'savename', savename,...
        'k_range',[0, 1e4],...
        'unknown_sat_constant_range',[0.01, 1],...
        'rate_constant_bounds',[0.1, 10],...
        'sat_constant_bounds',[0.1, 10],...
        'SwarmSize',n_best,...
        'opt_alg','patternsearch',...
        'best_n_ksets',[],...
        'unknown_Ki_range', [5e-1, 5e1],...
        'dataset_for_sampling',1,...
        'manual_k_ranges',[444,1,40;445,0.1,2;449,0.1,1;40,10000,1000000;19,200,4000;682,0.1,2;683,1,100;685,0.005,0.5;611,0.1,4;61,3,30],...
        'optim_rate_const_bounds', [0, 1e4],...
        'optim_sat_const_bounds',[1e-3, 1e2],...
        'optim_inhib_const_bounds',[1e-1, 1e2],...
        'best_n_ksets_as_init', n_best,...
        'best_init_kset_inds', counter,...
        'metab_fitness_multipliers',{'succ_c',6;'1btoh_c',20},...
        'exp_conditions', [3,4],...
        'optim_param_bounds', [],...
        'MaxTime',2*60,...
        'use_log_params',true,...
        'UseParallel',true);

        fprintf("Finished run %i\n\n\n\n\n\n", counter)

        counter = counter + 1;

    catch ME

        a=1;
        fprintf("Whoops - something failed\n\n\n\n\n")

    end

end
