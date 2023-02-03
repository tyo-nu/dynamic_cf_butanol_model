% CALCUALTE_PERTURBED_STEADY_STATES Top-level function to take in model
% struct, IE struct, and Experimental Data struct, and distribute each
% perturbation for calculation
%
% [fitnesses, flux_solutions, complete_times, ode_flags, slope_norms,...
%       elem_fluxes, final_metab_concs, metab_concs, exceptions, ...
%       Ksets_kept, perturbed_stoich_struct, alt_net_flux] = ...
%   calculate_perturbed_steady_states(model, Experimental_Data,
%       Initial_Ensemble, calculate_fitness, Options)

% TODO


function [fvals, all_complete_times, ...
    all_ode_warn_flags, timecourse_metab_concs, all_time_points] = ... 
    calculate_perturbed_steady_states(model, Experimental_Data, ...
    Initial_Ensemble, Options)

% Pull All_K, All_fractions from Initial_Ensemble
All_K = Initial_Ensemble.K_sets;

n_conditions = length(Experimental_Data.flux_rxns);
n_ksets = size(All_K,2);

all_complete_times = cell(n_conditions,1);
all_ode_warn_flags = cell(n_conditions,1);

for y = 1:n_conditions

    yth_enzymes_changed = Experimental_Data.enzymes_changed{y};
    yth_enzyme_levels = Experimental_Data.enzyme_levels{y};

    % Get solutions for this perturbation
    [all_complete_times{y}, ...
        all_ode_warn_flags{y}, yth_metab_concs, time_points] = ... 
    perturb_Ksets(model, Experimental_Data, Initial_Ensemble, ...
        y, yth_enzymes_changed, yth_enzyme_levels, Options);

    for i=1:n_ksets
        timecourse_metab_concs{y,i} = yth_metab_concs{1,i};
        all_time_points{y,i} = time_points{i};
    end

end

fvals = zeros(n_ksets, n_conditions);
all_complete_times = cell2mat(all_complete_times);

%Re-calculate fitness for remaining K-sets
for x = 1:n_ksets
    
    for y = 1:n_conditions

        fvals(x,y) = calculate_Kset_fitness(model,...
            Experimental_Data, y, Options,...
            timecourse_metab_concs{y,x}, all_time_points{y,x});

    end

end


end %end function