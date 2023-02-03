% CALCULATE_KSET_FITNESS Calculate fitness for a single kset, using one of
% several fitness functions
%
% [fitness] = calculate_Kset_fitness(model, Experimental_Data, solutions,...
%     pert_num, Options, elem_fluxes)

% TODO

function fitness = calculate_Kset_fitness(model, Experimental_Data, ...
    pert_num, Options, all_metab_concs, calc_timepoints)



    
%% Get experimental data and model predicted data vectors

if iscell(all_metab_concs)
    metab_concs = all_metab_concs{k};
else
    metab_concs = all_metab_concs;
end

 n_timepoints = size(metab_concs,1);

 assert(n_timepoints > 0, strcat("Option `use_metab_timepoints_for_fitness",...
     " cannot be used without field last_metab_concs\n"));

 timecourse = Experimental_Data.metab_timecourse{pert_num};

 if isfield(Options, 'metabs_for_fitness_fxn') && ...
         ~isempty(Options.metabs_for_fitness_fxn)

     [~, timecourse_met_inds] = ...
         ismember(Options.metabs_for_fitness_fxn, timecourse(:, 1));

     timecourse = timecourse(timecourse_met_inds, :);
 end

 n_metabs = size(timecourse,1);

 exp_vals = [];
 pred_vals = [];

 % Set `timecourse_weights` to track weight of each timepoint
 if isfield(Options,'timecourse_fitness_weights')
        timecourse_weights_per_data = [];
        weights_by_time = Options.timecourse_fitness_weights;
        % This field should have an entry for each element of the
        % predicted tspan, not just experimental timepoints
        assert(length(weights_by_time) == length(Options.t_interval),...
            strcat("`Options.timecourse_fitness_weights` should ",...
            "have ",num2str(length(Options.t_interval)),...
            "elements - one for each predicted timepoint") );
 elseif isfield(Options, 'final_timepoint_weight')
     final_weight = Options.final_timepoint_weight;
     assert((0 <= final_weight) && (final_weight <= 1),...
         "Options.final_timepoint_weight must be between 0 and 1\n")
     timecourse_weights_per_data = [];
 end

 % Set fitness multipliers
 met_fitness_multiplier_vals = [];

 for m = 1:n_metabs
     %Get index of this metab
     met_name = timecourse{m,1};
     met_ind = find(strcmpi(model.metabs_and_enzyme_complexes,...
         met_name));
     assert(length(met_ind) == 1,...
         strcat("Didn't find single metabolite: ",timecourse{m,1}));

     exp_timepoints = timecourse{m,2}(1,:);

     tspan_inds = zeros(size(exp_timepoints));

     for t=1:length(exp_timepoints)
         % get the index in t_interval that corresponds with each
         % experimental timepoint (for now, EXACTLY)
         ind = find(Options.t_interval == exp_timepoints(t));
         assert(length(ind) == 1,...
             strcat("No timepoint corresponds with experimental data"));

         tspan_inds(t) = ind;
     end

     exp_metab_concs = timecourse{m,2}(2,:);
     exp_metab_concs = exp_metab_concs(:)';

     % If predicted concs stop before end of experimental, set
     % additional (predicted) values to zero
     uncalc_t_inds = find(~ismember(exp_timepoints, calc_timepoints));
     tspan_inds(uncalc_t_inds) = [];
     pred_metab_concs = metab_concs(tspan_inds,met_ind);
     pred_metab_concs(uncalc_t_inds) = 0;
     pred_metab_concs = pred_metab_concs(:)';

     exp_vals = [exp_vals, exp_metab_concs];
     pred_vals = [pred_vals, pred_metab_concs];

     if isfield(Options,'timecourse_fitness_weights')

         % Concatenate
         timecourse_weights_per_data = ...
             [timecourse_weights_per_data, ...
             weights_by_time(tspan_inds)];

     end

     if isfield(Options, 'metab_fitness_multipliers') && ...
             ~isempty(Options.metab_fitness_multipliers)

         fit_mults = Options.metab_fitness_multipliers;

         % Get index of this metab, if there
         [~, mult_cell_ind] = ismember(met_name, fit_mults(:,1));

         % If not there, just make a one
         if mult_cell_ind == 0

             met_fitness_multiplier_vals = [met_fitness_multiplier_vals,...
                ones(1, length(tspan_inds) )];
         else
             met_fitness_multiplier_vals = [met_fitness_multiplier_vals,...
                repmat(fit_mults{mult_cell_ind, 2}, 1, length(tspan_inds) )];
         end

     else
         % If field doesn't exist, just make a one
         met_fitness_multiplier_vals = [met_fitness_multiplier_vals,...
                ones(1, length(tspan_inds) )];
     end

end

n_refs = length(exp_vals);        

%% Select fitness function and calcualte
% get fitness function, make lowercase
fitness_fxn = lower(Options.fitness_fxn);

switch fitness_fxn

    case 'nrmse'

        fitness = sqrt( sum( (pred_vals - exp_vals) .^ 2) ...
            ./ n_refs ) ./ mean(exp_vals);

    case 'old_abs_nrmse'
        % this is what was used before - not exactly normalized RMSE,
        % but should be proportional to it, not counting that the
        % 'mean' experimental value is the mean of absolute values
        fitness = 1/n_refs*sum(1./CV.*(abs(pred_vals - ...         
            exp_vals)./sum(abs(exp_vals))));

    case 'abs_nrmse'

        fitness = sqrt( sum( (pred_vals - exp_vals) .^ 2) ...
            ./ n_refs ) ./ mean(abs(exp_vals));

    case 'rmse'

        fitness = sqrt( sum( ((pred_vals - exp_vals) .* met_fitness_multiplier_vals) .^ 2) ...
            ./ n_refs );

    case 'exp_error'

        fitness = sum( abs(pred_vals - exp_vals)./ ...
            abs(exp_vals) ) ./ n_refs;

    case 'wgt_timecourse_rmse'

        assert(isfield(Options,'timecourse_fitness_weights'),...
            strcat("`Options` struct must have a [n_timepoints x 1] ", ...
            "field `timecourse_fitness_weights`\n") )

        fitness = sqrt( sum( timecourse_weights_per_data .* ...
            ((pred_vals - exp_vals) .^ 2) ) ./ n_refs );

    otherwise

        % This should throw an error
        error("Specified fitness function is not defined.")

        fprintf('Specified fitness function is not defined. Using RMSE as default.\n')

        fitness = sqrt( sum( (pred_vals - exp_vals) .^ 2) ...
            ./ n_refs );

end
                        

end