% FUNCTION POST_RECALC_KSET_FITNESS Recalculate fitness of a kset using
% different fitness function or method
%
% [Perturbation_Results] = ...
%   POST_recalc_kset_fitness(model,Experimental_Data,Perturbation_Results,...
%   varargin)

% Revision history:
%{
2019-09-18: jpm
    Script created
2020-08-07
    Updated to give calculate_Kset_fitness needed metab_concs input
2020-10-12
    Added option to only calculate certain ksets
%} 

function [new_PR, Options] = POST_recalc_kset_fitness(model,Experimental_Data,...
    Perturbation_Results,Options,varargin)

%% Parse varargin
p = inputParser;

checkName = @(x) ischar(x) || isstring(x);

addRequired(p,'model',@isstruct)
addRequired(p,'Experimental_Data',@isstruct)
addRequired(p,'Perturbation_Results',@isstruct)
addRequired(p,'Options',@isstruct)
addParameter(p,'fitness_fxn','',checkName)
addParameter(p,'rename_fitness_field',0,@isnumeric)
addParameter(p, 'ksets', [], @isnumeric)

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,model,Experimental_Data,Perturbation_Results,Options,varargin{:})

fitness_fxn = p.Results.fitness_fxn;
rename_fitness_field = p.Results.rename_fitness_field;
ksets = p.Results.ksets;

% save old fitness fxn
old_fitness_fxn = Options.fitness_fxn;

if ~isempty(fitness_fxn)
    Options.fitness_fxn = fitness_fxn;
end

Experimental_Data = check_fitness_function(model,Experimental_Data,...
    'Options',Options);

n_conds = length(Experimental_Data.flux_rxns);
if isempty(ksets)
    n_ksets = size(Perturbation_Results.fitness_values,1);
    ksets = [1:n_ksets];
    % Copy over all of PR
    new_PR = Perturbation_Results;
else
    n_prev_ksets = size(Perturbation_Results.fitness_values,1);
    n_ksets = length(ksets);
    ksets = ksets(:)';
     %Only copy those ksets
    PR_fields = fieldnames(Perturbation_Results);
    PR_fields = PR_fields(:)';
    for f = PR_fields
        old_field = Perturbation_Results.(f{1});
        % TODO - this could get confused on multiple conditions
        if iscell(old_field) && length(old_field) == n_prev_ksets
            new_PR.(f{1}) = old_field(ksets);
        end
    end
    
end

n_mets_and_enz_comps = size(model.S_f_b,1);



new_fits = zeros(n_ksets, n_conds);

% Get metab conc timecourses if available
if isfield(Perturbation_Results,'last_metab_concs') && ...
        ~isempty(Perturbation_Results.last_metab_concs)
    all_metab_concs = Perturbation_Results.last_metab_concs;
else
    all_metab_concs = repmat({zeros(0,0)},n_conditions,n_ksets);
end

for cond=1:n_conds
    parfor k = 1:n_ksets
       
        
        metab_concs = all_metab_concs{cond, ksets(k)};
        calc_timepoints = Perturbation_Results.timepoints{cond, ksets(k)}
        if isfield(Perturbation_Results, 'solutions')
            net_fluxes = Perturbation_Results.solutions{cond}(:,ksets(k));
        else
            net_fluxes = NaN;
        end
        if isfield(Perturbation_Results, 'elem_fluxes')
            elem_fluxes = Perturbation_Results.elem_fluxes{cond}(:,ksets(k));
        else
            elem_fluxes = NaN;
        end
        new_fits(k,cond) = ...
            calculate_Kset_fitness(model,Experimental_Data,...
                net_fluxes, cond, Options, elem_fluxes, metab_concs, calc_timepoints);
        
    end    
% 
%     metab_concs = all_metab_concs(cond, :);
%     net_fluxes = Perturbation_Results.solutions{cond};
%     elem_fluxes = Perturbation_Results.elem_fluxes{cond};
%     new_fits(:,cond) = ...
%         calculate_Kset_fitness(model,Experimental_Data,...
%             net_fluxes, cond, Options, elem_fluxes, metab_concs);
end

if rename_fitness_field
    new_fits_name = strcat('fitness_values_',Options.fitness_fxn);
    new_PR.(new_fits_name) = new_fits;
else
    old_fits_name = strcat('fitness_values_',old_fitness_fxn);
    new_PR.(old_fits_name) = Perturbation_Results.fitness_values;
    new_PR.fitness_values = new_fits;
    new_PR.avg_fitness_values = mean(new_fits,2);

end



