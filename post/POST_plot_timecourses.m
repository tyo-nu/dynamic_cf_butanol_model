% FUNCTION POST_PLOT_TIMECOURSES Plot experimental timecourses of
% metabolite (or enzyme and enzyme complexes) concentrations along with
% options for experimental data
%

% Revision history:
%{
2020-08-07: jpm
    Documentation created
    Added varargin to calculate fitness and show on screen
2020-08-10
    Fixed double axis and legend issue
2020-08-19
    Added features to group subplot by metabolite, show confidence
        intervals instead of individual ksets
    
%}

function POST_plot_timecourses(model,Experimental_Data,...
    Perturbation_Results,Options,varargin)

p = inputParser;

checkName = @(x) ischar(x) || isstring(x) || iscell(x);
checkNameOrEmpty = @(x) ischar(x) || isstring(x) || iscell(x) || isempty(x);
checkNameOrCell = @(x) ischar(x) || isstring(x) || ...
    (iscell(x) && ( ischar(x{1}) || isstring(x{1}) ) );
checkTspan = @(x) length(x) == 2;
default_metabs = {'glc__D_e','succ_c','lac__D_c','ac_c','etoh_c','1butanol','mal__L_c','cit_c','pyr_c'};
default_metab_names = ["glc__D_e","Glucose";"succ_c","Succinate";"lac__D_c","Lactate";"ac_c","Acetate";"etoh_c","Ethanol";"1btoh_c","1-Butanol";"mal__L_c","Malate";"cit_c","Citrate";"pyr_c","Pyruvate";"accoa_c","Acetyl-CoA";"atp_c","ATP";"nadh_c","NADH";"coa_c","CoA"];

default_subplot_titles = ["Glucose","Succinate","Lactate","Acetate","Ethanol","1-Butanol","Malate","Citrate","Pyruvate"];
default_ratios_w_nadp = {'nadh_c',{'nadh_c', 'nad_c'};...
                    'nadph_c',{'nadph_c', 'nadp_c'};...
                    'atp_c',{'atp_c', 'adp_c', 'amp_c'};...
                    'coa_c',{'coa_c','accoa_c','succoa_c','aacoa_c','3hbd_c','ctncoa_c','btcoa_c'} };
default_ratios_no_nadp = {'nadh_c',{'nadh_c', 'nad_c'};...
                    'atp_c',{'atp_c', 'adp_c', 'amp_c'};...
                    'coa_c',{'coa_c','accoa_c','succoa_c','aacoa_c','3hbd_c','ctncoa_c','btcoa_c'} };                
default_ratios = default_ratios_no_nadp;                
                
checkBool = @(x) islogical(x) || isnumeric(x);
isNatural = @(x) isnumeric(x) && mod(x,1) == 0 && x > 0;

addRequired(p, 'model',@isstruct)
addRequired(p, 'Experimental_Data',@isstruct)
addRequired(p, 'Perturbation_Results',@isstruct)
addRequired(p, 'Options',@isstruct)
addParameter(p, 'ksets',[],@isnumeric)
addParameter(p, 'n_rows',[],@isnumeric)
addParameter(p, 'n_cols',[],@isnumeric)
addParameter(p, 'show_legend',true,checkBool)
addParameter(p, 'metab_inds',[],@isnumeric)
addParameter(p, 'metab_names',default_metabs,checkNameOrCell)
addParameter(p, 'tspan',[],@isnumeric)
addParameter(p, 'conditions',[],@isnumeric)
addParameter(p, 'endtime',[],@isnumeric);
addParameter(p, 'met_inds_on_right_axis',[],@isnumeric)
addParameter(p, 'use_default_colors',0,checkBool)
addParameter(p, 'plot_exp_data',1,checkBool)
addParameter(p, 'plot_model_data',1,checkBool)
addParameter(p, 'y_log',0,checkBool)
addParameter(p, 'best_n_ksets',[],isNatural)
addParameter(p, 'use_true_ED',0,checkBool)
addParameter(p, 'plot_reference_conc',0,checkBool)
addParameter(p, 'Initial_Ensemble',[],@isstruct)
addParameter(p, 'show_fitness',1,checkBool)
addParameter(p, 'fitness_fxn',[],checkName)
addParameter(p, 'plot_per_species',1,checkBool)
addParameter(p, 'show_CI_and_best',0,checkBool)
addParameter(p, 'CI_int',.95, @isnumeric)
addParameter(p, 'show_all_true_ED', 0,checkBool)
addParameter(p, 'blank_plots',0,@isnumeric)
addParameter(p, 'show_true_fluxes', 0,checkBool)
addParameter(p, 'plot_per_condition', 0, checkBool)
addParameter(p, 'additional_metabs', [], checkNameOrEmpty)
addParameter(p, 'x_max', [], @isnumeric)
addParameter(p, 'show_grid', true, checkBool);
addParameter(p, 'subplot_titles', {}, @iscell)
addParameter(p, 'plot_default_ratios', true, checkBool)
addParameter(p, 'plot_ratios', {}, @iscell)
addParameter(p, 'show_axes_labels', true, checkBool)
addParameter(p, 'starting_figure', 1, @isnumeric)
addParameter(p, 'ksets_per_fig', [], @isnumeric)
addParameter(p, 'show_condition_in_legend', false, checkBool)
addParameter(p, 'color_inds', [], @isnumeric)
addParameter(p, 'exp_data_color', [], @(x) isnumeric(x) || ischar(x) )
addParameter(p, 'legend_condition_names', [], @isstring)
addParameter(p, 'plot_best_kset', true, checkBool)

p.KeepUnmatched = false;
p.CaseSensitive = false;

parse(p,model,Experimental_Data,Perturbation_Results,Options,varargin{:})

ksets = p.Results.ksets(:)';
n_rows = p.Results.n_rows;
n_cols = p.Results.n_cols;
n_blank_subplots = p.Results.blank_plots;
show_legend = p.Results.show_legend;
metab_inds = p.Results.metab_inds;
metab_names = p.Results.metab_names;
tspan_init = p.Results.tspan;
conditions = p.Results.conditions;
endtime = p.Results.endtime;
met_inds_on_right_axis = p.Results.met_inds_on_right_axis;
use_default_colors = p.Results.use_default_colors;
plot_exp_data = p.Results.plot_exp_data;
plot_model_data = p.Results.plot_model_data;
y_log = p.Results.y_log;
best_n_ksets = p.Results.best_n_ksets;
use_true_ED = p.Results.use_true_ED;
plot_reference_conc = p.Results.plot_reference_conc;
Initial_Ensemble = p.Results.Initial_Ensemble;
show_fitness = p.Results.show_fitness;
fitness_fxn = p.Results.fitness_fxn;
plot_per_species = p.Results.plot_per_species;
show_CI_and_best = p.Results.show_CI_and_best;
CI_int = p.Results.CI_int;
show_all_true_ED = p.Results.show_all_true_ED;
show_true_fluxes = p.Results.show_true_fluxes;
plot_per_condition = p.Results.plot_per_condition;
additional_metabs = p.Results.additional_metabs;
x_max = p.Results.x_max;
show_grid = p.Results.show_grid;
subplot_titles = p.Results.subplot_titles;
plot_default_ratios = p.Results.plot_default_ratios;
plot_ratios = p.Results.plot_ratios;
show_axes_labels = p.Results.show_axes_labels;
starting_figure = p.Results.starting_figure;
ksets_per_fig = p.Results.ksets_per_fig;
show_condition_in_legend = p.Results.show_condition_in_legend;
color_inds = p.Results.color_inds;
exp_data_color = p.Results.exp_data_color;
legend_condition_names = p.Results.legend_condition_names;
plot_best_kset = p.Results.plot_best_kset;


% TODO
% Check which figure (if plotting across many figures)
if ~isempty(ksets_per_fig)
    error("TODO")
    % Get index of this kset
    k_ind = find(ksets == k);
    fig_num = ceil(k_ind / ksets_per_fig) - 1;
    figure(starting_figure + fig_num)
end

% If plotting reference concs, make sure Initial_Ensemble is also given
if plot_reference_conc
    assert(~isempty(Initial_Ensemble),...
        "Must pass in 'Initial_Ensemble' varargin")
    all_ref_concs = [Initial_Ensemble.reference_metab_concs; ...
                    Initial_Ensemble.reference_enzyme_totals];
end

% Set true_ED if using
if use_true_ED
    assert(isfield(Experimental_Data, 'true_ED'),...
        "Experimental_Data does not have `true_ED` field")
    Experimental_Data = Experimental_Data.true_ED;
end

if isempty(conditions)
    conditions = 1:length(Experimental_Data.flux_rxns);
end

% Get fitness fxn in Options if empty
if isempty(fitness_fxn) 
    if ~isfield(Options,'fitness_fxn')
        Options.fitness_fxn = 'RMSE';
    end
    fitness_fxn = Options.fitness_fxn;
    
end

gray_rgb = [.6,.6,.6];
green_rgb = [0,.4,0];
yellow_rgb = [.7,.7,0];
black_rgb = [0,0,0];
red_rgb = [1 0 0];
blue_rgb = [0 0 1];

linewidth = 2;

line_styles = {'-','--',':','-.'};
n_lines = length(line_styles);
% Make the lines lighter
% what proportion towards white ([1 1 1]) should the color move? 
% % % Now uses exponential decay
% Now uses linear shift from 0 (pure color) to max_hue (1 would be white)
max_hue = 1;

if show_all_true_ED
    assert(isfield(Experimental_Data, 'true_ED'),...
        "Experimental_Data needs field `true_ED`")
    
    true_ED = Experimental_Data.true_ED;
end

if show_true_fluxes
    true_PR = model.true_fluxes.Perturbation_Results;
end

% Get new PR with new fitness function if changed)
if ~strcmpi(fitness_fxn, Options.fitness_fxn)
    Perturbation_Results = POST_recalc_kset_fitness(model,Experimental_Data,...
        Perturbation_Results,Options,'fitness_fxn',fitness_fxn);   
end

if ~isempty(best_n_ksets)
    
    [val_lowest_max, index_lowest_max, val_lowest_avg, index_lowest_avg]...
        edit POST= find_best_Ksets(Perturbation_Results, best_n_ksets);
        
    ksets = index_lowest_avg;
    ksets = ksets(:)';  
    
elseif isempty(ksets)
    
    ksets = 1:size(Perturbation_Results.fitness_values,1);
        
else
    
    ksets = ksets(:)';

end
avg_fitness_vals = Perturbation_Results.avg_fitness_values(ksets);
fitness_vals = Perturbation_Results.fitness_values(ksets,:);

if isempty(metab_inds) 
    
    if isequal(metab_names, default_metabs) && isempty(additional_metabs)
        subplot_titles = default_subplot_titles;
    end
    
    metab_inds = zeros(length(metab_names),1);
    for m=1:length(metab_names)
       
        met_ind = find(strcmpi(...
            model.metabs_and_enzyme_complexes, metab_names{m}));
        if length(met_ind) ~= 1
            
            if strcmpi(metab_names{m}, '1butanol')
                met_ind = find(strcmpi(...
                    model.metabs_and_enzyme_complexes, '1btoh_c'));
                assert(length(met_ind) == 1, "Fix this")
            else

                input_prompt = strcat("Could not find unique index of metabolite '",...
                    metab_names{m},"'.\nPlease input index of metabolite:\n");
                user_index = input(char(input_prompt),'s');
                met_ind = str2num(user_index);
            end
        end
        metab_inds(m) = met_ind;
        
    end
    
    if ~isempty(additional_metabs)
        if ~isnumeric(additional_metabs)
            [~, new_inds] = ismember(...
                additional_metabs, model.metabs_and_enzyme_complexes);
            if length(new_inds) ~= length(additional_metabs)
                error("Could not find all additional metabolites\n")
            end
        else
            new_inds = additional_metabs;
        end
        if any(new_inds == 0)
            error("One of the additional species is not in the model")
        end
        new_inds = new_inds(:);
        metab_inds = [metab_inds; new_inds(:)];
    end
end

% Get indices for metabs used in ratios
if ~isempty(plot_ratios)
    ratio_met_inds = cell(size(plot_ratios, 1), 1);
    for i = 1:size(plot_ratios, 1)
        [~, primary_ind] = ismember(plot_ratios{i, 1}, model.mets);
        [~, group_inds] = ismember(plot_ratios{i, 2}, model.mets);
        ratio_met_inds{i} = [primary_ind, setdiff(group_inds, primary_ind)];
    end
elseif plot_default_ratios
    ratio_met_inds = cell(size(default_ratios, 1), 1);
    for i = 1:size(default_ratios, 1)
        [~, primary_ind] = ismember(default_ratios{i, 1}, model.mets);
        [~, group_inds] = ismember(default_ratios{i, 2}, model.mets);
        
        % Make sure they're present
        if primary_ind == 0
            error("Fix me")
        end
        group_inds(group_inds == 0) = [];
        
        %Add to list
        ratio_met_inds{i} = [primary_ind, setdiff(group_inds, primary_ind)];
    end
else
    ratio_met_inds = {};
end
n_ratios = size(ratio_met_inds, 1);
n_metabs = length(metab_inds);
n_timecourses = n_ratios + n_metabs;
%% Get number of plots/rows/columns
n_conditions = length(conditions);
% If plotting by metabolite, each species is a subplot
if plot_per_species
    
    n_plots = length(metab_inds);
    if ~isempty(plot_ratios)
        n_plots = n_plots + size(plot_ratios, 1);
    elseif plot_default_ratios
        n_plots = n_plots + size(default_ratios, 1);
    end
    
else
    % If plotting by condition, each condition is a subplot
    n_plots = n_conditions;
    
    % Also make sure that if not plotting per species, that the confidence
    % interval isn't being done (plotting ksets individually) - would just
    % get messy and don't want to implement
    if show_CI_and_best
        fprintf("Turning off confidence interval view - not yet implemented\n")
        show_CI_and_best = 0;
    end
end

% Account for blank subplots
n_plots = n_plots + n_blank_subplots;

if isempty(n_rows) && isempty(n_cols)
    n_rows = floor(n_plots ^ 0.5);
    n_cols = ceil(n_plots / n_rows);
elseif ~isempty(n_rows)
    n_cols = ceil(n_plots / n_rows);
else
    n_rows = ceil(n_plots / n_cols);
end

%% Get axis limits
% Save as matrix - max across conditions and species
max_y = zeros(n_conditions, n_timecourses);
max_y_alt = zeros(n_conditions, n_timecourses);
min_y = zeros(n_conditions, n_timecourses);
min_y_alt = zeros(n_conditions, n_timecourses);

for c = 1:n_conditions
    
    for m = 1:length(metab_inds)
        
        for k = ksets
            
            % Get t-span needed for this kset & condition
            if ~isempty(endtime)
                t_end_index = find(Options.t_interval == endtime);
                assert(length(t_end_index) == 1,...
                    "Could not find time in t_interval matching endtime")
                tspan_inds = 1:t_end_index;

            elseif isempty(tspan_init)
%                 t_start_index = 1;
%                 t_end_index = size(Perturbation_Results.last_metab_concs{1,1},1);
%                 tspan = [t_start_index: t_end_index];
                tspan_inds = 1 : length(Perturbation_Results.timepoints{conditions(c), k});
            else
                assert(checkTspan(tspan_init),...
                    "tspan needs to be a 2-element integer array")
                tspan_inds = [tspan_init(1) : tspan_init(2)];
            end
            
            % Get max conc over this span
            max_kset_conc = max(Perturbation_Results.last_metab_concs{conditions(c),k}...
                (tspan_inds,metab_inds(m)));
            min_kset_conc = min(Perturbation_Results.last_metab_concs{conditions(c),k}...
                (tspan_inds,metab_inds(m)));
            
            if ~ismember(metab_inds(m),met_inds_on_right_axis)
                
                max_y(c,m) = max(max_y(c,m), max_kset_conc);
                min_y(c,m) = min(min_y(c,m), min_kset_conc);
            
            else
                max_y_alt(c,m) = max(max_y_alt(c,m), max_kset_conc);
                min_y_alt(c,m) = min(min_y_alt(c,m), min_kset_conc);
            end
            
        end
        
        % Add experimental data
        if plot_exp_data && isfield(Experimental_Data,'metab_timecourse') && ...
                ~isempty(Experimental_Data.metab_timecourse{conditions(c)})
            
            tc = Experimental_Data.metab_timecourse{conditions(c)};
            tc_met_ind = find(strcmpi(...
                tc(:,1), model.metabs_and_enzyme_complexes{metab_inds(m)}));

            if length(tc_met_ind) == 1

                % timecourse concentrations                 
                tc_concs = tc{tc_met_ind,2}(2,:);
                
                if ~ismember(metab_inds(m),met_inds_on_right_axis)
                    max_y(c,m) = max(max_y(c,m), max(tc_concs));
                    min_y(c,m) = min(min_y(c,m), min(tc_concs));
                else
                    max_y_alt(c,m) = max(max_y_alt(c,m),max(tc_concs));
                    min_y_alt(c,m) = min(min_y_alt(c,m),min(tc_concs));
                end
            end
            
            % Also look at true/underlying data if needed
            if show_all_true_ED
                
                tc_true = true_ED.metab_timecourse{conditions(c)};
                tc_true_met_ind = find(strcmpi(...
                    tc_true(:,1), model.metabs_and_enzyme_complexes{metab_inds(m)}));

                if length(tc_true_met_ind) == 1

                    % timecourse concentrations                 
                    tc_true_concs = tc_true{tc_true_met_ind,2}(2,:);

                    if ~ismember(metab_inds(m),met_inds_on_right_axis)
                        max_y(c,m) = max(max_y(c,m), max(tc_true_concs));
                        min_y(c,m) = min(min_y(c,m), min(tc_true_concs));
                    else
                        max_y_alt(c,m) = max(max_y_alt(c,m),max(tc_true_concs));
                        min_y_alt(c,m) = min(min_y_alt(c,m),min(tc_true_concs));
                    end
                end
                
            elseif show_true_fluxes
                
                true_flux_concs = true_PR.last_metab_concs...
                    {conditions(c), Options.true_fluxes_kset}...
                    (:, metab_inds(m));
                
                if ~ismember(metab_inds(m),met_inds_on_right_axis)
                    max_y(c,m) = max(max_y(c,m), max(true_flux_concs));
                    min_y(c,m) = min(min_y(c,m), min(true_flux_concs));
                else
                    max_y_alt(c,m) = max(max_y_alt(c,m), max(true_flux_concs));
                    min_y_alt(c,m) = min(min_y_alt(c,m), min(true_flux_concs));
                end
                
            end

        end
    end % Loop over metabs
    
    for r = 1:n_ratios
        primary_conc = Perturbation_Results.last_metab_concs{conditions(c),k}...
                (tspan_inds,ratio_met_inds{r}(1));
       
        all_ratio_concs = Perturbation_Results.last_metab_concs{conditions(c),k}...
                (tspan_inds,ratio_met_inds{r});
            
        rth_ratio = primary_conc ./ sum(all_ratio_concs, 2);
        
        max_y(c,m) = max(max_y(c,m), max(rth_ratio));
        min_y(c,m) = min(min_y(c,m), min(rth_ratio));    
    end
end

% Set colors
if use_default_colors && n_timecourses <= 6
    colors = [black_rgb; red_rgb; yellow_rgb; blue_rgb; gray_rgb; green_rgb];

% If plotting all conditions on same plot, use colors to distinguish
% between conditions
elseif ~plot_per_condition
    if ~isempty(color_inds)
        assert(max(color_inds) >= n_conditions)
        all_colors = linspecer(max(color_inds), 'qualitative');
        colors = all_colors(color_inds, :);
    else
        colors = linspecer(n_conditions, 'qualitative');
    end
    
% Otherwise use colors to distinguish between species
else
    if ~isempty(color_inds)
        assert(max(color_inds) >= n_timecourses)
        all_colors = linspecer(n_timecourses, 'qualitative');
        colors = all_colors(color_inds, :);
    else
        colors = linspecer(n_timecourses, 'qualitative');
    end
end

handles_for_met_legend = [];
handles_for_kset_legend = [];
met_legend_labels = {};
kset_legend_labels = {};
legend_handles = [];
legend_labels = {};

figure(starting_figure)
clf

for c = 1:n_conditions
    
    if ~plot_per_species
        subplot(n_rows, n_cols, c + n_blank_subplots);
        hold on
        if show_grid
            grid on
        end
    end
    
    for m = 1:n_timecourses
        
        if plot_per_species
            
            % If making a plot/figure for each condition, 
            if plot_per_condition
                figure(starting_figure - 1 + c)
                if m == 1
                    clf
                end
            end
            % Specify the subplot for each metabolite
            subplot(n_rows, n_cols, m + n_blank_subplots);
            hold on
            if show_grid
                grid on
            end

        elseif plot_per_condition
            error("TODO: implement this")           
        
        end
        
        % counter for line styles
        ls = 1;
        hue = 0;
        shift_counter = 0;
        
        if m <= n_metabs && ~isempty(met_inds_on_right_axis)
            if ismember(metab_inds(m),met_inds_on_right_axis)

                yyaxis right
                % If plotting by metabolite, use the max for this metab
                % across all conditions
                if plot_per_species
                    ylim([min(0, min(min_y_alt(:,m))), max(max_y_alt(:,m))])
                else
                    % Otherwise get max for this condition across all
                    % metabs
                    ylim([min(0, min(min_y_alt(c,:))), max(max_y_alt(c,:))])
                end

            else
                yyaxis left
                if plot_per_species
                    ylim([min(0, min(min_y(:,m))), max(max_y(:,m))])
                else
                    % Otherwise get max for this condition across all
                    % metabs
                    ylim([min(0, min(min_y(c,:))), max(max_y(c,:))])
                end
            end
        else
            if plot_per_species
               
                y_limits = [min(0, min(min_y(:,m))), max(max_y(:,m))];
                if y_limits(2) <= y_limits(1)
                    y_limits(2) = y_limits(1) + 1;
                end
                ylim(y_limits)
            else
                % Otherwise get max for this condition across all
                % metabs
                y_limits = [min(0, min(min_y(c,:))), max(max_y(c,:))];
                if y_limits(2) <= y_limits(1)
                    y_limits(2) = y_limits(1) + 1;
                end
                ylim(y_limits)
            end
        end
        
        % If plotting confidence interval, need to hold all ksets info
        all_kset_tc = [];
        
        for k = ksets
            
            % Get index of k in this list of ksets
            kset_ind = find(ksets == k);
            
            % If plotting each condition separately 
            % (using color for metabolites), make lighter each time
            if plot_per_condition
                % So each time a new line style is needed, the line gets a bit
                % lighter
                basecolor = colors(m,:);
            
            % Otherwise, use color for this condition
            else
                basecolor = colors(c,:);
            end
            % Each time a line style is repeated, the line gets a bit
            % lighter
            linecolor = basecolor + ([1,1,1] - basecolor) * hue;
            
            % Plot model-predicted timecourses
            if plot_model_data
                
                % Get t-span needed for this kset & condition
                if ~isempty(endtime)
                    t_end_index = find(Options.t_interval == endtime);
                    assert(length(t_end_index) == 1,...
                        "Could not find time in t_interval matching endtime")
                    tspan_inds = 1:t_end_index;

                elseif isempty(tspan_init)
%                     t_start_index = 1;
%                     t_end_index = size(Perturbation_Results.last_metab_concs{1,1},1);
%                     tspan = [t_start_index: t_end_index];
                    tspan_inds = 1 : length(Perturbation_Results.timepoints{conditions(c), k});
                else
                    assert(checkTspan(tspan_init),...
                        "tspan needs to be a 2-element integer array")
                    tspan_inds = [tspan_init(1) : tspan_init(2)];
                end
                
                % Turn off axes interpreter
                set(0,'DefaultLegendInterpreter','none')
                    
                % Get x values
                tc_x_vals = Perturbation_Results.timepoints...
                    {conditions(c),k}(tspan_inds);
                
                
                % Get this metab's (or ratios) timecourses
                if m <= n_metabs
                    kset_tc = Perturbation_Results.last_metab_concs...
                            {conditions(c),k}(tspan_inds,metab_inds(m));
                    tc_name = model.metabs_and_enzyme_complexes{metab_inds(m)};
                    
                    % Check if this has a "name" saved
                    [~,def_names_ind] = ismember(tc_name, default_metab_names);
                    if def_names_ind ~= 0
                        tc_name = default_metab_names(def_names_ind, 2);
                    end
                    
                else
                    r = m - n_metabs;
                    primary_conc = Perturbation_Results.last_metab_concs{conditions(c),k}...
                        (tspan_inds,ratio_met_inds{r}(1));
       
                    all_ratio_concs = Perturbation_Results.last_metab_concs{conditions(c),k}...
                            (tspan_inds,ratio_met_inds{r});
            
                    kset_tc = primary_conc ./ sum(all_ratio_concs, 2);
                    
                    
                    primary_name = model.mets{ratio_met_inds{r}(1)};
                    [~,alias_row] = ismember(primary_name, default_metab_names);
                    if alias_row ~= 0
                        primary_name = default_metab_names(alias_row, 2);
                    end
                    tc_name = strcat(primary_name, " ratio");
                end
                
                % If plotting confidence interval and best kset, get all
                % data    
                if show_CI_and_best
                    
                    all_kset_tc = [all_kset_tc; kset_tc(:)'];
                    
                                        
                % Otherwise, plot each metab individually
                else
                    h = plot(tc_x_vals, kset_tc,...
                        'LineWidth',linewidth,'Color',linecolor,...
                        'LineStyle',line_styles{ls},...
                        'Marker','none');

                    ax = h.Parent;
                
                    if y_log
                        set(gca,'YScale','log')
                    end

                    if c==1 && k==ksets(1) && show_legend
    % % %                     set(h,'DisplayName',...
    % % %                         model.metabs_and_enzyme_complexes{metab_inds(m)});
    % % %                     leg1 = legend(h);

                        handles_for_met_legend = [handles_for_met_legend, h];
                        met_legend_labels = [met_legend_labels, ...
                            char(tc_name)];

                    else
                        % Remove this entry from legend
    % % %                     set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                    end

                    % Make second legend for each kset (only first metabolite)
                    if m == 1

                        

                        h_kset = copyobj(h, ax);

                        handles_for_kset_legend = [handles_for_kset_legend, h_kset];
    %                     kset_legend_labels = [kset_legend_labels, char(...
    %                         strcat("K-set ",num2str(k)," - ",fitness_fxn," = ",...
    %                         num2str(fitness_vals(kset_ind) ) ) )];
                        legend_label = char(...
                            sprintf("k-set %i - %s = %0.2f",...
                                k, fitness_fxn, ...
                                fitness_vals(kset_ind, conditions(c))));
                        
                            
                        if ~isempty(legend_condition_names)
                            legend_label = sprintf('%s - %s',...
                                legend_condition_names(c),...
                                legend_label);
                            
                        elseif show_condition_in_legend && ...
                                isfield(Experimental_Data, 'condition_names')
                            legend_label = sprintf('%s - %s',...
                                Experimental_Data.condition_names{conditions(c)},...
                                legend_label);
                        end
                        
                        kset_legend_labels = [kset_legend_labels, legend_label];

                    end
                    
                end  %End condition for plotting as CI vs individual ksets
                
            end
            
            % Plot reference concs & Only plot if a metabolite or free enzyme
            if plot_reference_conc && metab_inds(m) <= length(all_ref_concs)
               
                % Get ref conc for this metab and kset
                ref_metab_conc = all_ref_concs(metab_inds(m), k);
                
                assert(isequal(size(ref_metab_conc), [1,1]), "Fix this")
                
                % Horizontal line for each metab being plotted
                y = yline(ref_metab_conc,...
                    'Color', linecolor,...
                    'LineWidth', linewidth,...
                    'LineStyle', line_styles{ls});
                
                set(y,'HandleVisibility','off');
                
            end

            ls = ls + 1;
            if ls > n_lines
                % reset line style counter
                ls = 1;
                
                % How many times does the color need to shift
                n_shifts = ceil(length(ksets) / n_lines);
                
                % Index hue shift counter                
                shift_counter = shift_counter + 1;
                
                % Get new hue
                % Hue starts at 0 (pure basecolor) and moves to 0.5 (1 is white)
                
                hue = max_hue * shift_counter / n_shifts;
                
            end

            if plot_per_species
                if ~isempty(subplot_titles) && m <= length(subplot_titles)
                    title(subplot_titles{m}, 'Interpreter','none')
                elseif exist('tc_name','var')
                    title(tc_name, 'Interpreter','none')
                end
            else
                title(strcat("Condition ",num2str(conditions(c))), 'Interpreter','none')
            end
            
        end % End loop over ksets
        
        % If plotting CI (not individual ksets), do that now
        if show_CI_and_best && plot_model_data
            
            if ~isempty(legend_condition_names)
                condition_name = legend_condition_names(c);

            elseif show_condition_in_legend && ...
                    isfield(Experimental_Data, 'condition_names')
                
                condition_name = Experimental_Data.condition_names{conditions(c)};
            else
                condition_name = '';
            end
            
            plot_timecourse_CI_intervals(tc_x_vals, all_kset_tc, ksets, colors(c, :),...
                'CI_int', CI_int,...
                'plot_best_kset', plot_best_kset,...
                'condition_name', condition_name);
            
% % %             % Just use first kset given (best or not) as line
% % %             y_best = all_kset_tc(1,:);
% % %             
% % %             % Allow multiple CI intervals - order from largest to smallest
% % %             % and then plot from lightest to darkest
% % %             assert(all(CI_int < 1) && all(CI_int > 0),...
% % %                 "CI interval values must be between 0 and 1")
% % %             sorted_CI_ints = sort(CI_int,'descend');
% % %             
% % %             if plot_per_condition
% % %                 % If having each condition on a different plot, use gray
% % %                 basecolor = [.9, .9, .9];
% % %                 best_k_color = 'k'; % black
% % %             else
% % %                 % Otherwise, use color for that condition
% % %                 basecolor = colors(c, :);
% % %                 best_k_color = colors(c, :);
% % %             end
% % %             
% % %             int_colors = repmat(basecolor, length(CI_int), 1);
% % %             % Get shades multipleirs ([0.95, 0.9, 0.85, ...])
% % %             shades_multiplier = 1 - 0.2 .* [1 : length(CI_int)]';
% % %             int_colors = int_colors .* shades_multiplier;
% % %                 
% % %             % Plot each confidence interval
% % %             for int = 1:length(sorted_CI_ints)
% % %                 CI = CI_matrix(all_kset_tc, sorted_CI_ints(int));
% % %                 y_upper_CI = CI(2,:);
% % %                 y_lower_CI = CI(1,:);
% % %                 % Don't let lower go below 0
% % %                 y_lower_CI = max([zeros(size(y_lower_CI)); y_lower_CI]);
% % %                 CI_color = int_colors(int, :);
% % % 
% % %                 CI_legend = char(strcat(num2str(sorted_CI_ints(int)*100),...
% % %                         "% CI"));
% % %                     
% % %                 if ~isempty(legend_condition_names)
% % %                     CI_legend = sprintf('%s - %s',...
% % %                         legend_condition_names(c),...
% % %                         CI_legend);
% % %                             
% % %                 elseif show_condition_in_legend && ...
% % %                         isfield(Experimental_Data, 'condition_names')
% % %                     CI_legend = sprintf('%s - %s',...
% % %                         Experimental_Data.condition_names{conditions(c)},...
% % %                         CI_legend);
% % %                 end
% % %                 
% % %                 h_CI = fill([tc_x_vals(:)', fliplr(tc_x_vals(:)')],...
% % %                     [y_upper_CI, fliplr(y_lower_CI)],...
% % %                     CI_color,...
% % %                     'EdgeAlpha',0,...
% % %                     'DisplayName', CI_legend,...
% % %                     'FaceAlpha', 0.3);
% % %                 
% % %             end
% % % 
% % %             h = plot(tc_x_vals, y_best,...
% % %                 'Color', best_k_color,...
% % %                 'LineWidth', 2,...
% % %                 'DisplayName', kset_legend);

            % Plot line for best kset
            kset_legend = char(sprintf("k-set %i - %s = %0.2f",...
                    ksets(1), fitness_fxn, fitness_vals(1, conditions(c))));

            
        end
        
        
        % Add experimental data
        if plot_exp_data && isfield(Experimental_Data,'metab_timecourse') && ...
                ~isempty(Experimental_Data.metab_timecourse{conditions(c)})
            
            % Set colors for points based on plot type
            % If each condition has its own figure, use metabolite colors
            % for "true" data, black for training data
            if plot_per_condition
                truedata_color = colors(m, :);
                training_color = 'k'; % black
            
            % If plotting all conditions on single figure, use lighter
            % color for "true" data if required, darker version of same
            % color for training data
            elseif ~isempty(exp_data_color)
                truedata_color = exp_data_color;
                training_color = exp_data_color;
                
            else
                training_color = colors(c, :);
                truedata_color = training_color + ...
                    ([1,1,1] - training_color) .* 0.5;
            end                
            
            % If option to show all "underlying"/true data and this metab has no
            % synthetic data, plot colored point first
            if show_all_true_ED && m <= n_metabs
                tc_true = true_ED.metab_timecourse{conditions(c)};
                tc_true_met_ind = find(strcmpi(...
                    tc_true(:,1), model.metabs_and_enzyme_complexes{metab_inds(m)}));
                
                if length(tc_true_met_ind) == 1

                    % timecourse times and concentrations
                    tc_true_times = tc_true{tc_true_met_ind,2}(1,:);                    
                    tc_true_concs = tc_true{tc_true_met_ind,2}(2,:);

                    % Set marker size and line width
                    mrk_sz = 30;
                    mrk_ln = 1.0;

                    truedata_label = "True Data";
                    if ~plot_per_condition
                        truedata_label = truedata_label + ...
                            ": Condition #" + conditions(c);
                    end
                    
                    scatter(tc_true_times, tc_true_concs,...
                        mrk_sz, truedata_color,...
                        'LineWidth', mrk_ln,...
                        'DisplayName', truedata_label)
                    
                    hold on
                    
                end

            elseif show_true_fluxes && m <= n_metabs
                true_flux_concs = true_PR.last_metab_concs...
                    {conditions(c), Options.true_fluxes_kset}...
                    (:, metab_inds(m));
                true_flux_times = true_PR.timepoints...
                    {conditions(c), Options.true_fluxes_kset};
                
                % Set marker size and line width
                mrk_sz = 30;
                mrk_ln = 1.0;
                
                truedata_label = "True Data";
                if ~plot_per_condition
                    truedata_label = truedata_label + ...
                        ": Condition #" + conditions(c);
                end

                scatter(true_flux_times, true_flux_concs,...
                    mrk_sz, truedata_color,...
                    'LineWidth', mrk_ln,...
                    'DisplayName', truedata_label)

                hold on                
                
            end
            
            % Then plot darker/black experimental (training) data
            if m <= n_metabs
                tc = Experimental_Data.metab_timecourse{conditions(c)};
                tc_met_ind = find(strcmpi(...
                    tc(:,1), model.metabs_and_enzyme_complexes{metab_inds(m)}));

                if length(tc_met_ind) == 1

                    % timecourse times and concentrations
                    tc_times = tc{tc_met_ind,2}(1,:);                    
                    tc_concs = tc{tc_met_ind,2}(2,:);

                    % Set marker size and line width
                    mrk_sz = 30;
                    mrk_ln = 1.0;

                    traindata_label = "Training Data";
                    if ~plot_per_condition
                        traindata_label = traindata_label + ...
                            ": Condition #" + conditions(c);
                    end

                    scatter(tc_times, tc_concs,...
                        mrk_sz, training_color,...
                        'LineWidth',mrk_ln,...
                        'DisplayName', traindata_label)

                end
            end
            
        end
        
        if show_axes_labels
            xlabel('Time (hrs)')
            ylabel('Conc (mM)')
        end
        
        if ~isempty(x_max)
            xlim([0, x_max])
        end
        
    end % Loop over metabolites
    
    
    
% % %     if show_legend && ~plot_per_species
% % %        
% % %         legend_handles = handles_for_met_legend;
% % %         legend_labels = [met_legend_labels];
% % %         
% % %     end
% % %     
% % %     % Option to add table for kset fitnesses
% % %     if show_fitness
% % %         
% % %         legend_handles = [legend_handles, handles_for_kset_legend];
% % %         legend_labels = [legend_labels, kset_legend_labels];
% % %         
% % %     end
% % %     
% % %     if show_legend || show_fitness
% % %         
% % %         lgd = legend(legend_handles, legend_labels);
% % %         lgd.AutoUpdate = 'off';
% % %         % Don't know why this is needed
% % %         lgd.String = legend_labels;
% % %         
% % %     end
        
    
%     ylim([0,max_y])

    
end

% Moved these outside of loop - 2020-10-23
if show_legend && ~plot_per_species
       
    legend_handles = handles_for_met_legend;
    legend_labels = [met_legend_labels];

end

% Option to add table for kset fitnesses
if show_fitness

    legend_handles = [legend_handles, handles_for_kset_legend];
    legend_labels = [legend_labels, kset_legend_labels];

end

if show_legend || show_fitness

    lgd = legend(legend_handles, legend_labels);
    lgd.AutoUpdate = 'off';
    % Don't know why this is needed
    lgd.String = legend_labels;

end

if ~show_legend
    legend('off')
end


end

% Moved to separate file

% % % function CI = CI_vector(x, p)
% % % 
% % % % Taken from Adam Danz - 
% % % % https://www.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
% % % 
% % % CI_alt = @(x,p) std(x(:),'omitnan')/sqrt(sum(~isnan(x(:)))) * ...
% % %     tinv(abs([0,1]-(1-p)/2),sum(~isnan(x(:)))-1) + mean(x(:),'omitnan'); 
% % % 
% % % % Function to get percentil (z-test) CI
% % % CI_z = @(x,p)prctile(x,abs([0,100]-(100-p*100)/2));
% % % 
% % % % Now taken from https://www.mathworks.com/help/stats/tinv.html
% % % 
% % % % Compute the sample mean, standard error, and degrees of freedom.
% % % xbar = mean(x);
% % % n = sum(~isnan(x(:)));
% % % se = std(x)/sqrt(n);
% % % nu = n - 1;
% % % % Find the upper and lower confidence bounds for the 95% confidence interval.
% % % alpha = 1 - p;
% % % alphaLo = alpha/2;
% % % alphaHi = 1 - alpha/2;
% % % % Compute the critical values for the confidence bounds.
% % % crit = tinv([alphaLo alphaHi], nu);
% % % % Determine the confidence interval for the population mean.
% % % CI = xbar + crit*se;
% % % 
% % % % Instead use bootstrapping - doesn't assume normal distribution
% % % % [CI_mean_std, bootstat] = bootci(1e2, {@(x)[mean(x) std(x)], x}, ...
% % % %     'Alpha', alpha,...
% % % %     'Type', 'bca',...
% % % %     'NBootStd',10);
% % % % 
% % % % CI = CI_mean_std(:, 1);
% % % 
% % % % Can also just do raw percentiles
% % % 
% % % CI = prctile(x, (100.*(1 - [alphaHi alphaLo])) );
% % % 
% % % end
% % % 
% % % function CI = CI_matrix(M, p, dim)
% % % % Calculate CI for each dim of M - 
% % % % dim==1 is CI for each column (default), dim==2 is CI for each row - 
% % % 
% % % if ~exist('dim','var')
% % %     dim = 1;
% % % end
% % % 
% % % if dim == 2
% % %     M = M';
% % % end
% % % 
% % % n_vectors = size(M,2);
% % % CI = zeros(2,n_vectors);
% % % 
% % % for i=1:n_vectors
% % %     single_CI = CI_vector(M(:,i), p);
% % %     CI(:,i) = single_CI(:);
% % % end
% % % 
% % % if dim == 2
% % %     CI = CI';
% % % end
% % % end

