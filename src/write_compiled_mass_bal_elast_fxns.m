function write_compiled_mass_bal_elast_fxns(model, Experimental_Data,...
    Options, kset_and_init_concs)

%% Set up Initial_Ensemble
Initial_Ensemble = parse_param_vector(model, Experimental_Data, Options,...
    kset_and_init_concs);

sat_constants = Initial_Ensemble.sat_constants;
S_f_b_trimmed = model.S_f_b;
S_trimmed = model.S;

%% Trim based on constant species
% Get fields for fixing metab timecourses across timepoints (i.e., pH)
if Options.fix_experimental_timecourses && ...
        ~isempty(Experimental_Data.fixed_metab_timecourse{1})
    fixed_met_inds = Experimental_Data.fixed_metab_timecourse{1}(:, 1);
    if iscell(fixed_met_inds)
        [~, fixed_met_inds] = ismember(fixed_met_inds, model.mets);
    end
    
else
    fixed_met_inds = [];
end
const_met_inds = [model.constant_metabs_list(:); fixed_met_inds(:)];

% Add in from Options (varargin) if there
if isfield(Options, 'const_metabs') && ~isempty(Options.const_metabs)
    if isnumeric(Options.const_metabs)
        vararg_const_met_inds = Options.const_metabs(:, 1);
        vararg_const_met_levels = Options.const_metabs(:, 2);
    else
        [~,vararg_const_met_inds] = ...
            ismember(string(Options.const_metabs(:,1)), model.mets);
        assert(~any(vararg_const_met_inds == 0),...
            "Invalid input for 'const_metabs'")
        vararg_const_met_levels = cell2mat(...
            Options.const_metabs(:, 2));
    end
    const_met_inds = [const_met_inds; vararg_const_met_inds(:,1)];
end

sat_constants(:, const_met_inds) = [];
S_f_b_trimmed(const_met_inds, :) = [];
S_trimmed(const_met_inds, :) = [];

% Do the same for model.lit_values.sat_cosntants_matrix_locs
sat_const_matrix_locs = model.lit_values.sat_constants_matrix_locs;
sat_const_matrix_trimmed = sat_const_matrix_locs;
% Sort descending so we can just lower any below in order and not miss some
const_met_inds = sort(const_met_inds, 'descend');
for cm_ind = 1:length(const_met_inds)
    % Remove any rows that use that
    rm_row = find(sat_const_matrix_locs(:, 2) == const_met_inds(cm_ind));
    sat_const_matrix_trimmed(rm_row, :) = 0; % Don't actually remove - just set to zero
    % Get rows whose metabolites had higher index than one removed
    higher_ind_rows = find(sat_const_matrix_trimmed(:,2) > const_met_inds(cm_ind));
    % Lower those indices by 1
    sat_const_matrix_trimmed(higher_ind_rows, 2) = ...
        sat_const_matrix_trimmed(higher_ind_rows, 2) - 1;
end

%% Save back into model struct
model.lit_values.sat_constants = sat_constants;
model.lit_values.sat_constants_matrix_locs_trimmed = sat_const_matrix_trimmed;

model.S_f_b_trimmed = S_f_b_trimmed;
model.S_trimmed = S_trimmed;

fprintf("WARNING: this will be incorrect if using Conservation Analysis\n")
% % % model.conserved_model_info.S_f_b_conserved = S_f_b_trimmed; % This was dumb
% Want this to throw an error if I try to use it in the next script
model.conserved_model_info.S_f_b_conserved = NaN; 

%% TODO - check somehow if already here


%% Write functions
[dv_dx_str, v_str, dv_dk_str] = zz_create_elasticity_coefficient_matrix(model, Options);

zz_write_elasticity_coefficient_function(model, dv_dx_str);

zz_write_mass_balance_function(model, v_str);

zz_write_elasticity_param_coeff_function(model, dv_dk_str);

end